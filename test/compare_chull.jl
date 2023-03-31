# https://juliapolyhedra.github.io/
using CDDLib, QHull
using Polyhedra: Polyhedra, polyhedron, vrep
using InnerProductMax: HullLite

const libs = [Polyhedra.default_library(3, T), CDDLib.Library(:float), QHull.Library()]

function get_volume(point_set::AbstractMatrix{T}, lib::Polyhedra.Library, remove_duplicates::Bool=true) where {T}
    if lib != QHull.Library() && remove_duplicates
        # see https://github.com/JuliaPolyhedra/Polyhedra.jl/issues/321
        point_set = unique(point_set, dims=1)
    end
    poly = polyhedron(vrep(point_set), lib)
    Polyhedra.volume(poly)
end

@testset "volumes" begin
    @testset "test_volume_tricky1" begin
        points = [-1.0 0.0 1.0
            0.0 0.0 1.0
            0.0 0.0 0.0
            -1.0 -1.0 0.0
            -1.0 -1.0 0.0]
        for lib in libs
            @test get_volume(points, lib) ≈ 1.0 / 6
        end
        @test volume(HullLite{T}(Point3{T}.(eachrow(points)))) ≈ 1.0 / 6
    end

    @testset "test_volume_tricky2" begin
        points = [-1.0 0.0 1.0
            0.0 0.0 1.0
            0.0 0.0 0.0
            -1.0 -1.0 0.0
            -1.0 -1.0 0.0]
        for lib in libs # note: will fail if remove_duplicates = false
            @test get_volume(points, lib, true) ≈ 1.0 / 6
        end
        @test volume(HullLite{T}(Point3{T}.(eachrow(points)))) ≈ 1.0 / 6
    end

    @test_skip @testset "test_volume_tricky3" begin
        points = [-1.0 1.0 1.0
            -1.0 -1.0 1.0
            1.0 0.0 -1.0
            1.0 1.0 1.0
            0.0 1.0 -1.0]
        for lib in libs # default -> 5.333333333333333, qhull -> 2.6666666666666665
            vol = get_volume(points, lib, false)
            @test vol ≈ 8.0 / 3
            if !(vol ≈ 8.0 / 3)
                println("failed $lib got $vol")
            end
        end
    end

    @test_skip @testset "test_volume_tricky3_noised" begin
        points = [-1.0 1.0 1.0
            -1.0 -1.0 1.0
            1.0 0.0 -1.0
            1.0 1.0 1.0
            0.0 1.0 -1.0]
        Random.seed!(1234)
        points = points .+ randn(T, (5, 3)) * 0.001
        println(points)
        for lib in libs # default -> 5.333333333333333, qhull -> 2.6666666666666665
            vol = get_volume(points, lib, false)
            println(vol)
            # @test vol ≈ 8.0 / 3
            # if !(vol ≈ 8.0 / 3)
            #     println("failed $lib got $vol")
            # end
        end
    end

    @testset "test_volume" begin
        points = [
            0.0 0.0 0.0
            1.0 0.0 0.0
            0.0 1.0 0.0
            0.0 0.0 1.0
        ]
        for lib in libs
            @test get_volume(points, lib) ≈ 1.0 / 6
        end
        @test volume(HullLite{T}(Point3{T}.(eachrow(points)))) ≈ 1.0 / 6
    end

    function get_volumes(point_sets::Vector{AbstractMatrix{T}}, lib::Polyhedra.Library) where {T}
        vols = []
        for (i, p) in enumerate(point_sets)
            push!(vols, get_volume(p, lib))
        end
        vols
    end

    function get_volumes_lite(point_sets::Vector{AbstractMatrix{T}}) where {T}
        vols = T[]
        for (i, p) in enumerate(point_sets)
            push!(vols, volume(HullLite{T}(Point3{T}.(eachrow(p)))))
        end
        vols
    end

    function compare(libs::Vector, gen)
        point_sets = AbstractMatrix{T}[transpose(gen()) for _ in 1:10000]
        @time volumes = get_volumes_lite(point_sets)
        for lib in libs
            # println(lib)
            @time lib_volumes = get_volumes(point_sets, lib)
            @test lib_volumes ≈ volumes
            if !(lib_volumes ≈ volumes)
                dif = abs.(lib_volumes - volumes)
                num_differences = sum(.!(dif .≈ 0))
                i = argmax(dif)
                display(point_sets[i])
                println("Failed $lib, i = $i $(dif[i]) $(volumes[i]) $(lib_volumes[i]) $(num_differences)")
            end
        end
    end

    @testset "compare_volume_1_5" begin
        Random.seed!(1234)
        compare(libs[end:end], () -> cube(T, 1, 5))
    end

    @test_skip @testset "compare_volume_10_30" begin
        Random.seed!(1234)
        compare(libs[end:end], () -> cube(T, 10, 30))
    end
end