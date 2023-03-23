using Random
using LinearAlgebra
P3 = Point3{Float64}

function check_keys(keys::Matrix{T}, returned_keys::Matrix{T}) where {T}
    s = Set{Point3{T}}(Point3{T}.(eachcol(keys)))
    @test all(k in s for k in eachcol(returned_keys))
end

"""returns: max dot products, query time"""
function best_answers_and_query_time(ds::AbstractInnerProductMax{T}, keys::Matrix{T}, queries::Matrix{T}) where {T}
    query_time = @timed(returned_keys = best_vecs(ds, queries)).time
    check_keys(keys, returned_keys)
    dot.(eachcol(returned_keys), eachcol(queries)), query_time
end

"""returns: max dot products, preprocess time, query time"""
function best_answers_and_times(t::DataType, hull::Hull{T}, keys::Matrix{T}, queries::Matrix{T}) where {T}
    preprocess_time = @timed(ds = t(hull)).time
    answers, query_time = best_answers_and_query_time(ds, keys, queries)
    answers, preprocess_time, query_time
end

function test_point_set(keys::Matrix{T}, values::Matrix{T}, t1::DataType, t2::DataType) where {T<:Real}
    hull = Hull{T}(keys)
    a1, p1, q1 = best_answers_and_times(t1, hull, keys, values)
    a2, p2, q2 = best_answers_and_times(t2, hull, keys, values)
    @test a1 ≈ a2
    (p1, q1), (p2, q2)
end

"""gen: takes as input an integer and outputs some number of points"""
function test_point_set_gen(gen, n::Int, q::Int, t1::DataType, t2::DataType)
    keys = gen(n)
    values = gen(q)
    test_point_set(keys, values, t1, t2)
end

@testset "InnerProductMax" begin
    """keys: 3xn"""

    MAX_N = 10

    function test_tetrahedron(t1::DataType, t2::DataType)
        keys = tetrahedron()
        values = cube(10, 100)
        test_point_set(keys, values, t1, t2)
    end

    function test_int(t1::DataType, t2::DataType)
        for d in 1:100
            for n in 4:MAX_N
                test_point_set_gen(n -> cube(d, n), n, 100, t1, t2)
            end
        end
    end

    function test_sphere(t1::DataType, t2::DataType)
        for n in 4:MAX_N
            for trial in 1:100
                test_point_set_gen(unit_sphere, n, 100, t1, t2)
            end
        end
    end

    function test_all_small(t1::DataType, t2::DataType)
        Random.seed!(1234) # for consistency
        test_tetrahedron(t1, t2)
        test_int(t1, t2)
        test_sphere(t1, t2)
    end

    function sanity_check(t::DataType)
        keys = tetrahedron()
        hull = Hull{Float64}(keys)
        # for simplex in eachcol(hull.simplices)
        #     println(hull.points[:, simplex[1]])
        #     println(hull.points[:, simplex[2]])
        #     println(hull.points[:, simplex[3]])
        #     println()
        # end
        ds = t(hull)
        if t == InnerProductMaxNaive
            @test ds.points == P3[[0.0, 0.0, 0.0], [0.0, 1.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 0.0]]
        end
        # @testset "Upward" begin
        @test query(ds, P3(1, 2, 3)) == P3(0, 1, 1)
        @test query(ds, P3(3, 1, 2)) == P3(1, 0, 1)
        @test query(ds, P3(0, 0, 1)) in [P3(1, 0, 1), P3(0, 1, 1)]
        # end
        # @testset "Downward" begin
        @test query(ds, P3(-1, -1, -1)) == P3(0, 0, 0)
        @test query(ds, P3(-1, -2, -3)) == P3(0, 0, 0)
        @test query(ds, P3(3, -2, -3)) == P3(1, 1, 0)
        @test query(ds, P3(1, -2, -3)) == P3(0, 0, 0)
        # end
        # @testset "Best Answers" begin
        @test best_answers_and_query_time(ds, keys, vec_to_matrix([
            P3(1, 2, 3),
            P3(-1, -2, -3),
            P3(4, 1, 2),
        ]))[1] == [5, 0, 6]
        # end
    end

    @testset "sanity_check_naive" begin
        sanity_check(InnerProductMaxNaive{Float64})
    end

    @testset "sanity_check_mine" begin
        sanity_check(InnerProductMaxMine{Float64})
    end

    @testset "sanity_check_classic" begin
        sanity_check(InnerProductMaxClassic{Float64})
    end

    @testset "small_naive" begin
        test_all_small(InnerProductMaxNaive{Float64}, InnerProductMaxNaive{Float64})
    end

    @testset "small_mine" begin
        test_all_small(InnerProductMaxNaive{Float64}, InnerProductMaxMine{Float64})
    end

    @testset "classic" begin
        test_all_small(InnerProductMaxNaive{Float64}, InnerProductMaxClassic{Float64})
    end
end