using Random
using LinearAlgebra
using Profile
using ProfileView

function check_keys(keys::Matrix{T}, returned_keys::Matrix{T}) where {T}
    s = Set{Point3{T}}(Point3{T}.(eachcol(keys)))
    @test all(k in s for k in eachcol(returned_keys))
end

"""returns: max dot products, query time"""
function best_answers_and_query_info(ds::AbstractInnerProductMax{T}, keys::Matrix{T}, queries::Matrix{T}) where {T}
    query_info = @timed(begin
        returned_keys = query_all(ds, queries)
        0
    end)
    check_keys(keys, returned_keys)
    dot.(eachcol(returned_keys), eachcol(queries)), query_info
end

"""returns: max dot products, preprocess time, query time"""
function best_answers_and_infos(t::DataType, hull::Hull{T}, keys::Matrix{T}, queries::Matrix{T}, profile::Bool=false) where {T}
    # @profview ds = t(hull)
    # @assert false
    if t == InnerProductMax.InnerProductMaxNested{T} && profile
        @profview preprocess_info = @timed(begin
            ds = t(hull)
            0
        end)
    else
        preprocess_info = @timed(begin
            ds = t(hull)
            0
        end)
    end
    # if t == InnerProductMax.InnerProductMaxNested{T} && profile
    #     # println("best_answers ", t)
    #     answers, query_info = best_answers_and_query_info(ds, keys, queries)
    #     # Profile.print()
    #     answers, preprocess_info, query_info
    # else
    answers, query_info = best_answers_and_query_info(ds, keys, queries)
    answers, preprocess_info, query_info
    # end
end

function test_point_set(keys::Matrix{T}, queries::Matrix{T}, t1::DataType, t2::DataType) where {T<:Real}
    hull = Hull(keys)
    a1, p1, q1 = best_answers_and_infos(t1, hull, keys, queries)
    a2, p2, q2 = best_answers_and_infos(t2, hull, keys, queries)
    @test a1 ≈ a2
    if !(a1 ≈ a2)
        dif = abs.(a1 - a2)
        amax = argmax(dif)
        println("failed test_point_set: size(keys) = $(size(keys)) keys = $(keys) dif = $(dif[amax]) query = $(queries[:, amax])} naive_ans = $(a1[amax]) other_ans = $(a2[amax])")
    end
    (p1, q1), (p2, q2)
end

"""gen: takes as input an integer and outputs some number of points"""
function test_point_set_gen(gen, n::Int, q::Int, t1::DataType, t2::DataType)
    keys = gen(T, n)
    values = gen(T, q)
    # println("TEST $n $q")
    test_point_set(keys, values, t1, t2)
end

const MAX_N = 10

@testset "InnerProductMax correctness" begin
    """keys: 3xn"""
    function test_tetrahedron(t1::DataType, t2::DataType)
        keys = tetrahedron(T)
        values = cube(T, 10, 100)
        test_point_set(keys, values, t1, t2)
    end

    """test all points in [-1, 0, 1]"""
    function test_int_tiny(t1::DataType, t2::DataType)
        d = 1
        for n in 4:MAX_N
            for _ in 1:100
                test_point_set_gen((T, n) -> cube(T, d, n), n, 1000, t1, t2)
            end
        end
    end

    function test_int(t1::DataType, t2::DataType)
        for n in 4:MAX_N
            for d in 1:100
                for _ in 1:3
                    test_point_set_gen((T, n) -> cube(T, d, n), n, 1000, t1, t2)
                end
            end
        end
    end

    function test_sphere(t1::DataType, t2::DataType)
        for n in 4:MAX_N
            for _ in 1:100
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
        keys = tetrahedron(T)
        hull = Hull(keys)
        ds = t(hull)
        if t == InnerProductMaxNaive
            @test ds.points == P3[[0.0, 0.0, 0.0], [0.0, 1.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 0.0]]
        end

        @test InnerProductMax.query(ds, P3(1, 2, 3)) == P3(0, 1, 1)
        @test InnerProductMax.query(ds, P3(3, 1, 2)) == P3(1, 0, 1)
        @test InnerProductMax.query(ds, P3(0, 0, 1)) in [P3(1, 0, 1), P3(0, 1, 1)]

        @test InnerProductMax.query(ds, P3(-1, -1, -1)) == P3(0, 0, 0)
        @test InnerProductMax.query(ds, P3(-1, -2, -3)) == P3(0, 0, 0)
        @test InnerProductMax.query(ds, P3(3, -2, -3)) == P3(1, 1, 0)
        @test InnerProductMax.query(ds, P3(1, -2, -3)) == P3(0, 0, 0)

        @test best_answers_and_query_info(ds, keys, vec_to_matrix([
            P3(1, 2, 3),
            P3(-1, -2, -3),
            P3(4, 1, 2),
        ]))[1] == [5, 0, 6]
        # end
    end

    @testset "sanity_check_naive" begin
        sanity_check(InnerProductMaxNaive{T})
    end

    @testset "sanity_check_mine_treap" begin
        sanity_check(InnerProductMaxMine{T,PointLocationDsTreap})
    end

    @testset "sanity_check_mine_rb" begin
        sanity_check(InnerProductMaxMine{T,PointLocationDsRB})
    end

    @testset "sanity_check_nested" begin
        sanity_check(InnerProductMaxNested{T})
    end

    @testset "tiny_treap" begin
        test_int_tiny(InnerProductMaxNaive{T}, InnerProductMaxMine{T,PointLocationDsTreap})
    end

    @testset "tiny_nested" begin
        test_int_tiny(InnerProductMaxNaive{T}, InnerProductMaxNested{T})
    end

    @testset "small_naive" begin
        test_all_small(InnerProductMaxNaive{T}, InnerProductMaxNaive{T})
    end

    @testset "small_mine_treap" begin
        test_all_small(InnerProductMaxNaive{T}, InnerProductMaxMine{T,PointLocationDsTreap})
    end

    @testset "small_mine_rb" begin
        test_all_small(InnerProductMaxNaive{T}, InnerProductMaxMine{T,PointLocationDsRB})
    end

    @testset "small_nested" begin
        test_all_small(InnerProductMaxNaive{T}, InnerProductMaxNested{T})
    end

    @testset "tricky_rb" begin # -0. vs 0.
        test_point_set([-1.0 -4.0 -3.0 0.0; -4.0 3.0 2.0 -3.0; 1.0 -4.0 -4.0 -4.0], hcat([0.0, 1.0, -1.0]),
            InnerProductMaxNaive{T}, InnerProductMaxMine{T,PointLocationDsRB})
    end

    @testset "tricky_tiny" begin # -0. vs 0.
        test_point_set([1.0 1.0 -1.0 1.0; -1.0 1.0 0.0 0.0; 1.0 0.0 1.0 1.0], hcat([1.0, -1.0, 1.0]),
            InnerProductMaxNaive{T}, InnerProductMaxMine{T,PointLocationDsTreap})
    end

    @testset "tricky_treap" begin # parallel rays
        test_point_set([-1.0 -4.0 -3.0 0.0; -4.0 3.0 2.0 -3.0; 1.0 -4.0 -4.0 -4.0], hcat([0.0, 1.0, -1.0]),
            InnerProductMaxNaive{T}, InnerProductMaxMine{T,PointLocationDsTreap})
    end

    @testset "tricky_nested" begin
        for it in 1:100
            # println(it)
            test_point_set([-1.0 1.0 0.0 0.0 -1.0 -1.0; 0.0 1.0 -1.0 -1.0 1.0 0.0; 0.0 1.0 1.0 -1.0 0.0 1.0], hcat([0.0, -1.0, 1.0]),
                InnerProductMaxNaive{T}, InnerProductMaxNested{T})
        end
    end
end