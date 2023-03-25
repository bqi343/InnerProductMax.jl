using Random
using LinearAlgebra

P3 = Point3{Float64}

function check_keys(keys::Matrix{T}, returned_keys::Matrix{T}) where {T}
    s = Set{Point3{T}}(Point3{T}.(eachcol(keys)))
    @test all(k in s for k in eachcol(returned_keys))
end

"""returns: max dot products, query time"""
function best_answers_and_query_info(ds::AbstractInnerProductMax{T}, keys::Matrix{T}, queries::Matrix{T}) where {T}
    query_info = @timed(begin
        returned_keys = best_vecs(ds, queries)
        0
    end)
    check_keys(keys, returned_keys)
    dot.(eachcol(returned_keys), eachcol(queries)), query_info
end

"""returns: max dot products, preprocess time, query time"""
function best_answers_and_infos(t::DataType, hull::Hull{T}, keys::Matrix{T}, queries::Matrix{T}) where {T}
    preprocess_info = @timed(begin
        ds = t(hull)
        0
    end)
    answers, query_info = best_answers_and_query_info(ds, keys, queries)
    answers, preprocess_info, query_info
end

function test_point_set(keys::Matrix{T}, queries::Matrix{T}, t1::DataType, t2::DataType) where {T<:Real}
    hull = Hull{T}(keys)
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
    keys = gen(n)
    values = gen(q)
    # println("TEST $n $q")
    test_point_set(keys, values, t1, t2)
end

@testset "InnerProductMax correctness" begin
    """keys: 3xn"""

    MAX_N = 10

    function test_tetrahedron(t1::DataType, t2::DataType)
        keys = tetrahedron()
        values = cube(10, 100)
        test_point_set(keys, values, t1, t2)
    end

    """test all points in [-1, 0, 1]"""
    function test_int_tiny(t1::DataType, t2::DataType)
        d = 1
        for n in 4:MAX_N
            for rep in 1:100
                test_point_set_gen(n -> cube(d, n), n, 1000, t1, t2)
            end
        end
    end

    function test_int(t1::DataType, t2::DataType)
        for n in 4:MAX_N
            for d in 1:100
                for rep in 1:3
                    test_point_set_gen(n -> cube(d, n), n, 1000, t1, t2)
                end
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
        @test InnerProductMax.query(ds, P3(1, 2, 3)) == P3(0, 1, 1)
        @test InnerProductMax.query(ds, P3(3, 1, 2)) == P3(1, 0, 1)
        @test InnerProductMax.query(ds, P3(0, 0, 1)) in [P3(1, 0, 1), P3(0, 1, 1)]
        # end
        # @testset "Downward" begin
        @test InnerProductMax.query(ds, P3(-1, -1, -1)) == P3(0, 0, 0)
        @test InnerProductMax.query(ds, P3(-1, -2, -3)) == P3(0, 0, 0)
        @test InnerProductMax.query(ds, P3(3, -2, -3)) == P3(1, 1, 0)
        @test InnerProductMax.query(ds, P3(1, -2, -3)) == P3(0, 0, 0)
        # end
        # @testset "Best Answers" begin
        @test best_answers_and_query_info(ds, keys, vec_to_matrix([
            P3(1, 2, 3),
            P3(-1, -2, -3),
            P3(4, 1, 2),
        ]))[1] == [5, 0, 6]
        # end
    end

    @testset "sanity_check_naive" begin
        sanity_check(InnerProductMaxNaive{Float64})
    end

    @testset "sanity_check_mine_treap" begin
        sanity_check(InnerProductMaxMine{Float64,PointLocationDsTreap})
    end

    @testset "sanity_check_mine_rb" begin
        sanity_check(InnerProductMaxMine{Float64,PointLocationDsRB})
    end

    @testset "sanity_check_nesting" begin
        sanity_check(InnerProductMaxNested{Float64})
    end

    @testset "tiny_treap" begin
        test_int_tiny(InnerProductMaxNaive{Float64}, InnerProductMaxMine{Float64,PointLocationDsTreap})
    end

    @testset "tiny_nesting" begin
        test_int_tiny(InnerProductMaxNaive{Float64}, InnerProductMaxNested{Float64})
    end

    @testset "small_naive" begin
        test_all_small(InnerProductMaxNaive{Float64}, InnerProductMaxNaive{Float64})
    end

    @testset "small_mine_treap" begin
        test_all_small(InnerProductMaxNaive{Float64}, InnerProductMaxMine{Float64,PointLocationDsTreap})
    end

    @testset "small_mine_rb" begin
        test_all_small(InnerProductMaxNaive{Float64}, InnerProductMaxMine{Float64,PointLocationDsRB})
    end

    @testset "small_nesting" begin
        test_all_small(InnerProductMaxNaive{Float64}, InnerProductMaxNested{Float64})
    end

    @testset "tricky_rb" begin # -0. vs 0.
        test_point_set([-1.0 -4.0 -3.0 0.0; -4.0 3.0 2.0 -3.0; 1.0 -4.0 -4.0 -4.0], hcat([0.0, 1.0, -1.0]),
            InnerProductMaxNaive{Float64}, InnerProductMaxMine{Float64,PointLocationDsRB})
    end

    @testset "tricky_tiny" begin # -0. vs 0.
        test_point_set([1.0 1.0 -1.0 1.0; -1.0 1.0 0.0 0.0; 1.0 0.0 1.0 1.0], hcat([1.0, -1.0, 1.0]),
            InnerProductMaxNaive{Float64}, InnerProductMaxMine{Float64,PointLocationDsTreap})
    end

    @testset "tricky_treap" begin # parallel rays
        test_point_set([-1.0 -4.0 -3.0 0.0; -4.0 3.0 2.0 -3.0; 1.0 -4.0 -4.0 -4.0], hcat([0.0, 1.0, -1.0]),
            InnerProductMaxNaive{Float64}, InnerProductMaxMine{Float64,PointLocationDsTreap})
    end

    @testset "tricky_nesting" begin
        for it in 1:100
            println(it)
            test_point_set([-1.0 1.0 0.0 0.0 -1.0 -1.0; 0.0 1.0 -1.0 -1.0 1.0 0.0; 0.0 1.0 1.0 -1.0 0.0 1.0], hcat([0.0, -1.0, 1.0]),
            InnerProductMaxNaive{Float64}, InnerProductMaxNested{Float64})
        end
    end
end