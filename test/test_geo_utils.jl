using GeometryBasics
using InnerProductMax: collinear, coplanar, span_len, compute_span, diameter

@testset "geo_utils" begin
    @testset "collinear" begin
        @test !collinear(Point3(1, 0, 0), Point3(0, 1, 0), Point3(0, 0, 1))
        @test collinear(Point3(1, 0, 0), Point3(0, 1, 2), Point3(-1, 2, 4))
    end

    @testset "coplanar" begin
        @test !coplanar(Point3(0, 0, 0), Point3(1, 0, 0), Point3(0, 1, 0), Point3(0, 0, 1))
        @test coplanar(Point3(0, 0, 10), Point3(1, 0, 12), Point3(0, 1, 13), Point3(1, 1, 15))
    end

    @testset "span" begin
        @test span_len([Point3(1, 0, 0)]) == 1
        @test span_len([Point3(1, 0, 0), Point3(1, 0, 0)]) == 1
        @test compute_span([Point3(1, 0, 0), Point3(1, 0, 0)]) == ([Point3(1, 0, 0)], [1])
        @test span_len([Point3(1, 0, 0), Point3(0, 1, 0), Point3(0, 0, 1)]) == 3
        @test span_len([Point3(1, 0, 0), Point3(0, 1, 2), Point3(-1, 2, 4)]) == 2
        @test span_len([Point3(0, 0, 0), Point3(1, 0, 0), Point3(0, 1, 0), Point3(0, 0, 1)]) == 4
        @test span_len([Point3(0, 0, 10), Point3(1, 0, 12), Point3(0, 1, 13), Point3(1, 1, 15)]) == 3
    end

    @testset "diameter" begin
        @test diameter(vec_to_matrix([Point3(1, 0, 3), Point3(2, 0, 13)])) == 10
        @test diameter(vec_to_matrix([Point3(2, 0, 13), Point3(1, 0, 3)])) == 10
        @test diameter(vec_to_matrix([Point3(1, 0, 3), Point3(2, 0, 3)])) == 1
    end
end