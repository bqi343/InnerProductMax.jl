using GeometryBasics
using InnerProductMax: collinear, coplanar, num_indep

@testset "geo_utils" begin
    @testset "collinear" begin
        @test !collinear(Point3(1, 0, 0), Point3(0, 1, 0), Point3(0, 0, 1))
        @test collinear(Point3(1, 0, 0), Point3(0, 1, 2), Point3(-1, 2, 4))
    end

    @testset "coplanar" begin
        @test !coplanar(Point3(0, 0, 0), Point3(1, 0, 0), Point3(0, 1, 0), Point3(0, 0, 1))
        @test coplanar(Point3(0, 0, 10), Point3(1, 0, 12), Point3(0, 1, 13), Point3(1, 1, 15))
    end

    @testset "num_indep" begin
        @test num_indep([Point3(1, 0, 0)]) == 1
        @test num_indep([Point3(1, 0, 0), Point3(1, 0, 0)]) == 1
        @test num_indep([Point3(1, 0, 0), Point3(0, 1, 0), Point3(0, 0, 1)]) == 3
        @test num_indep([Point3(1, 0, 0), Point3(0, 1, 2), Point3(-1, 2, 4)]) == 2
        @test num_indep([Point3(0, 0, 0), Point3(1, 0, 0), Point3(0, 1, 0), Point3(0, 0, 1)]) == 4
        @test num_indep([Point3(0, 0, 10), Point3(1, 0, 12), Point3(0, 1, 13), Point3(1, 1, 15)]) == 3
    end
end