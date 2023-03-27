using GeometryBasics
using InnerProductMax: Segment, Ray, get_y_coord, min_x, max_x, normalized

@testset "edge" begin
    @testset "segment" begin
        s = Segment{T}(Point2{T}(3, 4), Point2{T}(5, 5))
        @test get_y_coord(s, 6.0) ≈ 5.5
        @test min_x(s) == 3
        @test max_x(s) == 5
    end

    @testset "ray" begin
        dir = Point2{T}(2, 1)
        dir = normalized(dir)
        r = Ray{T}(Point2{T}(3, 4), dir)

        @test get_y_coord(r, 6.0) ≈ 5.5
        @test min_x(r) == 3
        @test max_x(r) == Inf

        r2 = Ray{T}(Point2{T}(3, 4), -dir)
        @test min_x(r2) == -Inf
        @test max_x(r2) == 3

        dir = Point2{T}(2, -1)
        dir = normalized(dir)
        r3 = Ray{T}(Point2{T}(3, 4), dir)
        @test r3 < r
        @test !(r < r3)
        @test !(r3 < r3)
    end
end