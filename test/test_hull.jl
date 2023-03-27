using GeometryBasics
using InnerProductMax: negated_copy

@testset "hull" begin
    @testset "basic" begin
        points = tetrahedron(T)
        hull = Hull(points)
        neg_hull = negated_copy(hull)
        @test hull != neg_hull
    end
end