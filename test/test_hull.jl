using GeometryBasics
using InnerProductMax: negated_copy

@testset "hull" begin
    @testset "basic" begin
        points = tetrahedron()
        hull = Hull{Float64}(points)
        # @show(hull)
        # println()
        neg_hull = negated_copy(hull)
        # @show(hull)
        # println()
        # @show(neg_hull)
        @test hull != neg_hull
        # println()
    end
end