using InnerProductMax: HullLite, volume

@testset "hull_lite" begin
    points = P3[P3(-1, -1, 0), P3(1, -1, 0), P3(0, 1, 0), P3(0, 0, 1), P3(0, 0, 7)] .+ P3(3, 5, 10)
    h = HullLite{T}(points)
    @test h.points == points
    @test size(h.simplices, 2) == 4
    @test volume(h) ≈ 4.0 * 7 / 6
end