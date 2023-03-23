using Random

@testset "distributions_3d" begin
    @testset "sphere" begin
        Random.seed!(1234)
        @test isapprox(unit_sphere(5), [0.589154 -0.0209547 0.869464 0.426918 0.577249
                -0.594351 -0.383788 0.489562 -0.672442 -0.0358788
                0.547398 -0.923183 -0.0660363 0.604618 0.81578], atol=1e-4)
    end

    @testset "cube" begin
        Random.seed!(1234)
        @test cube(1, 5) == [-1 1 1 1 -1
            0 0 1 0 -1
            -1 0 0 1 0]
    end
end