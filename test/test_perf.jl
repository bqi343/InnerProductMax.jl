import Random

"""
n = q = 100000:
InnerProductMax.InnerProductMaxNaive{Float64} p1 0.001996282 q1 18.183047406
InnerProductMax.InnerProductMaxMine{Float64} p2 2.51332117 q2 0.444145011
"""
function sphere_perf(t1, t2)
    n = 100000
    q = 100000
    (p1, q1), (p2, q2) = test_point_set_gen(unit_sphere, n, q, t1, t2)
    println(t1, " p1 ", p1, " q1 ", q1)
    println(t2, " p2 ", p2, " q2 ", q2)
end

@testset "sphere_perf" begin
    Random.seed!(1234)
    sphere_perf(InnerProductMaxNaive{Float64}, InnerProductMaxMine{Float64})
end