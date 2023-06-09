using Random
using Pandas: DataFrame

"""Returns a Pandas DataFrame comparing the implementations wrt time and memory"""
function sphere_perf(t1::DataType, t2::DataType, nqs::Vector{Tuple{Int,Int}})
    N, Q, method, p_time, p_bytes, q_time = [], [], [], [], [], []
    for (n, q) in nqs
        function add_row(t, pi, qi)
            cmethod = string(t)[length("InnerProductMax.InnerProductMax")+1:end]
            mb = pi.bytes * 1e-6
            push!(N, n)
            push!(Q, q)
            push!(method, cmethod)
            push!(p_time, pi.time)
            push!(p_bytes, mb)
            push!(q_time, qi.time)
            println("Added Row $n $q $cmethod $(pi.time) $(mb) $(qi.time)")
        end
        (p1, q1), (p2, q2) = test_point_set_gen(unit_sphere, n, q, t1, t2)
        add_row(t1, p1, q1)
        add_row(t2, p2, q2)
    end
    cols = ["N", "Q", "Method", "Preproc Time (s)", "Preproc Memory (MB)", "Query Time (s)"]
    d = Dict(k => v for (k, v) in zip(cols, [N, Q, method, p_time, p_bytes, q_time]))
    df = DataFrame(d) # for some reason the column order isn't preserved ...
    df[cols]
end

@testset "sphere_perf_treap" begin
    Random.seed!(1234)
    nqs = [(n, n) for n in Int.([1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5])]
    df = sphere_perf(InnerProductMaxNaive{T}, InnerProductMaxMine{T,PointLocationDsTreap}, nqs)
    println(df)
end

@testset "sphere_perf_rb" begin
    Random.seed!(1234)
    nqs = [(n, n) for n in Int.([1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5])]
    df = sphere_perf(InnerProductMaxNaive{T}, InnerProductMaxMine{T,PointLocationDsRB}, nqs)
    println(df)
end

@testset "sphere_perf_nested" begin
    Random.seed!(1234)
    nqs = [(n, n) for n in Int.([1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5])]
    df = sphere_perf(InnerProductMaxNaive{T}, InnerProductMaxNested{T}, nqs)
    println(df)
end

@test_skip @testset "sphere_perf_5e4" begin
    Random.seed!(1234)
    nqs = [(n, n) for n in Int.([5e4])]
    df = sphere_perf(InnerProductMaxNaive{T}, InnerProductMaxMine{T,PointLocationDsRB}, nqs)
    println(df)
end