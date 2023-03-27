"""Convert 3D into 2D before running qhull. Returns: vertices of hull in CCW order"""
function hull2d(points::Vector{Point3{T}}) where {T}
    span, _ = compute_span(points)
    @assert length(span) == 3
    d1 = normalized(span[1] - span[2])
    d2 = span[1] - span[3]
    d2 -= dot(d2, d1) * d1
    d2 = normalized(d2)
    points = [[dot(p, d1), dot(p, d2)] for p in points]
    points = Matrix{T}(transpose(reduce(hcat, points))) # n x 2
    chull(points).vertices
end

struct Level2D{T}
    s::Subhull{T}
    in_order::Vector{V_id}
    Level2D{T}(s::Subhull{T}) where {T} = new(s, hull2d(s.points))
end

function init_m_index_naive(level::Level2D{T}, z_dir::Point3{T}) where {T}
    n = length(level.in_order)
    eval_at(i) = dot(z_dir, level.s.points[level.in_order[mod1(i, n)]])
    best = 1
    for i in 1:n
        if eval_at(i) > eval_at(best)
            best = i
        end
    end
    best
end

"""returns: index of any maximum"""
function init_m_index(level::Level2D{T}, z_dir::Point3{T}, compare_against_naive::Bool=false) where {T}
    n = length(level.in_order)
    eval_at(i) = dot(z_dir, level.s.points[level.in_order[mod1(i, n)]])
    best_1 = if eval_at(1) < eval_at(2)
        m = findlast(i -> eval_at(i) >= eval_at(1), 1:n)
        # println("AA ", m)
        findfirst(i -> i == m || eval_at(i) >= eval_at(i + 1), 1:m)
    else
        m = findfirst(i -> i > n || eval_at(i) > eval_at(1), 1:n+1)
        # println("BB ", m)
        res = findfirst(i -> i == n + 1 || eval_at(i) >= eval_at(i + 1), m:n+1) + m-1
        # println("HA ", res)
        mod1(res, n)
    end
    if compare_against_naive
        best_2 = init_m_index_naive(level, z_dir)
        eval_res = [eval_at(i) for i in 1:n]
        @assert eval_at(best_2) <= eval_at(best_1) "wrong best $(level.s.points) $(level.in_order) $(z_dir) $(best_1) $(best_2) $(eval_res)"
    end
    best_1
end

"""returns: l, m, r"""
function init_triple(level::Level2D{T}, z_dir::Point3{T}, x_dir::Point3{T}) where {T}
    m_index = init_m_index(level, z_dir)
    n = length(level.in_order)
    l, m, r = level.in_order[mod1(m_index - 1, n)], level.in_order[m_index], level.in_order[mod1(m_index + 1, n)]
    brute_force_r(level.s, z_dir, -x_dir, m, [l, r]), m, brute_force_r(level.s, z_dir, x_dir, m, [l, r])
end
