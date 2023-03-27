using DataStructures
include("hull_lite.jl")

"""args: points in ccw order"""
function outward_simplices(subset_points::Vector{Point3{T}}) where {T}
    @assert length(subset_points) >= 3
    n_span = span_len(subset_points)
    @assert n_span >= 3
    if n_span == 3
        return [[1, i, i + 1] for i in 2:length(subset_points)-1]
    end
    subset_hull = HullLite{T}(subset_points)
    n_simplices = size(subset_hull.simplices, 2)
    start_simplex = 0
    for i in 1:n_simplices
        for j in 1:3
            if (subset_hull.simplices[j, i], subset_hull.simplices[mod1(j + 1, 3), i]) == (1, 2)
                start_simplex = i
            end
        end
    end
    @assert start_simplex != 0 "didn't find start simplex"
    q = Queue{Int}()
    visited = zeros(Bool, n_simplices)
    function visit(s)
        if visited[s]
            return
        end
        visited[s] = true
        enqueue!(q, s)
    end
    visit(start_simplex)
    simplices = []
    while !isempty(q)
        s = dequeue!(q)
        push!(simplices, subset_hull.simplices[:, s])
        for j in 1:3
            a = subset_hull.simplices[j, s]
            b = subset_hull.simplices[mod1(j + 1, 3), s]
            if b == mod(a, length(subset_points)) + 1 # on boundary
                continue
            end
            visit(subset_hull.adj_simplices[j, s])
        end
    end
    simplices
end

function remove_indep_set(subhull::Subhull{T}) where {T}
    # println("remove_indep_set ", vertices(subhull))
    n = length(subhull.points)
    marked = zeros(Bool, n)
    simplices = Vector{V_id}[]
    umbrella_pairs = Pair{Tuple{V_id,V_id},V_id}[]
    for v in 1:n
        adj_verts = neighbors(subhull, V_id(v))
        if length(adj_verts) > 9
            marked[v] = true
        end
        if marked[v]
            continue
        end
        # construct an independent set
        for w in adj_verts
            marked[w] = true
        end
        subset_points = subhull.points[adj_verts]
        new_simplices = outward_simplices(subset_points)
        @assert length(new_simplices) == length(adj_verts) - 2 "wrong number of simplices"
        for simplex in new_simplices
            push!(simplices, [adj_verts[simplex[1]], adj_verts[simplex[2]], adj_verts[simplex[3]]])
            for j in 1:3
                a, b = simplex[j], simplex[mod1(j + 1, 3)]
                if b == mod1(a + 1, length(subset_points)) # edge is on boundary
                    continue
                end
                # @assert !haskey(umbrella, (adj_verts[a], adj_verts[b])) "points = $(subhull.points)\nvertices = $(vertices(subhull))\nkeys(umbrella) = $(keys(umbrella))\nduplicated edge = $(adj_verts[a]) $(adj_verts[b])\nsimplices = $(simplices)\nstart_simplices = $(start_simplices)\nin_indep_set = $(in_indep_set)\nedges = $(keys(subhull.pair_to_edge))"
                # note: key may already exist in umbrella if all points coplanar!
                # umbrella[(adj_verts[a], adj_verts[b])] = v
                push!(umbrella_pairs, (adj_verts[a], adj_verts[b]) => v)
            end
        end
    end
    # add unchanged simplices
    for ((a, b), e) in subhull.pair_to_edge
        c = e.next.next.v
        if c < min(a, b) && all(marked[x] for x in (a, b, c))
            push!(simplices, [a, b, c])
        end
    end
    marked, reduce(hcat, simplices), Dict{Tuple{V_id,V_id},V_id}(umbrella_pairs)
end

struct NestedLevel{T}
    outer_subhull::Subhull{T}
    umbrella::Dict{Tuple{V_id,V_id},V_id}
    new_to_old::Vector{V_id}
end

"""returns: nested level, inner subhull"""
function make_nested_level(outer_subhull::Subhull{T}) where {T}
    marked, simplices, umbrella = remove_indep_set(outer_subhull)
    n = length(outer_subhull.points)
    new_to_old = [i for i in 1:n if marked[i]]
    old_to_new = zeros(V_id, n)
    for (i, v) in enumerate(new_to_old)
        old_to_new[v] = i
    end
    NestedLevel{T}(outer_subhull, umbrella, new_to_old), Subhull(outer_subhull.points[new_to_old], old_to_new[simplices])
end

function parents(level::NestedLevel, a::V_id, b::V_id)::Vector{V_id}
    if has_edge(level.outer_subhull, a, b)
        adj_verts(level.outer_subhull, a, b)
    else
        # @assert haskey(level.umbrella, (a, b)) "edge $((a,b)) doesn't exist\nouter: $(keys(level.outer_subhull.pair_to_edge))\ninner: $(keys(level.inner_subhull.pair_to_edge))"
        [level.umbrella[(a, b)]]
    end
end

get_value(level::NestedLevel{T}, z_dir::Point3{T}, m::V_id) where {T} = dot(level.outer_subhull.points[m], z_dir)

function advance(level::NestedLevel{T}, z_dir::Point3{T}, x_dir::Point3{T}, l::V_id, m::V_id, r::V_id) where {T}
    l = level.new_to_old[l]
    m = level.new_to_old[m]
    r = level.new_to_old[r]
    best = (get_value(level, z_dir, m), m)
    parents_l = parents(level, l, m)
    parents_r = parents(level, m, r)
    for v in vcat(parents_l, parents_r)
        cur = get_value(level, z_dir, v)
        if cur > best[1]
            best = (cur, v)
        end
    end
    if best[2] != m
        m = best[2]
        # global num_a += 1
        brute_force_triple(level.outer_subhull, z_dir, x_dir, m)
    else
        # global num_b += 1
        # if num_b % 50000 == 0
        #     println("num_a = $num_a num_b = $num_b")
        # end
        # brute_force_r(level.outer_subhull, z_dir, -x_dir, m, filter_verts(l, parents_l))
        # filter_verts(v, par_v) = [v for v in verts if haskey(level.outer_subhull.pair_to_edge, (m, v))]
        function find_best(x_dir::Point3{T}, v::V_id, par_v::Vector{V_id})
            if length(par_v) == 1 # umbrella parent
                par_v[1]
            else # two adjacent parents, should compare both against v
                brute_force_r(level.outer_subhull, z_dir, x_dir, m, vcat(par_v, [v]))
            end
        end
        best_l = find_best(-x_dir, l, parents_l)
        best_r = find_best(x_dir, r, parents_r)

        # @assert haskey(level.outer_subhull.pair_to_edge, (best_l, m))
        # @assert haskey(level.outer_subhull.pair_to_edge, (m, best_r)) "z_dir = $z_dir x_dir = $x_dir\nright edge $((m, best_r[2])) doesn't exist\nparents = $(parents(level, m, r))\npoints = $(level.outer_subhull.points)\n(l, m, r) = ($l, $m, $r)\nouter = $(vertices(level.outer_subhull)), $(keys(level.outer_subhull.pair_to_edge))\ninner = $(vertices(level.inner_subhull)), $(keys(level.inner_subhull.pair_to_edge))"
        best_l, m, best_r
    end
end