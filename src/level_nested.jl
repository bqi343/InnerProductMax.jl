using DataStructures

"""args: points in ccw order"""
function outward_simplices(subset_points::Vector{Point3{T}}) where {T}
    @assert length(subset_points) >= 3
    n_span = length(compute_span(subset_points))
    @assert n_span >= 3
    if n_span == 3
        return [[1, i, i + 1] for i in 2:length(subset_points)-1]
    end
    subset_hull = Hull{T}(vec_to_matrix(subset_points))
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
    marked = Set{V_id}()
    simplices = Vector{V_id}[]
    umbrella = Dict{Tuple{V_id,V_id},V_id}()
    in_indep_set = V_id[]
    for v in vertices(subhull)
        adj_verts = neighbors(subhull, v)
        if length(adj_verts) > 9
            push!(marked, v)
        end
        if v in marked
            continue
        end
        for w in adj_verts
            push!(marked, w)
        end
        push!(in_indep_set, v)
        # println("construct indep ", v, " ", adj_verts)
        # otherwise: construct an independent set
        subset_points = points(subhull, adj_verts)
        new_simplices = outward_simplices(subset_points)
        # println("new_simplices ", new_simplices)
        @assert length(new_simplices) == length(adj_verts) - 2 "wrong number of simplices"
        start_simplices = copy(simplices)
        for simplex in new_simplices
            # println("new simplex ", [adj_verts[simplex[1]], adj_verts[simplex[2]], adj_verts[simplex[3]]])
            push!(simplices, [adj_verts[simplex[1]], adj_verts[simplex[2]], adj_verts[simplex[3]]])
            for j in 1:3
                a, b = simplex[j], simplex[mod1(j + 1, 3)]
                if b == mod1(a + 1, length(subset_points)) # edge is on boundary
                    continue
                end
                # @assert !haskey(umbrella, (adj_verts[a], adj_verts[b])) "points = $(subhull.points)\nvertices = $(vertices(subhull))\nkeys(umbrella) = $(keys(umbrella))\nduplicated edge = $(adj_verts[a]) $(adj_verts[b])\nsimplices = $(simplices)\nstart_simplices = $(start_simplices)\nin_indep_set = $(in_indep_set)\nedges = $(keys(subhull.pair_to_edge))"
                # note: key may already exist in umbrella if all points coplanar!
                umbrella[(adj_verts[a], adj_verts[b])] = v
            end
        end
    end
    # add unchanged simplices
    for ((a, b), e) in subhull.pair_to_edge
        c = e.next.next.v
        # println("Checking simplex ", a, " ", b, " ", c)
        if c < min(a, b) && all(x in marked for x in (a, b, c))
            # println("Adding simplex")
            push!(simplices, [a, b, c])
        end
    end
    # println("done")
    reduce(hcat, simplices), umbrella
end

struct NestedLevel{T}
    outer_subhull::Subhull{T}
    inner_subhull::Subhull{T}
    umbrella::Dict{Tuple{V_id,V_id},V_id}
    function NestedLevel{T}(outer_subhull::Subhull{T}) where {T}
        simplices, umbrella = remove_indep_set(outer_subhull)
        new{T}(outer_subhull, Subhull{T}(outer_subhull.points, simplices), umbrella)
    end
end

function parents(level::NestedLevel, a::V_id, b::V_id)
    if has_edge(level.outer_subhull, a, b)
        adj_verts(level.outer_subhull, a, b)
    else
        @assert haskey(level.umbrella, (a, b)) "edge $((a,b)) doesn't exist\nouter: $(keys(level.outer_subhull.pair_to_edge))\ninner: $(keys(level.inner_subhull.pair_to_edge))"
        [level.umbrella[(a, b)]]
    end
end

get_value(level::NestedLevel{T}, z_dir::Point3{T}, m::V_id) where {T} = dot(point(level.outer_subhull, m), z_dir)

num_a = 0
num_b = 0

function advance(level::NestedLevel{T}, z_dir::Point3{T}, x_dir::Point3{T}, l::V_id, m::V_id, r::V_id) where {T}
    best = (get_value(level, z_dir, m), m)
    for v in vcat(parents(level, l, m), parents(level, m, r))
        cur = get_value(level, z_dir, v)
        if cur > best[1]
            best = (cur, v)
        end
    end
    if best[2] != m
        m = best[2]
        global num_a += 1
        brute_force_triple(level.outer_subhull, z_dir, x_dir, m)
    else
        global num_b += 1
        if num_b % 50000 == 0
            println("num_a = $num_a num_b = $num_b")
        end
        filter_verts(verts) = [v for v in verts if haskey(level.outer_subhull.pair_to_edge, (m, v))]
        best_l = brute_force_r(level.outer_subhull, z_dir, -x_dir, m, filter_verts(vcat([l], parents(level, l, m))))
        best_r = brute_force_r(level.outer_subhull, z_dir, x_dir, m, filter_verts(vcat([r], parents(level, m, r))))

        # @assert haskey(level.outer_subhull.pair_to_edge, (best_l, m))
        # @assert haskey(level.outer_subhull.pair_to_edge, (m, best_r)) "z_dir = $z_dir x_dir = $x_dir\nright edge $((m, best_r[2])) doesn't exist\nparents = $(parents(level, m, r))\npoints = $(level.outer_subhull.points)\n(l, m, r) = ($l, $m, $r)\nouter = $(vertices(level.outer_subhull)), $(keys(level.outer_subhull.pair_to_edge))\ninner = $(vertices(level.inner_subhull)), $(keys(level.inner_subhull.pair_to_edge))"
        best_l, m, best_r
    end
end