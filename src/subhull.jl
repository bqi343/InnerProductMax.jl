using GeometryBasics

"""TODO: add reference for half-edge"""
mutable struct HalfEdge # a -> b
    const v::V_id
    flip::Union{HalfEdge,Nothing} # b -> a
    next::Union{HalfEdge,Nothing} # b- > c
    HalfEdge(v::V_id) = new(v, nothing, nothing)
end

"""Data structure for a subset of a 3D convex hull."""
struct Subhull{T}
    points::Vector{Point3{T}}
    pair_to_edge::Dict{Tuple{V_id,V_id},HalfEdge}
    vert_to_any_edge::Vector{HalfEdge}
    function Subhull(points::Vector{Point3{T}}, simplices::Matrix{V_id}) where {T}
        pair_to_edge = Dict{Tuple{V_id,V_id},HalfEdge}()
        # create halfedges and next pointers
        for (a, b, c) in eachcol(simplices)
            ha, hb, hc = HalfEdge(a), HalfEdge(b), HalfEdge(c)
            ha.next = hb
            hb.next = hc
            hc.next = ha
            pair_to_edge[(a, b)] = ha
            pair_to_edge[(b, c)] = hb
            pair_to_edge[(c, a)] = hc
        end
        vert_to_any_edge = Vector{HalfEdge}(undef, length(points))
        # create flip pointers and vert_to_any_edge
        for ((a, b), e) in pair_to_edge
            e.flip = pair_to_edge[(b, a)]
            vert_to_any_edge[a] = e
        end
        new{T}(points, pair_to_edge, vert_to_any_edge)
    end
end

has_edge(subhull::Subhull, a::V_id, b::V_id) = haskey(subhull.pair_to_edge, (a, b))
"""return: the two vertices that share a face with this edge"""
function adj_verts(subhull::Subhull, a::V_id, b::V_id)::Vector{V_id}
    edge = subhull.pair_to_edge[(a, b)]
    [edge.next.next.v, edge.flip.next.next.v]
end

"""list of all neighbors of v"""
function neighbors(s::Subhull, m::V_id)::Vector{V_id}
    res = V_id[]
    start_edge = s.vert_to_any_edge[m]
    cur_edge = start_edge
    while true
        push!(res, cur_edge.flip.v) # TODO: speedup?
        cur_edge = cur_edge.next.next.flip
        cur_edge !== start_edge || break
    end
    res
end

"""find r by iterating through all cands"""
function brute_force_r(s::Subhull{T}, z_dir::Point3{T}, x_dir::Point3{T}, m::V_id, cands::Vector{V_id}) where {T}
    best = (Point2{T}(0, 0), 0)
    for r in cands
        offset = s.points[r] - s.points[m]
        x = dot(offset, x_dir)
        z = dot(offset, z_dir)
        EPS = 1e-12
        @assert z < EPS "m is not optimal $z"
        cur = Point2{T}(x, z)
        if GeometryBasics.cross(best[1], cur) >= 0
            best = (cur, r)
        end
    end
    best[2]
end

"""Returns l, m, r"""
function brute_force_triple(s::Subhull{T}, z_dir::Point3{T}, x_dir::Point3{T}, m::V_id) where {T}
    cands = neighbors(s, m)
    @assert length(cands) <= 9 "adj list too large to brute force"
    brute_force_r(s, z_dir, -x_dir, m, cands), m, brute_force_r(s, z_dir, x_dir, m, cands)
end