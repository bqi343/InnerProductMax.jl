export InnerProductMaxNested

include("subhull.jl")
include("level_2d.jl")
include("level_nested.jl")

"""
Preproc Time: O(n)
Preproc Mem: O(n)
Query Time: O(log n)

References: 
 - Edelsbrunner, Herbert, and Hermann A. Maurer. "Finding extreme points in three dimensions and solving the post-office problem in the plane." Information processing letters 21.1 (1985): 39-47.
 - o'Rourke, Joseph. Computational geometry in C. Cambridge university press, 1998.
"""
struct InnerProductMaxNested{T} <: AbstractInnerProductMax{T}
    levels::Vector{NestedLevel{T}}
    last_level::Level2D{T}
    function InnerProductMaxNested{T}(hull::Hull{T}) where {T}
        n = size(hull.points, 2)
        old_to_new = zeros(V_id, n)
        for (i, v) in enumerate(hull.vertices)
            old_to_new[v] = i
        end
        cur_subhull = Subhull(Point3.(eachcol(hull.points[:, hull.vertices])), old_to_new[hull.simplices])
        levels = NestedLevel{T}[]
        while true
            # println(length(vertices(cur_subhull)))
            span, _ = compute_span(cur_subhull.points)
            if length(span) < 4
                @assert length(span) == 3
                last_level = Level2D{T}(cur_subhull)
                return new(levels, last_level)
            end
            nested_level, inner_subhull = make_nested_level(cur_subhull)
            push!(levels, nested_level)
            cur_subhull = inner_subhull
        end
    end
end

function query(ds::InnerProductMaxNested{T}, p::Point3{T}) where {T}
    if norm(p) == 0
        p = Point3{T}(0, 0, 1)
    end
    z_dir = normalized(p)
    x_dir = rand_perp(z_dir)
    l, m, r = init_triple(ds.last_level, z_dir, x_dir)
    for level in 0:length(ds.levels)-1
        l, m, r = advance(ds.levels[end-level], z_dir, x_dir, l, m, r)
    end
    ds.levels[1].outer_subhull.points[m]
end