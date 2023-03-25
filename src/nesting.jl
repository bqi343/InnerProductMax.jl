export InnerProductMaxNested

include("subhull.jl")
include("level_2d.jl")
include("level_nested.jl")

"""
Preproc Time: O(n)
Preproc Mem: O(n)
Query Time: O(log n)

References: Edelsbrunner  & Maurer  (1985); 
Edelsbrunner (1987, Section 9.5.3); 
Computational Geometry in C 7.10.5
"""
struct InnerProductMaxNested{T} <: AbstractInnerProductMax{T}
    hull::Hull{T}
    levels::Vector{NestedLevel{T}}
    last_level::Level2D{T}
    function InnerProductMaxNested{T}(hull::Hull{T}) where {T}
        cur_subhull = Subhull{T}(hull.points, hull.simplices)
        levels = NestedLevel{T}[]
        while true
            # println(length(vertices(cur_subhull)))
            span = compute_span(hull.points[:, Int.(vertices(cur_subhull))])
            if length(span) < 4
                @assert length(span) == 3
                last_level = Level2D{T}(cur_subhull)
                return new(hull, levels, last_level)
            end
            push!(levels, NestedLevel{T}(cur_subhull))
            cur_subhull = levels[end].inner_subhull
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
    # println("AA ", typeof(l))
    # println("LEN ", length(ds.levels))
    # println(ds.hull.points)
    # println("start query")
    # println("init ", l, " ", m, " ", r)
    # println("z_dir = ", z_dir, " x_dir = ", x_dir)
    for level in 0:length(ds.levels)-1
        # println("before ", vertices(ds.levels[end-level].inner_subhull))
        l, m, r = advance(ds.levels[end-level], z_dir, x_dir, l, m, r)
        # println("after ", vertices(ds.levels[end-level].outer_subhull))
        # println("advance ", l, " ", m, " ", r)
    end
    # println("DONE")
    point(ds.levels[1].outer_subhull, m)
end