include("persistent_treap.jl")

module PLTreap

import ..query_pl
export PointLocationDsTreap

using InnerProductMax: PointLocationDs, AugmentedEdge, min_x, max_x
using InnerProductMax.PersistentTreap, GeometryBasics

TNodeEdge{T} = TNode{AugmentedEdge{T}}

"""
Point Location Data Structure on 2D Mesh
"""
struct PointLocationDsTreap{T} <: PointLocationDs{T}
    roots::Vector{Tuple{T,TNodeEdge{T}}}
    function PointLocationDsTreap{T}(edges::Vector{AugmentedEdge{T}}) where {T}
        events = Tuple{T,Int,AugmentedEdge{T}}[]
        for e in edges
            push!(events, (min_x(e), 1, e))
            push!(events, (max_x(e), -1, e))
        end
        sort!(events, by=x -> (x[1], x[2]), lt=(x, y) -> begin
            if x[1] != y[1]
                x[1] < y[1]
            else
                x[2] < y[2]
            end
        end) # see https://discourse.julialang.org/t/negative-zeros-and-sorting/85252
        i = 1
        root::TNodeEdge{T} = nothing
        roots = Tuple{T,TNodeEdge{T}}[(-Inf, root)]
        while i <= length(events)
            j = i
            while i <= length(events) && events[i][1] == events[j][1]
                if events[i][2] == -1
                    root = delete(root, events[i][3])
                else
                    root = insert(root, events[i][3])
                end
                i += 1
            end
            push!(roots, (events[j][1], root))
        end
        new{T}(roots)
    end
end

query_pl(ds::PointLocationDsTreap{T}, p::Point2{T}) where {T} = begin
    x, root = ds.roots[searchsortedlast(ds.roots, p[1], by=x -> x[1], lt=(x, y) -> x < y)]
    @assert root !== nothing
    edge = lower_bound(root, p)
    if edge === nothing
        get_last(root).vert_above
    else
        edge.vert_below
    end
end

end