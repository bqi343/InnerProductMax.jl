include("persistent_rb.jl")

module PLRB

import ..query_pl
export PointLocationDsRB

using InnerProductMax: PointLocationDs, AugmentedEdge, min_x, max_x
using InnerProductMax.PersistentRB, GeometryBasics
RBNodeEdge = RBNode{T,AugmentedEdge{T}} where {T}

"""
Point Location Data Structure on 2D Mesh, using red-black tree
"""
struct PointLocationDsRB{T} <: PointLocationDs{T}
    roots::Vector{Tuple{T,RBNodeEdge{T}}}
    function PointLocationDsRB{T}(edges::Vector{AugmentedEdge{T}}) where {T}
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
        root::RBNodeEdge{T} = nothing
        roots = Tuple{T,RBNodeEdge}[(-Inf64, root)]
        # println("START")
        while i <= length(events)
            j = i
            # println("OOP ", events[i][1])
            while i <= length(events) && events[i][1] == events[j][1]
                # println("GOT ", events[i][1], " ", events[i][2], " ", events[i][3])
                if events[i][2] == -1
                    # println("DEL")
                    root = delete_root(root, events[i][1], events[i][3])
                else
                    # println("INS")
                    root = insert_root(root, events[i][1], events[i][3])
                end
                # println("CHECKING")
                # if root !== nothing
                #     println(get_colors(root))
                #     tour(root, events[j][1], true)
                # else
                #     println("NO ROOT")
                # end
                i += 1
            end
            if isempty(roots) || roots[end][2] !== root # optimization: root doesn't change that often
                push!(roots, (events[j][1], root))
            else
                # println("saved ...")
            end
        end
        new{T}(roots)
    end
end

query_pl(ds::PointLocationDsRB{T}, p::Point2{T}) where {T} = begin
    # println("Query RB ", p)
    x, root = ds.roots[searchsortedlast(ds.roots, p[1], by=x -> x[1], lt=(x, y) -> x < y)]
    # println("x = ", x)
    # println("tour = ", tour(root, p[1]))
    @assert root !== nothing
    edge = lower_bound(root, p)
    # println("edge = ", edge)
    if edge === nothing
        # println("Case A ", get_last(root, p[1]).vert_above)
        get_last(root, p[1]).vert_above
    else
        # println("Case B ", edge.vert_below)
        edge.vert_below
    end
end

end