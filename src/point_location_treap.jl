include("persistent_treap.jl")

module PLTreap

import ..query_pl
export PointLocationDsTreap

using InnerProductMax: PointLocationDs, AugmentedEdge, min_x, max_x
using InnerProductMax.PersistentTreap, GeometryBasics
TNodeEdge = TNode{AugmentedEdge{T}} where {T}

"""
Point Location Data Structure on 2D Mesh
"""
struct PointLocationDsTreap{T} <: PointLocationDs{T}
    roots::Vector{Tuple{T,TNodeEdge}}
    function PointLocationDsTreap{T}(edges::Vector{AugmentedEdge{T}}) where {T}
        # println("PointLocationDsTreap")
        # println(edges)
        events = Tuple{T,Int,AugmentedEdge{T}}[]
        # println("Init")
        for e in edges
            # println(e.e, " ", e.vert_below, " ", e.vert_above, " ", min_x(e), " ", max_x(e))
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
        # println("Events")
        # for e in events
        #     println(e[1]," ",e[2])
        # end
        i = 1
        root::TNodeEdge = nothing
        roots = Tuple{T,TNodeEdge}[(-Inf64, root)]
        # for i in 1:length(events)
        #     if events[i][1] == 0
        #         println(i, " ", events[i][1], " ", events[i][2])
        #     end
        #     for j in 1:length(events)
        #         if events[i][1] == 0 && events[j][1] == 0
        #             println(i, " ", j, " ", (events[i][1], events[i][2]) < (events[j][1], events[j][2]))
        #         end
        #     end
        # end
        # for (x, y) in zip(events, events[2:end])
        #     println(x[1] == 0)
        #     println((x[1], x[2]), " ", (y[1], y[2]), " ", x[1] < y[1], " ", x[2] < y[2], " ", (x[1], x[2]) < (y[1], y[2]), " ", (y[1], y[2]) < (x[1], x[2]))
        # end
        # println("Events")
        # for (a, b, c) in events
        #     println(a, " ", b)
        # end
        # println("----")
        while i <= length(events)
            j = i
            # println("Processing ", events[i][1])
            # println("Event ", events[j][1])
            while i <= length(events) && events[i][1] == events[j][1]
                # println("Event ", events[i])
                if events[i][2] == -1
                    root = delete(root, events[i][3])
                else
                    root = insert(root, events[i][3])
                end
                i += 1
            end
            push!(roots, (events[j][1], root))
            # println("State at")
            # println(events[j][1])
            # println(tour(root))
        end
        new{T}(roots)
    end
end

query_pl(ds::PointLocationDsTreap{T}, p::Point2{T}) where {T} = begin
    # println("Query ", p)
    x, root = ds.roots[searchsortedlast(ds.roots, p[1], by=x -> x[1], lt=(x, y) -> x < y)]
    # println("x = ", x)
    # println("tour = ", tour(root))
    @assert root !== nothing
    edge = lower_bound(root, p)
    # println("edge = ", edge)
    if edge === nothing
        get_last(root).vert_above
    else
        edge.vert_below
    end
end

end