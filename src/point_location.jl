include("bst.jl")
# include("persistent_red_black_tree.jl")
include("persistent_treap.jl")
include("edge.jl")

BNode = Union{BstNode{AugmentedEdge{T}},Nothing} where {T};

"""
Point Location Data Structure on 2D Mesh
"""
struct PointLocationDS{T<:Real}
    roots::Vector{Tuple{T,BNode}}
    function PointLocationDS(edges::Vector{AugmentedEdge{T}}) where {T}
        # println("PointLocationDS")
        # println(edges)
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
        # println("Events")
        # for e in events
        #     println(e[1]," ",e[2])
        # end
        i = 1
        root::TNode{AugmentedEdge{T}} = nothing
        roots = Tuple{T,BNode}[(-Inf64, root)]
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
        # println(-0.0 < 0.0)
        # println((-0.0, 1) < (0.0, -1))
        # println((0.0, -1) < (-0.0, 1))
        # for (a, b, c) in events
        #     println(a, " ", b)
        # end
        # println("----")
        while i <= length(events)
            j = i
            # println("Event ", events[j][1])
            while i <= length(events) && events[i][1] == events[j][1]
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

function query(ds::PointLocationDS{T}, p::Point2) where {T}
    # println("Query ", p)
    x, root = ds.roots[searchsortedlast(ds.roots, p[1], by=x -> x[1])]
    # println("x = ", x)
    # println("tour = ", tour(root))
    @assert root !== nothing
    edge = lower_bound(root, p)
    # println("edge = ", edge)
    if edge === nothing
        last(root).vert_above
    else
        edge.vert_below
    end
end