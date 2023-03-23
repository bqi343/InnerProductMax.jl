include("point_location.jl")

export InnerProductMaxMine

function upward(dir::Point3)
    dir[3] > 0
end

function extract_2d(dir::Point3{T}) where {T}
    if !upward(dir)
        throw(DomainError("must point upward"))
    end
    Point2{T}(dir[1] / dir[3], dir[2] / dir[3])
end

"""Answers maximum inner product queries in the upper half-plane."""
struct UpperDS{T}
    pl_ds::PointLocationDS{T}
    function UpperDS(hull::Hull{T}) where {T}
        # iterate over upward-pointing facets
        augmented_edges = AugmentedEdge{T}[]

        function add_two_way_edge(p1::Point2{T}, p2::Point2{T}, a::Int32, b::Int32, duplicated=true)
            if p1[1] ≈ p2[1]
                return
            end
            if !(p1[1] < p2[1])
                if duplicated
                    return
                else
                    p1, p2 = p2, p1
                    a, b = b, a
                end
            end
            push!(augmented_edges, AugmentedEdge{T}(Segment{T}(p1, p2), hull.points[:, a], hull.points[:, b]))
        end

        function add_one_way_edge(p1::Point3, p2::Point3, a::Int32, b::Int32)
            @assert upward(p1) && !upward(p2)
            point_at_infinity = p1[3] * p2 - p2[3] * p1
            point_at_infinity = point_at_infinity ./ norm(point_at_infinity)
            @assert point_at_infinity[3] ≈ 0
            point_at_infinity = Point3(point_at_infinity[1:2]..., 1e-12) # TODO: adjust?
            add_two_way_edge(extract_2d(p1), extract_2d(point_at_infinity), a, b, false)
        end

        for i in 1:size(hull.simplices, 2)
            cur_facet = Point3(hull.facets[1:3, i])
            if upward(cur_facet)
                # println("Found Facet")
                # println(hull.points[:, hull.simplices[1, i]])
                # println(hull.points[:, hull.simplices[2, i]])
                # println(hull.points[:, hull.simplices[3, i]])
                # println(extract_2d(cur_facet))
                for j in 1:3
                    a = hull.simplices[j, i]
                    b = hull.simplices[rem(j, 3)+1, i]
                    # println("considering edge ", hull.points[:, a], " ", hull.points[:, b])
                    adj_simplex = hull.adj_simplices[j, i]
                    adj_facet = Point3(hull.facets[1:3, adj_simplex])
                    if upward(adj_facet)
                        # println("adj is upward")
                        add_two_way_edge(extract_2d(cur_facet), extract_2d(adj_facet), a, b)
                    else
                        # println("adj is downward")
                        add_one_way_edge(cur_facet, adj_facet, a, b)
                    end
                end
            end
        end
        new{T}(PointLocationDS(augmented_edges))
    end
end

function query(ds::UpperDS, p::Point2)
    query(ds.pl_ds, p)
end

struct InnerProductMaxMine{T} <: AbstractInnerProductMax{T}
    upper::UpperDS{T}
    lower::UpperDS{T}
    function InnerProductMaxMine{T}(hull::Hull{T}) where {T}
        new{T}(UpperDS(hull), UpperDS(negated_copy(hull)))
    end
end

function query(ds::InnerProductMaxMine{T}, p::Point3{T}) where {T}
    if p == zeros(3) # if zero, replace with arbitrary point
        p = Point3{T}(0, 0, 1)
    end
    p = p ./ norm(p)
    eps = 1e-6
    if p[3] >= 0
        x, y, z = p
        z = max(z, eps)
        query(ds.upper, Point2{T}(x / z, y / z))
    else
        x, y, z = -p
        z = max(z, eps)
        -query(ds.lower, Point2{T}(x / z, y / z))
    end
end