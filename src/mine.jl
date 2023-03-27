include("point_location.jl")

"""whether vector points in positive z-dir"""
upward(dir::Point3) = dir[3] > 0
function extract_2d(dir::Point3{T}) where {T}
    if !upward(dir)
        throw(DomainError("must point upward"))
    end
    Point2{T}(dir[1] / dir[3], dir[2] / dir[3])
end

"""Answers maximum inner product queries in the upper half-plane."""
struct UpperDS{T,P}
    pl_ds::P
    function UpperDS{T,P}(hull::Hull{T}) where {T,P}
        # iterate over upward-pointing facets
        augmented_edges = AugmentedEdge{T}[]

        function add_two_way_edge(p1::Point2{T}, p2::Point2{T}, a::V_id, b::V_id)
            if !(p1[1] < p2[1]) || p1[1] ≈ p2[1]
                return
            end
            # println("Adding Two-Way Edge $(p1) $(p2) $(hull.points[:, a]) $(hull.points[:, b])")
            s = Segment{T}(p1, p2)
            # println("Seg ", s)
            push!(augmented_edges, AugmentedEdge{T}(s, hull.points[:, a], hull.points[:, b]))
        end

        function add_one_way_edge(p1::Point3{T}, p2::Point3{T}, a::V_id, b::V_id)
            # println("One Way ", p1, " ", p2)
            @assert upward(p1) && !upward(p2)
            point_at_infinity = p1[3] * p2 - p2[3] * p1
            point_at_infinity = point_at_infinity ./ norm(point_at_infinity)
            @assert point_at_infinity[3] ≈ 0
            dir = Point2{T}(point_at_infinity[1:2])
            dir = normalized(dir)
            if dir[1] ≈ 0
                return
            elseif dir[1] < 0
                a, b = b, a
            end
            r = Ray{T}(extract_2d(p1), Point2{T}(dir))
            # println("Ray ", r)
            push!(augmented_edges, AugmentedEdge{T}(r, hull.points[:, a], hull.points[:, b]))
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
                    b = hull.simplices[mod1(j + 1, 3), i]
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
        new{T,P}(P{T}(augmented_edges))
    end
end

query(ds::UpperDS{T,P}, p::Point2{T}) where {P,T} = query_pl(ds.pl_ds, p)

"""
Preproc Time: O(n log n)
Preproc Mem: O(n) if red-black tree, O(n log n) if treap
Query Time: O(log n)
"""
struct InnerProductMaxMine{T,P} <: AbstractInnerProductMax{T}
    upper::UpperDS{T,P}
    lower::UpperDS{T,P}
    function InnerProductMaxMine{T,P}(hull::Hull{T}) where {T,P}
        new{T,P}(UpperDS{T,P}(hull), UpperDS{T,P}(negated_copy(hull)))
    end
end

function query(ds::InnerProductMaxMine{T,P}, p::Point3{T}, eps::T=1e-9) where {T,P}
    if p == zeros(3) # if zero, replace with arbitrary point
        p = Point3{T}(0, 0, 1)
    end
    p = normalized(p)
    if p[3] >= 0
        x, y, z = p
        z = max(z, eps)
        # println("Query Upper")
        query(ds.upper, Point2{T}(x / z, y / z))
    else
        x, y, z = -p
        z = max(z, eps)
        # println("Query Lower")
        -query(ds.lower, Point2{T}(x / z, y / z))
    end
end