export Hull, negated_copy

using QHull
const V_id = Int32

"""
ensures all faces point outwards, in-place
"""
function orient_faces!(hull::Chull)
    for (i, (a,b,c)) in enumerate(eachrow(hull.simplices))
        pa, pb, pc = Point3.([hull.points[a, :], hull.points[b, :], hull.points[c, :]])
        if dot(cross(pa, pb, pc), hull.facets[i, 1:3]) < 0
            hull.simplices[i, :] = reverse(hull.simplices[i, :])
        end
    end
end

struct Hull{T<:Real}
    points::Matrix{T}
    vertices::Vector{V_id}
    simplices::Matrix{V_id}
    facets::Matrix{T}
    adj_simplices::Matrix{V_id}
    Hull{T}(p::Matrix{T}, v::Vector{V_id}, s::Matrix{V_id}, f::Matrix{T}, a::Matrix{V_id}) where T = new{T}(p, v, s, f, a)
    """points: 3xn"""
    function Hull{T}(points::Matrix{T}) where {T}
        hull = chull(Matrix{T}(transpose(points)))
        orient_faces!(hull)
        d = Dict{Tuple{V_id,V_id},Int}()
        function associate(t, val)
            @assert !haskey(d, t)
            d[t] = val
        end
        for (i, (a, b, c)) in enumerate(eachrow(hull.simplices))
            associate((a, b), i)
            associate((b, c), i)
            associate((c, a), i)
        end
        adj = zeros(V_id, 3, size(hull.simplices, 1))
        for (i, (a, b, c)) in enumerate(eachrow(hull.simplices))
            adj[:, i] = [d[(b, a)], d[(c, b)], d[(a, c)]]
        end
        new{T}(transpose(hull.points), hull.vertices, transpose(hull.simplices), transpose(hull.facets), adj)
    end
end

"""
Returns a copy of hull with all points negated, and facets adjusted
accordingly.
"""
function negated_copy(hull::Hull{T}) where T
    facets = copy(hull.facets)
    facets[1:3, :] .*= -1
    Hull{T}(hull.points .* -1, hull.vertices, hull.simplices[[1, 3, 2], :], facets, hull.adj_simplices[[3, 2, 1], :])
end
