export Hull

using QHull
const V_id = Int32

"""
ensures all faces point outwards, in-place
"""
function orient_faces!(hull::Chull)
    for (i, (a, b, c)) in enumerate(eachrow(hull.simplices))
        pa, pb, pc = Point3.([hull.points[a, :], hull.points[b, :], hull.points[c, :]])
        if dot(cross(pa, pb, pc), hull.facets[i, 1:3]) < 0
            hull.simplices[i, :] = reverse(hull.simplices[i, :])
        end
    end
end

"""simplices: 3 x n"""
function construct_adj_simplices(simplices::AbstractMatrix{V_id})
    d = Dict{Tuple{V_id,V_id},Int}()
    function associate(t, val)
        @assert !haskey(d, t) "failed $simplices"
        d[t] = val
    end
    for (i, (a, b, c)) in enumerate(eachcol(simplices))
        associate((a, b), i)
        associate((b, c), i)
        associate((c, a), i)
    end
    n = size(simplices, 2)
    adj = zeros(V_id, 3, n)
    for (i, (a, b, c)) in enumerate(eachcol(simplices))
        adj[:, i] = [d[(b, a)], d[(c, b)], d[(a, c)]]
    end
    adj
end

function construct_adj_simplices_n2(n::Int, simplices::AbstractMatrix{V_id})
    d = zeros(V_id, n, n)
    function associate(x, y, val)
        @assert d[x, y] == 0 "failed $simplices"
        d[x, y] = val
    end
    for (i, (a, b, c)) in enumerate(eachcol(simplices))
        associate(a, b, i)
        associate(b, c, i)
        associate(c, a, i)
    end
    adj = zeros(V_id, 3, size(simplices, 2))
    for (i, (a, b, c)) in enumerate(eachcol(simplices))
        adj[:, i] = [d[b, a], d[c, b], d[a, c]]
    end
    adj
end

struct Hull{T<:Real}
    points::Matrix{T}
    vertices::Vector{V_id}
    simplices::Matrix{V_id}
    facets::Matrix{T}
    adj_simplices::Matrix{V_id}
    Hull(p::Matrix{T}, v::Vector{V_id}, s::Matrix{V_id}, f::Matrix{T}, a::Matrix{V_id}) where {T} = new{T}(p, v, s, f, a)
    """points: 3xn"""
    function Hull(points::AbstractMatrix{T}) where {T}
        hull = chull(Matrix{T}(points'))
        orient_faces!(hull)
        simplices = hull.simplices'
        new{T}(hull.points', hull.vertices, simplices, hull.facets', construct_adj_simplices(simplices))
    end
end

"""
Returns a copy of hull with all points negated, and facets adjusted
accordingly.
"""
function negated_copy(hull::Hull{T}) where {T}
    facets = copy(hull.facets)
    facets[1:3, :] .*= -1
    Hull(hull.points .* -1, hull.vertices, hull.simplices[[1, 3, 2], :], facets, hull.adj_simplices[[3, 2, 1], :])
end
