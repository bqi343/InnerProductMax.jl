using GeometryBasics
using LinearAlgebra

export vec_to_matrix

"""Returns: 3xn"""
function vec_to_matrix(x::Vector{Point3{T}}) where T
    reduce(hcat, x)
end

function cross(a::Point3{T}, b::Point3{T}, c::Point3{T}) where {T}
    LinearAlgebra.cross(b - a, c - a)
end

function collinear(a::Point3{T}, b::Point3{T}, c::Point3{T}) where {T}
    cross(a, b, c) == zeros(3)
end

function coplanar(a::Point3{T}, b::Point3{T}, c::Point3{T}, d::Point3{T}) where {T}
    dot(cross(a, b, c), d - a) == 0
end

"""Returns dimension of point set plus one. Expect this to be 4 in general case.

ponts: 3xn"""
function num_indep(points::Vector{Point3{T}}) where {T<:Real}
    n = length(points)
    if n == 0
        throw(DomainError("points cannot be empty"))
    end
    a = 1
    b = a + 1
    while b <= n && points[a] == points[b]
        b += 1
    end
    if b > n
        return 1
    end
    c = b + 1
    while c <= n && collinear(points[a], points[b], points[c])
        c += 1
    end
    if c > n
        return 2
    end
    d = c + 1
    while d <= n && coplanar(points[a], points[b], points[c], points[d])
        d += 1
    end
    if d > n
        return 3
    end
    return 4
end

# m: 3 x n
function num_indep(m::Matrix{T}) where {T<:Real}
    return num_indep(Point3.(eachcol(m)))
end