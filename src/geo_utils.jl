using GeometryBasics
using LinearAlgebra

export vec_to_matrix


"""Returns: 3xn"""
vec_to_matrix(x::Vector{Point3{T}}) where {T} = reduce(hcat, x)

cross(a::Point3{T}, b::Point3{T}, c::Point3{T}) where {T} = LinearAlgebra.cross(b - a, c - a)

collinear(a::Point3{T}, b::Point3{T}, c::Point3{T}) where {T} = cross(a, b, c) ≈ zeros(3)

coplanar(a::Point3{T}, b::Point3{T}, c::Point3{T}, d::Point3{T}) where {T} = dot(cross(a, b, c), d - a) ≈ 0

normalized(x::Point) = x ./ norm(x)

"""Returns dimension of point set plus one. Expect this to be 4 in general case.

ponts: 3xn"""
function compute_span(points::Vector{Point3{T}}) where {T<:Real}
    res = []
    function in_span(p)
        if length(res) == 0
            false
        elseif length(res) == 1
            res[1] ≈ p
        elseif length(res) == 2
            collinear(res[1], res[2], p)
        elseif length(res) == 3
            coplanar(res[1], res[2], res[3], p)
        else
            true
        end
    end
    for p in points
        if !in_span(p)
            push!(res, p)
        end
    end
    res
end

"""m: 3 x n"""
compute_span(m::Matrix{T}) where {T<:Real} = compute_span(Point3.(eachcol(m)))

"""Returns random unit vector perpendicular to a unit vector p"""
function rand_perp(p::Point3)
    dir = randn(3)
    dir = dir - dot(dir, p) * p
    normalized(dir)
end