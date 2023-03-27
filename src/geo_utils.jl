using GeometryBasics
using LinearAlgebra

export vec_to_matrix


"""Returns: 3xn"""
vec_to_matrix(x::Vector{Point3{T}}) where {T} = reduce(hcat, x)
cross(a::Point3{T}, b::Point3{T}, c::Point3{T}) where {T} = LinearAlgebra.cross(b - a, c - a)
collinear(a::Point3{T}, b::Point3{T}, c::Point3{T}) where {T} = cross(a, b, c) ≈ zeros(T, 3)
dot_cross(a::Point3{T}, b::Point3{T}, c::Point3{T}, d::Point3{T}) where {T} = dot(cross(a, b, c), d - a)
coplanar(a::Point3{T}, b::Point3{T}, c::Point3{T}, d::Point3{T}) where {T} = dot_cross(a, b, c, d) ≈ 0

normalized(x::Point) = x ./ norm(x)

"""Returns dimension of point set plus one. Expect this to be 4 in general case.

ponts: 3xn"""
function compute_span(points::Vector{Point3{T}}) where {T<:Real}
    res = Point3{T}[]
    ids = V_id[]
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
    for (i, p) in enumerate(points)
        if !in_span(p)
            push!(ids, i)
            push!(res, p)
        end
    end
    res, ids
end

"""m: 3 x n"""
span_len(m::Vector{Point3{T}}) where {T} = length(compute_span(m)[1])
span_len(m::Matrix{T}) where {T<:Real} = span_len(Point3.(eachcol(m)))

"""Returns random unit vector perpendicular to a unit vector p"""
function rand_perp(p::Point3)
    dir = randn(3)
    dir = dir - dot(dir, p) * p
    normalized(dir)
end