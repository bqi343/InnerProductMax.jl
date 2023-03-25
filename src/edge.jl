"""A non-vertical line segment"""
abstract type Edge{T} end

struct Segment{T} <: Edge{T}
    l::Point2{T}
    r::Point2{T}
    function Segment{T}(l::Point2{T}, r::Point2{T}) where {T}
        if !(l[1] < r[1])
            throw(DomainError("cannot construct segment not pointing to left"))
        end
        new{T}(l, r)
    end
end

min_x(s::Segment) = s.l[1]
max_x(s::Segment) = s.r[1]
function get_y_coord(s::Segment{T}, x::T) where {T}
    ((x - s.l[1]) * s.r[2] + (s.r[1] - x) * s.l[2]) / (s.r[1] - s.l[1])
end

struct Ray{T} <: Edge{T}
    start::Point2{T}
    dir::Point2{T}
    function Ray{T}(start::Point2{T}, dir::Point2{T}) where {T}
        if dir[1] ≈ 0
            throw(DomainError("ray cannot be vertical"))
        end
        @assert norm(dir) ≈ 1 "expected dir to be normalized"
        new{T}(start, dir)
    end
end

min_x(s::Ray) =
    if s.dir[1] < 0
        -Inf
    else
        s.start[1]
    end
max_x(s::Ray) =
    if s.dir[1] > 0
        Inf
    else
        s.start[1]
    end

get_y_coord(r::Ray{T}, x::T) where {T} = r.start[2] + (x - r.start[1]) / r.dir[1] * r.dir[2]

struct AugmentedEdge{T}
    e::Edge{T}
    vert_below::Point3
    vert_above::Point3
end

min_x(a::AugmentedEdge{T}) where {T} = min_x(a.e)
max_x(a::AugmentedEdge{T}) where {T} = max_x(a.e)
Base.isless(p::Point2{T}, e::Edge{T}) where {T} = p[2] < get_y_coord(e, p[1])
Base.isless(p::Point2{T}, a::AugmentedEdge{T}) where {T} = Base.isless(p, a.e)
Base.isless(a::AugmentedEdge{T}, b::AugmentedEdge{T}) where {T} = a.e < b.e

function isless_at(a::Edge{T}, b::Edge{T}, x::T) where {T}
    if abs(x) == Inf
        throw(DomainError("x cannot be Inf"))
    end
    get_y_coord(a, x) < get_y_coord(b, x)
end

function Base.isless(a::Edge{T}, b::Edge{T}) where {T}
    l = max(min_x(a), min_x(b))
    r = min(max_x(a), max_x(b))
    @assert l < r "edge comparison failed: l = $l r = $r $a $b"
    if l == -Inf
        if a.dir ≈ b.dir
            isless_at(a, b, r)
        else
            a.dir[2] < b.dir[2]
        end
    elseif r == Inf
        if a.dir ≈ b.dir
            isless_at(a, b, l)
        else
            a.dir[2] < b.dir[2]
        end
    else
        m = (l + r) / 2
        isless_at(a, b, m)
    end
end