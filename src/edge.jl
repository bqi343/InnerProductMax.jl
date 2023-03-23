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

function min_x(s::Segment{T}) where {T}
    s.l[1]
end

function max_x(s::Segment{T}) where {T}
    s.r[1]
end

function get_y_coord(s::Segment{T}, x::T) where {T}
    return ((x - s.l[1]) * s.r[2] + (s.r[1] - x) * s.l[2]) / (s.r[1] - s.l[1])
end

function Base.isless(a::Segment{T}, b::Segment{T}) where {T}
    l = max(min_x(a), min_x(b))
    r = min(max_x(a), max_x(b))
    @assert l < r "whoops $l $r $(min_x(a)) $(max_x(a)) $(min_x(b)) $(max_x(b)) $(a.l) $(a.r) $(b.l) $(b.r)"
    m = (l + r) / 2
    get_y_coord(a, m) < get_y_coord(b, m)
end

function is_above(s::Segment{T}, p::Point2{T}) where {T}
    get_y_coord(s, p[1]) > p[2]
end

struct AugmentedEdge{T}
    e::Edge{T}
    vert_below::Point3
    vert_above::Point3
end

min_x(a::AugmentedEdge{T}) where {T} = min_x(a.e)
max_x(a::AugmentedEdge{T}) where {T} = max_x(a.e)
Base.isless(p::Point2{T}, a::AugmentedEdge{T}) where {T} = is_above(a.e, p)
Base.isless(a::AugmentedEdge{T}, b::AugmentedEdge{T}) where {T} = a.e < b.e