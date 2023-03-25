using QHull
using GeometryBasics

export InnerProductMaxNaive, query

"""Preproc Time: O(N)
Preproc Mem: O(N)
Query Time: O(N)"""
struct InnerProductMaxNaive{T} <: AbstractInnerProductMax{T}
    points::Vector{Point3{T}}
    function InnerProductMaxNaive{T}(hull::Hull{T}) where {T}
        new{T}(Point3.(eachcol(hull.points[:, hull.vertices])))
    end
end

function query(ds::InnerProductMaxNaive{T}, q::Point3{T}) where {T<:Real}
    best = ds.points[1]
    best_val = dot(best, q)
    for k in ds.points
        cur_val = dot(k, q)
        if cur_val > best_val
            best = k
            best_val = cur_val
        end
    end
    best
end