"""
Abstract type for persistently storing a set of unique keys in sorted order.
Supports insertion and deletion in O(log n) time each.
"""

abstract type BstNode{K} end
abstract type PointLocationDs{T<:Real} end

function query_pl end

include("edge.jl")
include("point_location_treap.jl")
include("point_location_rb.jl")

using InnerProductMax.PLTreap
using InnerProductMax.PLRB
export InnerProductMaxMine, PointLocationDsRB, PointLocationDsTreap