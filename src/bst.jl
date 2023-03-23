export tour

"""
Abstract type for persistently storing a set of unique keys in sorted order.
Supports insertion and deletion in O(log n) time each.
"""
abstract type BstNode{K} end

"""
Pushes all elements of the binary tree to v
"""
function tour(t::Union{BstNode{K},Nothing}, v::Vector{K}) where {K}
    if t === nothing
        return
    end
    tour(t.l, v)
    if !isempty(v)
        @assert v[end] < t.data
    end
    push!(v, t.data)
    tour(t.r, v)
end

"""
Returns result of traversing BST from smallest to largest
"""
function tour(t::BstNode{K}) where {K}
    v = K[]
    tour(t, v)
    v
end