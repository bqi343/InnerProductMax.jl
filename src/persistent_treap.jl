"""
O(log n) fully-persistent BST
References:
- https://en.wikipedia.org/wiki/Treap
- https://github.com/JuliaCollections/DataStructures.jl/blob/master/src/avl_tree.jl
- https://github.com/yurivish/Treaps.jl/blob/master/src/Treaps.jl
"""

struct TreapNode{K} <: BstNode{K}
    priority::Int32
    data::K
    l::Union{TreapNode{K},Nothing}
    r::Union{TreapNode{K},Nothing}
    TreapNode{K}(d::K) where {K} = new{K}(rand(Int32), d, nothing, nothing)
    TreapNode{K}(priority::Int32, data::K, l::Union{TreapNode{K},Nothing}, r::Union{TreapNode{K},Nothing}) where {K} = new{K}(priority, data, l, r)
end

const TNode = Union{TreapNode{K},Nothing} where {K};

"""Initialize new treap node from old one."""
function with_children(tnode::TreapNode{K}, l::TNode, r::TNode) where {K}
    TreapNode{K}(tnode.priority, tnode.data, l, r)
end

"""
splits t such that values >= d go to right
"""
function split(t::TNode, d::K)::Tuple{TNode,TNode} where {K}
    if t === nothing
        return t, t
    end
    if t.data < d
        r_l, r_r = split(t.r, d)
        with_children(t, t.l, r_l), r_r
    else
        l_l, l_r = split(t.l, d)
        l_l, with_children(t, l_r, t.r)
    end
end

function merge(l::TNode{K}, r::TNode{K}) where {K}
    if l === nothing
        return r
    end
    if r === nothing
        return l
    end
    if l.priority > r.priority
        with_children(l, l.l, merge(l.r, r))
    else
        with_children(r, merge(l, r.l), r.r)
    end
end

function delete(t::TNode, d::K) where {K}
    if t === nothing
        throw(DomainError("key not in treap"))
    end
    if t.data < d
        with_children(t, t.l, delete(t.r, d))
    elseif d < t.data
        with_children(t, delete(t.l, d), t.r)
    else
        merge(t.l, t.r)
    end
end

function insert(t::TNode, d::K) where {K}
    l, r = split(t, d)
    merge(merge(l, TreapNode{K}(d)), r)
end

"""returns first data greater than p"""
function lower_bound(t::TNode, p)
    if t === nothing
        return nothing
    end
    if is_above(t.data, p)
        ret = lower_bound(t.l, p)
        if ret === nothing
            t.data
        else
            ret
        end
    else
        lower_bound(t.r, p)
    end
end

function last(t::TreapNode{K})::K where {K}
    if t.r === nothing
        t.data
    else
        last(t.r)
    end
end