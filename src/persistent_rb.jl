module PersistentRB

using InnerProductMax: BstNode

export RedBlackNode, RBNode, insert_root, delete_root, tour, lower_bound, get_colors, get_last

operations = 0
rotations = 0
allocations = 0

"""
O(log n)-time, O(1) amortized memory partially-persistent BBST
References:
- https://dl.acm.org/doi/pdf/10.1145/6138.6151
- https://github.com/JuliaCollections/DataStructures.jl/blob/master/src/red_black_tree.jl
"""
mutable struct RedBlackNode{T<:Real,K} <: BstNode{K}
    const data::K
    const c::Tuple{Union{RedBlackNode{T,K},Nothing},Union{RedBlackNode{T,K},Nothing}}
    red::Bool # only matters on latest nodes
    # extra pointer for persistence
    time::T
    extra_dir::Int
    extra_c::Union{RedBlackNode{T,K},Nothing}
    function RedBlackNode{T,K}(d::K, c::Tuple{Union{RedBlackNode{T,K},Nothing},Union{RedBlackNode{T,K},Nothing}}, red::Bool) where {T,K}
        global allocations += 1
        if allocations % 100000 == 0
            println("ops = $operations rotations = $rotations allocations = $allocations")
        end
        new{T,K}(d, c, red, Inf, 0, nothing)
    end
    RedBlackNode{T,K}(d::K) where {T,K} = RedBlackNode{T,K}(d, (nothing, nothing), true)
end

const RBNode = Union{RedBlackNode{T,K},Nothing} where {T,K}

# getters
"""return whether t is red. null nodes are black."""
is_red(t::RBNode) = t !== nothing && t.red
"""return child in specified direction and time"""
function get_c(t::RedBlackNode{T,K}, dir::Int, time::T) where {T,K}
    # @assert dir in [1, 2] "dir must be 1 or 2"
    if time >= t.time && t.extra_dir == dir
        return t.extra_c
    end
    t.c[dir]
end
"""return latest child in specified direction"""
function get_c_latest(t::RedBlackNode{T,K}, dir::Int) where {T,K}
    if t.extra_dir == dir
        t.extra_c
    else
        t.c[dir]
    end
end
"""returns children as list"""
function get_children_latest(t::RedBlackNode{T,K}) where {T,K}
    RBNode{T,K}[get_c_latest(t, 1), get_c_latest(t, 2)]
end

function any_red_children(t::RedBlackNode{T,K}) where {T,K}
    any(is_red(get_c_latest(t, i)) for i in 1:2)
end

"""returns t or a copy of t with the specified child"""
function with_child(t::RedBlackNode{T,K}, dir::Int, time::T, c::RBNode{T,K}) where {T,K}
    if get_c_latest(t, dir) == c # no change
        return t
    end
    if t.extra_dir == 0 # use the extra pointer
        t.time = time
        t.extra_dir = dir
        t.extra_c = c
        return t
    end
    new_c = get_children_latest(t)
    new_c[dir] = c
    RedBlackNode{T,K}(t.data, Tuple(new_c), t.red)
end

"""returns a copy of t with new data"""
function with_new_data(t::RedBlackNode{T,K}, new_data::K) where {T,K}
    RedBlackNode{T,K}(new_data, Tuple(get_children_latest(t)), t.red)
end

"""persistently rotate up child of t in direction dir"""
function rotate(t::RedBlackNode{T,K}, dir::Int, time::T) where {T,K}
    global rotations += 1
    a = get_c_latest(t, dir)
    b = get_c_latest(a, dir ⊻ 3)
    t = with_child(t, dir, time, b)
    with_child(a, dir ⊻ 3, time, t)
end

"""figure 4"""
function insert_int(t::RBNode{T,K}, time::T, d::K) where {T,K}
    if t === nothing
        return RedBlackNode{T,K}(d)
    end
    dir = if d < t.data
        1
    else
        2
    end
    c = insert_int(get_c_latest(t, dir), time, d)
    # @assert c !== nothing
    t = with_child(t, dir, time, c)
    # @assert get_c_latest(t, dir) === c
    if c.red
        if any_red_children(c) # found two consecutive red nodes, need to fix
            # @assert !t.red
            if is_red(get_c_latest(t, dir ⊻ 3))
                # case (a): only change colors
                for i in 1:2
                    get_c_latest(t, i).red = false
                end
                t.red = true
                return t
            end
            # @assert !is_red(get_c_latest(c, 1)) || !is_red(get_c_latest(c, 2)) "children can't both be red"
            if is_red(get_c_latest(c, dir ⊻ 3)) # case (d): rotate
                root = rotate(c, dir ⊻ 3, time)
                t = with_child(t, dir, time, root)
            end
            # case (c), (d): change colors and rotate
            root = rotate(t, dir, time)
            root.red = false
            get_c_latest(root, dir ⊻ 3).red = true
            return root
        end
    end
    return t
end

function insert_root(t::RBNode{T,K}, time::T, d::K) where {T,K}
    global operations += 1
    t = insert_int(t, time, d)
    t.red = false # figure 4 case (b): root is always black
    t
end

"""Figure 5 cases (c), (d), (e). (TODO: check)"""
function delete_fixup_final(t::RedBlackNode{T,K}, dir::Int, time::T) where {T,K}
    # @assert !is_red(get_c_latest(t, dir))
    c_other = get_c_latest(t, dir ⊻ 3)
    # @assert other_child !== nothing
    root_red = t.red
    if root_red && !any_red_children(c_other) # case (c)
        # @assert !is_red(other_child)
        t.red = false
        c_other.red = true
        return t
    end
    # remaining: (d) or (e)
    # @assert !other_child.red
    # @assert any_red_children(other_child)
    if !is_red(get_c_latest(c_other, dir ⊻ 3)) # case (e)
        # @assert is_red(get_c_latest(other_child, dir))
        t = with_child(t, dir ⊻ 3, time, rotate(c_other, dir, time))
    end
    root = rotate(t, dir ⊻ 3, time) # cases (d), (e)
    root.red = root_red
    for i in 1:2
        get_c_latest(root, i).red = false
    end
    root
end

"""Figure 5. dir is direction of short node"""
function delete_fixup(t::RedBlackNode{T,K}, dir::Int, time::T) where {T,K}
    if !is_red(t)
        other_child = get_c_latest(t, dir ⊻ 3)
        if !is_red(other_child)
            if !any_red_children(other_child) # case (a): bubble up
                other_child.red = true
                return t, true
            end
        else # case (b)
            root = rotate(t, dir ⊻ 3, time)
            root.red = false
            get_c_latest(root, dir).red = true
            new_c = delete_fixup_final(get_c_latest(root, dir), dir, time)
            return with_child(root, dir, time, new_c), false
        end
    end
    delete_fixup_final(t, dir, time), false
end

function delete_rightmost(t::RedBlackNode{T,K}, time::T) where {T,K}
    dir = 2
    c = get_c_latest(t, dir)
    if c === nothing
        c_other = get_c_latest(t, dir ⊻ 3)
        if is_red(c_other) # not explicitly described?
            c_other.red = false
            return (c_other, false, t.data)
        end
        return (c_other, !t.red, t.data)
    end
    (c, is_short, new_data) = delete_rightmost(c, time)
    t = with_child(t, dir, time, c)
    if !is_short
        return (t, false, new_data)
    end
    delete_fixup(t, dir, time)..., new_data
end

function delete_int(t::RedBlackNode{T,K}, time::T, d::K) where {T,K}
    dir = if d < t.data
        1
    elseif t.data < d
        2
    else
        0
    end
    # println("delete_int ", dir)
    # println("delete_int $(t.data) $(d) $(dir)")
    if dir == 0 # data is at current node
        c = get_c_latest(t, 1)
        if c === nothing
            c_other = get_c_latest(t, 2)
            if is_red(c_other) # not explicitly described?
                c_other.red = false
                return (c_other, false)
            end
            return (c_other, !t.red)
        end
        dir = 1
        (c, is_short, new_data) = delete_rightmost(c, time)
        t = with_new_data(t, new_data) # swap data
    else
        (c, is_short) = delete_int(get_c_latest(t, dir), time, d)
    end
    t = with_child(t, dir, time, c)
    if !is_short
        return (t, false)
    end
    delete_fixup(t, dir, time)
end

function delete_root(t::RedBlackNode{T,K}, time::T, d::K) where {T,K}
    global operations += 1
    res = delete_int(t, time, d)
    # make root red
    if res[1] !== nothing
        res[1].red = false
    end
    res[1]
end

function lower_bound(t::RBNode, p)
    ret = nothing
    while t !== nothing
        if p < t.data
            ret = t.data
            t = get_c(t, 1, p[1])
        else
            t = get_c(t, 2, p[1])
        end
    end
    ret
end

function get_last(t::RedBlackNode{T,K}, time::T)::K where {T,K}
    while true
        r = get_c(t, 2, time)
        if r === nothing
            return t.data
        end
        t = r
    end
end

"""
Pushes all elements of the binary tree to v
"""
function tour(t::RBNode, time::T, v::Vector{K}, check::Bool) where {K,T}
    if t === nothing
        return 0
    end
    if check
        if t.red # no two consecutive red nodes
            @assert !any_red_children(t) "found two consecutive red children"
        end
    end
    path_1 = tour(get_c(t, 1, time), time, v, check)
    if !isempty(v)
        @assert v[end] < t.data
    end
    push!(v, t.data)
    path_2 = tour(get_c(t, 2, time), time, v, check)
    if check
        @assert path_1 == path_2 "$path_1 != $path_2, but paths must have equal number of black nodes"
    end
    path_2 + !t.red
end

"""
Returns result of traversing BST from smallest to largest
"""
function tour(t::RedBlackNode{T,K}, time::T, check::Bool=false) where {T,K}
    v = K[]
    @assert !t.red
    tour(t, time, v, check)
    v
end

function get_colors(t::RBNode)
    if t === nothing
        return 'x'
    end
    # $(t.data)
    "($(if t.red 'R' else 'B' end) $(get_colors(get_c_latest(t, 1))) $(get_colors(get_c_latest(t, 2))))"
end

end