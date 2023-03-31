# https://github.com/MakieOrg/Makie.jl/issues/797
using GLMakie
using Polyhedra: Mesh, polyhedron, vrep
using LazySets: Hyperplane, BallInf, intersection, translate, plot3d!
using InnerProductMax

export plot_point_location_interactive, plot_3d_hull_interactive

"""Source: https://discourse.julialang.org/t/plotting-a-distribution-of-planes/53689/2"""
function get_plane(dir::Point3{T}, extreme_point::Point3{T}, scale::T) where {T}
    plane = Hyperplane(T.(dir), 0.0)
    # intersect with bounding box
    box = BallInf(zeros(3), scale)
    patch = intersection(plane, box)
    translate(patch, extreme_point)
end

function get_range(points::AbstractMatrix, i)
    mn = minimum(points[i, :])
    mx = maximum(points[i, :])
    0.3 * (mn - mx) + mn, 0.3 * (mx - mn) + mx
end

"""points: 3 x n"""
function plot_3d_hull(points::AbstractMatrix{T}, grid, dir::Observable{Point3{T}}, extreme_point::Observable{Point3{T}}) where {T}
    # https://juliapolyhedra.github.io/Polyhedra.jl/stable/plot/#D-plotting-with-Plots-2
    m = Mesh{3}(polyhedron(vrep(points')))
    ax = Axis3(grid, aspect=(1, 1, 1))
    wireframe!(m, color=:red)

    scale = 0.2 # * diameter(points)
    plane = @lift get_plane($dir, $extreme_point, scale)
    prev_plane = plot3d!(to_value(plane), color=:azure)
    on(plane) do plane
        try
            GLMakie.delete!(ax.scene, prev_plane)
            prev_plane = nothing
            prev_plane = plot3d!(plane, color=:azure)
        catch e
        end
    end
    extreme_points = @lift [$extreme_point]
    dirs = @lift [$dir]
    arrows!(extreme_points, dirs, lengthscale=1.0 * scale, linewidth=0.05 * scale, arrowsize=Vec3f(0.2, 0.2, 0.3) * scale)
    scatter!(extreme_points, color=:blue, markersize=20)
    xlims!(get_range(points, 1)...)
    ylims!(get_range(points, 2)...)
    zlims!(get_range(points, 3)...)
end

"""Creates hull visualization and sliders."""
function plot_3d_hull_interactive(iprod_max::DataType, points::AbstractMatrix{T}) where {T}
    fig = Figure(resolution=(1000, 1000))
    fig[0, :] = Label(fig, "Inner Product Maximization", fontsize=50)
    sg = SliderGrid(fig[2, 1],
        (label="qx", range=-1:0.1:1, startvalue=0.5),
        (label="qy", range=-1:0.1:1, startvalue=0.2),
        (label="qz", range=-1:0.1:1, startvalue=0.3)
    )
    dir = @lift Point3{T}($(sg.sliders[1].value), $(sg.sliders[2].value), $(sg.sliders[3].value))
    h = Hull(points)
    ds = iprod_max(h)
    extreme_point = @lift query(ds, $dir)
    plot_3d_hull(points, fig[1, 1], dir, extreme_point)
    fig
end

function points_from_segments(segs::Vector{Tuple{Point2{T},Point2{T}}}) where {T}
    prv = Dict{Point2{T},Point2{T}}()
    nxt = Dict{Point2{T},Point2{T}}()
    for (a, b) in segs
        @assert !haskey(nxt, a) && !haskey(prv, b)
        nxt[a] = b
        prv[b] = a
    end
    start_points = []
    for a in keys(nxt)
        if !haskey(prv, a)
            push!(start_points, a)
        end
    end
    @assert length(start_points) <= 1
    if length(start_points) == 0
        for a in keys(nxt)
            push!(start_points, a)
            break
        end
    end
    res = [start_points[1]]
    while true
        if !haskey(nxt, res[end])
            break
        end
        next_point = nxt[res[end]]
        if next_point == res[1]
            break
        end
        push!(res, next_point)
    end
    res
end

function plot_point_location_interactive(iprod_max::DataType, points::AbstractMatrix{T}) where {T}
    h = Hull(points)
    ds = iprod_max(h)

    segments = Tuple{Point2{T},Point2{T}}[]
    for aug_edge in ds.upper.all_augmented_edges
        edge = aug_edge.e
        push!(segments, endpoints(edge))
    end
    polygon_segments = Dict{Point3{T},Vector{Tuple{Point2{T},Point2{T}}}}()
    for aug_edge in ds.upper.all_augmented_edges
        edge = aug_edge.e
        ends = endpoints(edge)
        # println("ends = ", ends)
        # println(aug_edge.vert_above)
        # println(aug_edge.vert_below)
        if !haskey(polygon_segments, aug_edge.vert_above)
            polygon_segments[aug_edge.vert_above] = Vector{Tuple{Point2{T},Point2{T}}}()
        end
        push!(polygon_segments[aug_edge.vert_above], ends)
        if !haskey(polygon_segments, aug_edge.vert_below)
            polygon_segments[aug_edge.vert_below] = Vector{Tuple{Point2{T},Point2{T}}}()
        end
        push!(polygon_segments[aug_edge.vert_below], (ends[2], ends[1]))
    end
    polygons = Dict{Point3{T},Vector{Point2{T}}}(k => points_from_segments(v) for (k, v) in polygon_segments)

    # println("segments = ", segments)
    # println("polygon_segments = ", polygon_segments)
    # println("polygons = ", polygons)

    fig = Figure(resolution=(2000, 1000))
    fig[0, 1:2] = Label(fig, "Inner Product Maximization via Point Location", fontsize=50)
    sg = SliderGrid(fig[2, 1:2],
        (label="qx", range=-2:0.1:2, startvalue=0.5),
        (label="qy", range=-2:0.1:2, startvalue=0.2),
    )
    ax = Axis(fig[1, 2], aspect=1, xlabel="qx", ylabel="qy")

    dir = @lift Point3{T}($(sg.sliders[1].value), $(sg.sliders[2].value), 1)
    dirs = @lift [Point2{T}(($dir)[1], ($dir)[2])]
    extreme_point = @lift query(ds, $dir)
    polygon = @lift polygons[$extreme_point]
    # println("polygon = ", polygon)
    # println(typeof(polygon))
    # poly!(to_value(polygon), color=:blue)
    # p = Point2f[(0.75, 0.25), (1.75, 0.25), (2.25, 0.75), (1.25, 0.75)]
    # println(typeof(p)) # Vector{GeometryBasics.Point{2, Float32}}
    xlims!(ax, -2, 2)
    ylims!(ax, -2, 2)
    poly!(polygon, color=:lightblue)
    linesegments!(segments, color=:red, linewidth=5)
    scatter!(dirs, markersize=20)
    plot_3d_hull(points, fig[1, 1], dir, extreme_point)
    fig
end