# https://github.com/MakieOrg/Makie.jl/issues/797
using GLMakie, Polyhedra, LazySets
using LazySets: Hyperplane, translate
using InnerProductMax

export make_interactive_plot

"""Source: https://discourse.julialang.org/t/plotting-a-distribution-of-planes/53689/2"""
function get_plane(dir::Point3{T}, extreme_point::Point3{T}, scale::T) where {T}
    plane = Hyperplane(T.(dir), 0.0)
    # intersect with bounding box
    box = BallInf(zeros(3), scale)
    patch = intersection(plane, box)
    translate(patch, extreme_point)
end

"""ps: 3 x n"""
function make_interactive_plot(iprod_max::DataType, ps::AbstractMatrix{T}) where {T}
    # https://juliapolyhedra.github.io/Polyhedra.jl/stable/plot/#D-plotting-with-Plots-2
    m = Polyhedra.Mesh{3}(polyhedron(vrep(ps')))
    fig = Figure(resolution=(1000, 1000))
    ax = Axis3(fig[1, 1])
    # mesh!(m, color=:red)
    wireframe!(m, color=:red)
    sg = SliderGrid(fig[2, 1],
        (label="x", range=-1:0.01:1, startvalue=0.5),
        (label="y", range=-1:0.01:1, startvalue=0.2),
        (label="z", range=-1:0.01:1, startvalue=0.3)
    )

    h = Hull(ps)
    ds = iprod_max(h)

    dir = @lift Point3{T}($(sg.sliders[1].value), $(sg.sliders[2].value), $(sg.sliders[3].value))
    extreme_point = @lift query(ds, $dir)
    scale = 0.2
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
    fig
end

