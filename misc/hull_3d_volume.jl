using Polyhedra: polyhedron, vrep, volume
using LinearAlgebra

points = [-1.0 1.0 1.0
    -1.0 -1.0 1.0
    1.0 0.0 -1.0
    1.0 1.0 1.0
    0.0 1.0 -1.0]

function volume_simplex(points)
    A = Matrix{T}(undef, 3, 3)
    for i in 1:3
        A[i, :] = points[i+1, :] - points[1, :]
    end
    return abs(det(A)) / 6
end

p = polyhedron(vrep(points))
println(volume(p)) # 5.3..., wrong
vol_1 = volume_simplex(points[[1, 2, 4, 5], :])
vol_2 = volume_simplex(points[[2, 4, 5, 3], :])
println("$vol_1 $vol_2 $(vol_1 + vol_2)") # 1.3..., 1.3..., 2.6..., as expected