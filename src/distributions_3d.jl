using LinearAlgebra
export unit_sphere, cube, tetrahedron

"""Returns: A tetrahedron with an additional point inside"""
function tetrahedron()
    vec_to_matrix(Point3{Float64}[(0.5, 0.5, 0.5), (0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0)])
end

"""
n real points uniformly on the unit sphere. note that julia is column-major
"""
function unit_sphere(n::Int)
    x = randn(3, n)
    col_norms = norm.(eachcol(x))
    x ./ reshape(col_norms, (1, size(col_norms)...))
end

"""
n integer points from the unit cube
"""
function cube(d::Int, n::Int)
    if d <= 0
        throw(DomainError("d too small"))
    end
    if n < 4
        throw(DomainError("n too small"))
    end
    while true
        points = Float64.(rand(-d:d, (3, n)))
        if num_indep(points) >= 4
            return points
        end
    end
end