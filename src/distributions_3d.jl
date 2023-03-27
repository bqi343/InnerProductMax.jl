using LinearAlgebra
export tetrahedron, unit_sphere, cube

"""Returns: 3x5 matrix representing a tetrahedron with an additional point inside."""
tetrahedron(T::Type) = vec_to_matrix(Point3{T}[(0.5, 0.5, 0.5), (0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0)])

"""n real points uniformly on the unit sphere. note that Julia is column-major

Returns: 3xn"""
function unit_sphere(T::Type, n::Int)
    x = randn(T, 3, n)
    col_norms = norm.(eachcol(x))
    x ./ reshape(col_norms, (1, size(col_norms)...))
end

"""n integer points from the unit cube

Returns: 3xn"""
function cube(T::Type, d::Int, n::Int)
    if d <= 0
        throw(DomainError("d too small"))
    end
    if n < 4
        throw(DomainError("n too small"))
    end
    while true
        points = T.(rand(-d:d, (3, n)))
        if span_len(points) >= 4
            return points
        end
    end
end