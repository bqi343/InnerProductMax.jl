module InnerProductMax

export AbstractInnerProductMax, query_all
abstract type AbstractInnerProductMax{T<:Real} end

include("geo_utils.jl")
include("distributions_3d.jl")
include("hull.jl")

"""Returns: 3xn"""
function query_all(ds::AbstractInnerProductMax{T}, queries::Matrix{T}) where {T}
    res = zero(queries)
    for (i, q) in enumerate(eachcol(queries))
        res[:, i] = query(ds, Point3{T}(q))
    end
    res
end

include("naive.jl")
include("nested.jl")
include("mine.jl")
include("plot.jl")

end # module InnerProductMax
