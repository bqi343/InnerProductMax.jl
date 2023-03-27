module InnerProductMaxTests

using InnerProductMax, ReTest
using GeometryBasics

const T = Float64
const P3 = Point3{T}

include("compare_chull.jl")
include("test_geo_utils.jl")
include("test_distributions_3d.jl")
include("test_hull.jl")
include("test_hull_lite.jl")
include("test_edge.jl")
include("test_persistent_treap.jl")
include("test_persistent_rb.jl")
include("test_correctness.jl")
include("test_perf.jl")

end