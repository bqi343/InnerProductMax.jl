module InnerProductMaxTests

using InnerProductMax, ReTest

include("test_geo_utils.jl")
include("test_distributions_3d.jl")
include("test_hull.jl")
include("test_edge.jl")
include("test_persistent_treap.jl")
include("test_persistent_rb.jl")
include("test_correctness.jl")
include("test_perf.jl")

end