using Test
using NURBS

@testset "NURBS.jl" begin
    include("test_basis.jl")
    include("test_curves.jl")
    include("test_surfaces.jl")
    include("test_rational.jl")
    include("test_knots.jl")
    include("test_degree.jl")
    include("test_advanced.jl")
    include("test_conics.jl")
    include("test_fitting.jl")
    include("test_modification.jl")
end
