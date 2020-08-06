using DataInterpolations, Test

@testset "Interface" begin include("interface.jl") end
@testset "Interpolation Tests" begin include("interpolation_tests.jl") end
@testset "Derivative Tests" begin include("derivative_tests.jl") end
