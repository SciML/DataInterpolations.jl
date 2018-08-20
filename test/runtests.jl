using DataInterpolations, Test

@testset "Interface" begin include("interface.jl") end
@testset "Linear Interpolation" begin include("linear.jl") end
