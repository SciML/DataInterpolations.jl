using SafeTestsets

@safetestset "Quality Assurance" include("qa.jl")
@safetestset "Interface" include("interface.jl")
@safetestset "Parameter Tests" include("parameter_tests.jl")
@safetestset "Interpolation Tests" include("interpolation_tests.jl")
@safetestset "Derivative Tests" include("derivative_tests.jl")
@safetestset "Integral Tests" include("integral_tests.jl")
@safetestset "Online Tests" include("online_tests.jl")
@safetestset "Regularization Smoothing" include("regularization.jl")
@safetestset "Show methods" include("show.jl")
