using SafeTestsets

@safetestset "Quality Assurance" include("qa.jl")
@safetestset "Interface" include("interface.jl")
@safetestset "Parameter Tests" include("parameter_tests.jl")
@safetestset "Interpolation Tests" include("interpolation_tests.jl")
@safetestset "Derivative Tests" include("derivative_tests.jl")
@safetestset "Integral Tests" include("integral_tests.jl")
@safetestset "Integral Inverse Tests" include("integral_inverse_tests.jl")
@safetestset "Extrapolation Tests" include("extrapolation_tests.jl")
@safetestset "Online Tests" include("online_tests.jl")
@safetestset "Regularization Smoothing Tests" include("regularization.jl")
@safetestset "Show methods Tests" include("show.jl")
@safetestset "Zygote support Tests" include("zygote_tests.jl")
