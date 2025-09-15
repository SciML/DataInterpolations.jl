using SafeTestsets, Test

const GROUP = get(ENV, "GROUP", "All")

@time begin
    if GROUP == "All" || GROUP == "Core"
        @testset "Core" begin
            @safetestset "Interface" include("interface.jl")
            @safetestset "Parameter Tests" include("parameter_tests.jl")
            @safetestset "Interpolation Tests" include("interpolation_tests.jl")
            @safetestset "Extrapolation Tests" include("extrapolation_tests.jl")
        end
    end

    if GROUP == "All" || GROUP == "Methods"
        @testset "Methods" begin
            @safetestset "Derivative Tests" include("derivative_tests.jl")
            @safetestset "Integral Tests" include("integral_tests.jl")
            @safetestset "Integral Inverse Tests" include("integral_inverse_tests.jl")
            @safetestset "Online Tests" include("online_tests.jl")
            @safetestset "Regularization Smoothing Tests" include("regularization.jl")
        end
    end

    if GROUP == "All" || GROUP == "Extensions"
        @testset "Extensions" begin
            @safetestset "SparseConnectivityTracer Tests" include("sparseconnectivitytracer_tests.jl")
            @safetestset "Zygote support Tests" include("zygote_tests.jl")
        end
    end

    if GROUP == "All" || GROUP == "Misc"
        @testset "Misc" begin
            @safetestset "Show methods Tests" include("show.jl")
        end
    end

    if GROUP == "QA"
        @safetestset "Quality Assurance" include("qa.jl")
    end
end
