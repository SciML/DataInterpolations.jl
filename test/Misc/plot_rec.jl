using DataInterpolations
using DataInterpolations: to_plottable, _plot_reshape

# `Plots` is deliberately not a test dependency (heavy, slow to precompile), and
# `RecipesBase.apply_recipe` requires a registered backend to run at all (it calls
# `is_key_supported`, which only a concrete plotting package provides). So these tests
# exercise `to_plottable`/`_plot_reshape` directly — the plain-Julia functions that
# actually reshape interpolation output into the `(nt, ndim)` orientation Plots expects —
# rather than the thin `@recipe`/`@series` declarations that wrap them.

@testset "to_plottable" begin
    t = [0.0, 1.0, 2.0, 3.0]

    @testset "Scalar (Vector)" begin
        A = LinearInterpolation([1.0, 3.0, 2.0, 5.0], t)
        plott, output = to_plottable(A; plotdensity = 50)
        @test length(plott) == 50
        @test output isa AbstractVector{<:Number}
        @test length(output) == 50
        @test output[1] ≈ A(plott[1])
        @test output[end] ≈ A(plott[end])
    end

    @testset "Matrix" begin
        A = LinearInterpolation([1.0 3.0 2.0 5.0; 2.0 1.0 4.0 3.0], t)
        plott, output = to_plottable(A; plotdensity = 50)
        @test size(output) == (50, 2)
        @test output[1, :] ≈ A(plott[1])
        @test output[end, :] ≈ A(plott[end])
    end

    @testset "Vector of Vectors" begin
        A = LinearInterpolation([[1.0, 2.0], [3.0, 1.0], [2.0, 4.0], [5.0, 3.0]], t)
        plott, output = to_plottable(A; plotdensity = 50)
        @test size(output) == (50, 2)
        @test output[1, :] ≈ A(plott[1])
        @test output[end, :] ≈ A(plott[end])
    end

    @testset "String data (non-Number scalar)" begin
        A = ConstantInterpolation(["A", "B", "C"], [0.0, 1.0, 2.0])
        plott, output = to_plottable(A; denseplot = false)
        @test plott == A.t
        @test output == A.u
    end

    @testset "SmoothArcLengthInterpolation (Matrix, (ndim, npoints) storage)" begin
        u = [0.3 -1.5 3.1; -0.2 0.2 -1.5; 10.4 -37.2 -5.8]
        A = SmoothArcLengthInterpolation(u; m = 5)
        plott, output = to_plottable(A; plotdensity = 50)
        @test size(output) == (50, 3)
        @test output[1, :] ≈ A(plott[1])
        @test output[end, :] ≈ A(plott[end])
    end

    @testset "Vector of Matrices (higher-rank per point)" begin
        u = [
            [1.0 2.0; 3.0 4.0], [5.0 6.0; 7.0 8.0],
            [9.0 10.0; 11.0 12.0], [1.0 1.0; 1.0 1.0],
        ]
        A = LinearInterpolation(u, t)
        plott, output = to_plottable(A; plotdensity = 50)
        @test size(output) == (50, 4)
    end
end

@testset "_plot_reshape for raw data points (scatter series)" begin
    t = [0.0, 1.0, 2.0, 3.0]

    A = LinearInterpolation([1.0 3.0 2.0 5.0; 2.0 1.0 4.0 3.0], t)
    reshaped = _plot_reshape(A.(A.t))
    @test size(reshaped) == (4, 2)
    @test reshaped ≈ [1.0 2.0; 3.0 1.0; 2.0 4.0; 5.0 3.0]

    u_vov = [[1.0, 2.0], [3.0, 1.0], [2.0, 4.0], [5.0, 3.0]]
    A_vov = LinearInterpolation(u_vov, t)
    @test _plot_reshape(A_vov.(A_vov.t)) ≈ reshaped
end
