using DataInterpolations
using ForwardDiff
using Mooncake
using Test

# Differentiates a scalar accumulation sum over t w.r.t. u values.
# sum(A(_t)) handles both scalar u (A returns Number) and matrix u (A returns Vector).
function test_mooncake_ugrad(
        method, u, t; args = [], args_after = [], kwargs = [], name::String
    )
    trange = collect(range(minimum(t) - 5.0, maximum(t) + 5.0, step = 0.1))
    trange_exclude = filter(x -> !in(x, t), trange)
    return @testset "$name, derivatives w.r.t. u" begin
        function f(u)
            A = method(
                args..., u, t, args_after...; kwargs...,
                extrapolation = ExtrapolationType.Extension
            )
            out = zero(eltype(u))
            for _t in trange_exclude
                out += sum(A(_t))
            end
            out
        end
        cache = Mooncake.prepare_gradient_cache(f, u)
        _, (_, mgrad) = Mooncake.value_and_gradient!!(cache, f, u)
        fgrad = ForwardDiff.gradient(f, u)
        @test mgrad ≈ fgrad
    end
end

# Differentiates A(t) w.r.t. the scalar t and compares with DataInterpolations.derivative.
# Builds the rule once and reuses it across all test points.
# Subsamples to at most 50 points so CI time is bounded regardless of interpolation range.
function test_mooncake_tgrad(A, t; name::String)
    trange = collect(range(minimum(t) - 5.0, maximum(t) + 5.0, step = 0.1))
    trange_exclude = filter(x -> !in(x, t), trange)
    n = length(trange_exclude)
    stride = max(1, n ÷ 50)
    trange_test = trange_exclude[1:stride:end]
    return @testset "$name, derivatives w.r.t. t" begin
        f_t = _t -> A(_t)
        cache = Mooncake.prepare_gradient_cache(f_t, first(trange_test))
        for _t in trange_test
            _, (_, mgrad) = Mooncake.value_and_gradient!!(cache, f_t, _t)
            @test mgrad ≈ DataInterpolations.derivative(A, _t) atol = 1.0e-10
        end
    end
end

@testset "LinearInterpolation" begin
    u = vcat(collect(1.0:5.0), 2 * collect(6.0:10.0))
    t = collect(1.0:10.0)
    test_mooncake_ugrad(LinearInterpolation, u, t; name = "Linear Interpolation")
    A = LinearInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
    test_mooncake_tgrad(A, t; name = "Linear Interpolation")
    u_mat = Matrix(hcat(u, u)')
    test_mooncake_ugrad(
        LinearInterpolation, u_mat, t; name = "Linear Interpolation (matrix u)"
    )
end

@testset "QuadraticInterpolation" begin
    u = [1.0, 4.0, 9.0, 16.0]
    t = [1.0, 2.0, 3.0, 4.0]
    test_mooncake_ugrad(QuadraticInterpolation, u, t; name = "Quadratic Interpolation")
    A = QuadraticInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
    test_mooncake_tgrad(A, t; name = "Quadratic Interpolation")
    u_mat = Matrix(hcat(u, u)')
    test_mooncake_ugrad(
        QuadraticInterpolation, u_mat, t; name = "Quadratic Interpolation (matrix u)"
    )
end

@testset "ConstantInterpolation" begin
    u = [1.0, 2.0, 3.0, 4.0, 5.0]
    t = [0.0, 1.0, 2.0, 3.0, 4.0]
    test_mooncake_ugrad(ConstantInterpolation, u, t; name = "Constant Interpolation")
    A = ConstantInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
    test_mooncake_tgrad(A, t; name = "Constant Interpolation")
end

@testset "LagrangeInterpolation" begin
    u = [1.0, 4.0, 9.0, 16.0]
    t = [1.0, 2.0, 3.0, 4.0]
    test_mooncake_ugrad(LagrangeInterpolation, u, t; name = "Lagrange Interpolation")
    A = LagrangeInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
    test_mooncake_tgrad(A, t; name = "Lagrange Interpolation")
end

@testset "CubicSpline" begin
    u = [1.0, 4.0, 9.0, 16.0]
    t = [1.0, 2.0, 3.0, 4.0]
    test_mooncake_ugrad(CubicSpline, u, t; name = "Cubic Spline")
    A = CubicSpline(u, t; extrapolation = ExtrapolationType.Extension)
    test_mooncake_tgrad(A, t; name = "Cubic Spline")
end

@testset "AkimaInterpolation" begin
    # u-gradient excluded: abs() in slope computation makes ForwardDiff return NaN
    u = [1.0, 2.0, 3.0, 4.0, 5.0]
    t = [0.0, 1.0, 2.0, 3.0, 4.0]
    A = AkimaInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
    test_mooncake_tgrad(A, t; name = "Akima Interpolation")
end

@testset "CubicHermiteSpline" begin
    du = [-0.047, -0.058, 0.054, 0.012, -0.068, 0.0]
    u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
    t = [0.0, 62.25, 109.66, 162.66, 205.8, 252.3]
    test_mooncake_ugrad(CubicHermiteSpline, u, t; args = [du], name = "Cubic Hermite Spline")
    A = CubicHermiteSpline(du, u, t; extrapolation = ExtrapolationType.Extension)
    test_mooncake_tgrad(A, t; name = "Cubic Hermite Spline")
end

@testset "QuinticHermiteSpline" begin
    ddu = [0.0, -0.00033, 0.0051, -0.0067, 0.0029, 0.0]
    du = [-0.047, -0.058, 0.054, 0.012, -0.068, 0.0]
    u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
    t = [0.0, 62.25, 109.66, 162.66, 205.8, 252.3]
    test_mooncake_ugrad(
        QuinticHermiteSpline, u, t; args = [ddu, du], name = "Quintic Hermite Spline"
    )
    A = QuinticHermiteSpline(ddu, du, u, t; extrapolation = ExtrapolationType.Extension)
    test_mooncake_tgrad(A, t; name = "Quintic Hermite Spline")
end

@testset "BSplineInterpolation" begin
    u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
    t = [0.0, 62.25, 109.66, 162.66, 205.8, 252.3]
    test_mooncake_ugrad(
        BSplineInterpolation, u, t;
        args_after = [2, :Uniform, :Uniform], name = "BSpline Interpolation"
    )
    A = BSplineInterpolation(u, t, 2, :Uniform, :Uniform; extrapolation = ExtrapolationType.Extension)
    test_mooncake_tgrad(A, t; name = "BSpline Interpolation")
end

@testset "BSplineApprox" begin
    u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
    t = [0.0, 62.25, 109.66, 162.66, 205.8, 252.3]
    test_mooncake_ugrad(
        BSplineApprox, u, t;
        args_after = [2, 4, :Uniform, :Uniform], name = "BSpline Approximation"
    )
    A = BSplineApprox(u, t, 2, 4, :Uniform, :Uniform; extrapolation = ExtrapolationType.Extension)
    test_mooncake_tgrad(A, t; name = "BSpline Approximation")
end
