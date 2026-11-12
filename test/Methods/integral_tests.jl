using DataInterpolations, Test
using QuadGK
using DataInterpolations: integral
using CurveFit, ForwardDiff
using StableRNGs
using Unitful

function test_integral(
        method; args = [], kwargs = [], name::String,
        test_cache_parameters::Bool = true
    )
    func = method(args...; kwargs..., extrapolation = ExtrapolationType.Extension)
    (; t) = func
    t1 = minimum(t)
    t2 = maximum(t)
    @testset "$name" begin
        # integral(A, t1, t2)
        qint, err = quadgk(func, t1, t2; atol = 1.0e-12, rtol = 1.0e-12)
        aint = integral(func, t1, t2)
        @test isapprox(qint, aint, atol = 1.0e-6, rtol = 1.0e-8)

        # integral(A, t)
        qint, err = quadgk(func, t1, (t1 + t2) / 2; atol = 1.0e-12, rtol = 1.0e-12)
        aint = integral(func, (t1 + t2) / 2)
        @test isapprox(qint, aint, atol = 1.0e-6, rtol = 1.0e-8)

        # integral(A, t1, t), integral(A, t, t2), integral(A, t, t)
        ts = range(t1, t2; length = 100)
        for t in ts
            qint, err = quadgk(func, t1, t; atol = 1.0e-12, rtol = 1.0e-12)
            aint1 = integral(func, t1, t)
            @test isapprox(qint, aint1, atol = 1.0e-5, rtol = 1.0e-8)
            aint2 = integral(func, t)
            @test aint1 == aint2

            qint, err = quadgk(func, t, t2; atol = 1.0e-12, rtol = 1.0e-12)
            aint = integral(func, t, t2)
            @test isapprox(qint, aint, atol = 1.0e-5, rtol = 1.0e-8)

            aint = integral(func, t, t)
            @test iszero(aint)
        end

        # integrals with extrapolation
        qint, err = quadgk(func, t1 - 5.0, (t1 + t2) / 2; atol = 1.0e-12, rtol = 1.0e-12)
        aint = integral(func, t1 - 5.0, (t1 + t2) / 2)
        @test isapprox(qint, aint, atol = 1.0e-6, rtol = 1.0e-8)

        qint, err = quadgk(func, (t1 + t2) / 2, t2 + 5.0; atol = 1.0e-12, rtol = 1.0e-12)
        aint = integral(func, (t1 + t2) / 2, t2 + 5.0)
        @test isapprox(qint, aint, atol = 1.0e-6, rtol = 1.0e-8)

        # Integrate intervals fully outside data
        qint, err = quadgk(func, t1 - 5.0, t1 - 2.5; atol = 1.0e-12, rtol = 1.0e-12)
        aint = integral(func, t1 - 5.0, t1 - 2.5)
        @test isapprox(qint, aint, atol = 1.0e-6, rtol = 1.0e-8)

        qint, err = quadgk(func, t2 + 2.5, t2 + 5.0; atol = 1.0e-12, rtol = 1.0e-12)
        aint = integral(func, t2 + 2.5, t2 + 5.0)
        @test isapprox(qint, aint, atol = 1.0e-6, rtol = 1.0e-8)
    end
    func = method(args...; kwargs...)
    @test_throws DataInterpolations.LeftExtrapolationError integral(func, t[1] - 1.0)
    @test_throws DataInterpolations.RightExtrapolationError integral(func, t[end] + 1.0)
    @test_throws DataInterpolations.LeftExtrapolationError integral(func, t[1] - 1.0, t[2])
    @test_throws DataInterpolations.RightExtrapolationError integral(
        func, t[1], t[end] + 1.0
    )

    return if test_cache_parameters
        # Test integration with cached parameters
        func = method(
            args...; kwargs..., cache_parameters = true,
            extrapolation = ExtrapolationType.Extension
        )
        qint, err = quadgk(func, t1 - 1, t1; atol = 1.0e-12, rtol = 1.0e-12)
        aint = integral(func, t1 - 1, t1)
        @test isapprox(qint, aint, atol = 1.0e-6, rtol = 1.0e-8)
    end
end

@testset "LinearInterpolation" begin
    u = 2.0collect(1:10)
    t = 1.0collect(1:10)
    test_integral(
        LinearInterpolation; args = [u, t], name = "Linear Interpolation (Vector)"
    )
    u = vcat((2.0collect(1:10))', (3.0collect(1:10))')
    test_integral(
        LinearInterpolation; args = [u, t], name = "Linear Interpolation (Matrix)"
    )
    u = [[2.0i, 3.0i] for i in 1:10]
    test_integral(
        LinearInterpolation; args = [u, t],
        name = "Linear Interpolation (Vector of Vectors)"
    )
    u = round.(rand(100), digits = 5)
    t = 1.0collect(1:100)
    test_integral(
        LinearInterpolation; args = [u, t],
        name = "Linear Interpolation (Vector) with random points"
    )
end

@testset "SmoothedConstantInterpolation" begin
    u = [1.0, 4.0, 9.0, 16.0]
    t = [1.0, 2.0, 3.0, 4.0]
    test_integral(
        SmoothedConstantInterpolation; args = [u, t], name = "Smoothed constant interpolation"
    )
    u2 = vcat(u', u')
    test_integral(
        SmoothedConstantInterpolation; args = [u2, t],
        name = "Smoothed constant interpolation (Matrix)"
    )
    u_vov = [[u_, u_] for u_ in u]
    test_integral(
        SmoothedConstantInterpolation; args = [u_vov, t],
        name = "Smoothed constant interpolation (Vector of Vectors)"
    )

    A_constant = ConstantInterpolation(u, t)
    I_ref = DataInterpolations.integral(A_constant, first(t), last(t))
    I_smoothed = [
        DataInterpolations.integral(
            SmoothedConstantInterpolation(u, t; d_max), first(t), last(t)
        )
            for d_max in 0.0:0.1:1.0
    ]
    @test all(I_smoothed .≈ I_ref)

    A = SmoothedConstantInterpolation(
        u, t; extrapolation = ExtrapolationType.Constant, d_max = 0
    )
    @test DataInterpolations.integral(A, 5.0, 6.0) == 16.0
end

@testset "ConstantInterpolation" begin
    u = [1.0, 4.0, 9.0, 16.0]
    t = [1.0, 2.0, 3.0, 4.0]
    test_integral(
        ConstantInterpolation; args = [u, t], name = "Constant Interpolation (Vector)"
    )
    u2 = vcat(u', u')
    test_integral(
        ConstantInterpolation; args = [u2, t], name = "Constant Interpolation (Matrix)"
    )
    u_vov = [[u_, u_] for u_ in u]
    test_integral(
        ConstantInterpolation; args = [u_vov, t],
        name = "Constant Interpolation (Vector of Vectors)"
    )
end

@testset "QuadraticInterpolation" begin
    u = [1.0, 4.0, 9.0, 16.0]
    t = [1.0, 2.0, 3.0, 4.0]
    test_integral(
        QuadraticInterpolation; args = [u, t], name = "Quadratic Interpolation (Vector)"
    )
    u2 = [1.0 4.0 9.0 16.0; 1.0 4.0 9.0 16.0]
    test_integral(
        QuadraticInterpolation; args = [u2, t], name = "Quadratic Interpolation (Matrix)"
    )
    u_vov = [[1.0, 1.0], [4.0, 4.0], [9.0, 9.0], [16.0, 16.0]]
    test_integral(
        QuadraticInterpolation; args = [u_vov, t],
        name = "Quadratic Interpolation (Vector of Vectors)"
    )
    u = [3.0, 0.0, 3.0, 0.0]
    t = [1.0, 2.0, 3.0, 4.0]
    test_integral(
        QuadraticInterpolation;
        args = [u, t, :Backward],
        name = "Quadratic Interpolation (Vector)"
    )
    u = round.(rand(100), digits = 5)
    t = 1.0collect(1:100)
    test_integral(
        QuadraticInterpolation; args = [u, t],
        name = "Quadratic Interpolation (Vector) with random points"
    )
end

@testset "LagrangeInterpolation" begin
    u = [1.0, 4.0, 9.0]
    t = [1.0, 2.0, 6.0]
    A = LagrangeInterpolation(u, t)
    @test_throws DataInterpolations.IntegralNotFoundError integral(A, 1.0, 2.0)
    @test_throws DataInterpolations.IntegralNotFoundError integral(A, 5.0)
end

@testset "QuadraticSpline" begin
    u = [0.0, 1.0, 3.0]
    t = [-1.0, 0.0, 1.0]
    test_integral(QuadraticSpline; args = [u, t], name = "Quadratic Spline (Vector)")
    u2 = [0.0 1.0 3.0; 0.0 1.0 3.0]
    test_integral(QuadraticSpline; args = [u2, t], name = "Quadratic Spline (Matrix)")
    u_vov = [[u_, u_] for u_ in u]
    test_integral(
        QuadraticSpline; args = [u_vov, t],
        name = "Quadratic Spline (Vector of Vectors)"
    )
    u = round.(rand(100), digits = 5)
    t = 1.0collect(1:100)
    test_integral(
        QuadraticSpline; args = [u, t], name = "Quadratic Spline (Vector) with random points"
    )
end

@testset "CubicSpline" begin
    u = [0.0, 1.0, 3.0]
    t = [-1.0, 0.0, 1.0]
    test_integral(CubicSpline; args = [u, t], name = "Cubic Spline (Vector)")
    u2 = [0.0 1.0 3.0; 0.0 1.0 3.0]
    test_integral(CubicSpline; args = [u2, t], name = "Cubic Spline (Matrix)")
    u_vov = [[u_, u_] for u_ in u]
    test_integral(
        CubicSpline; args = [u_vov, t], name = "Cubic Spline (Vector of Vectors)"
    )
    u = round.(rand(100), digits = 5)
    t = 1.0collect(1:100)
    test_integral(
        CubicSpline; args = [u, t], name = "Cubic Spline (Vector) with random points"
    )
end

@testset "AkimaInterpolation" begin
    u = [0.0, 2.0, 1.0, 3.0, 2.0, 6.0, 5.5, 5.5, 2.7, 5.1, 3.0]
    t = collect(0.0:10.0)
    test_integral(AkimaInterpolation; args = [u, t], name = "Akima Interpolation (Vector)")
    u2 = vcat(u', u')
    test_integral(
        AkimaInterpolation; args = [u2, t], name = "Akima Interpolation (Matrix)"
    )
    u_vov = [[u_, u_] for u_ in u]
    test_integral(
        AkimaInterpolation; args = [u_vov, t],
        name = "Akima Interpolation (Vector of Vectors)"
    )
    u = round.(rand(100), digits = 5)
    t = 1.0collect(1:100)
    test_integral(
        AkimaInterpolation; args = [u, t], name = "Akima Interpolation (Vector) with random points"
    )
end

@testset "CubicHermiteSpline" begin
    du = [-0.047, -0.058, 0.054, 0.012, -0.068, 0.0]
    u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
    t = [0.0, 62.25, 109.66, 162.66, 205.8, 252.3]
    test_integral(
        CubicHermiteSpline; args = [du, u, t],
        name = "Cubic Hermite Spline (Vector)"
    )
    du2 = vcat(du', du')
    u2 = vcat(u', u')
    test_integral(
        CubicHermiteSpline; args = [du2, u2, t],
        name = "Cubic Hermite Spline (Matrix)"
    )
    du_vov = [[du_, du_] for du_ in du]
    u_vov = [[u_, u_] for u_ in u]
    test_integral(
        CubicHermiteSpline; args = [du_vov, u_vov, t],
        name = "Cubic Hermite Spline (Vector of Vectors)"
    )

    u = round.(rand(100), digits = 5)
    t = 1.0collect(1:100)
    du = diff(u) ./ diff(t)
    push!(du, 0)
    test_integral(
        CubicHermiteSpline; args = [du, u, t],
        name = "Cubic Hermite Spline (Vector) with random points"
    )
end

@testset "QuinticHermiteSpline" begin
    ddu = [0.0, -0.00033, 0.0051, -0.0067, 0.0029, 0.0]
    du = [-0.047, -0.058, 0.054, 0.012, -0.068, 0.0]
    u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
    t = [0.0, 62.25, 109.66, 162.66, 205.8, 252.3]
    test_integral(
        QuinticHermiteSpline; args = [ddu, du, u, t],
        name = "Quintic Hermite Spline (Vector)"
    )
    ddu2 = vcat(ddu', ddu')
    du2 = vcat(du', du')
    u2 = vcat(u', u')
    test_integral(
        QuinticHermiteSpline; args = [ddu2, du2, u2, t],
        name = "Quintic Hermite Spline (Matrix)"
    )
    ddu_vov = [[ddu_, ddu_] for ddu_ in ddu]
    du_vov = [[du_, du_] for du_ in du]
    u_vov = [[u_, u_] for u_ in u]
    test_integral(
        QuinticHermiteSpline; args = [ddu_vov, du_vov, u_vov, t],
        name = "Quintic Hermite Spline (Vector of Vectors)"
    )

    u = round.(rand(100), digits = 5)
    t = 1.0collect(1:100)
    du = diff(u) ./ diff(t)
    push!(du, 0)
    ddu = diff(du) ./ diff(t)
    push!(ddu, 0)
    test_integral(
        QuinticHermiteSpline; args = [ddu, du, u, t],
        name = "Quintic Hermite Spline (Vector) with random points"
    )
end

@testset "Curvefit" begin
    rng = StableRNG(12345)
    model(x, p) = @. p[1] / (1 + exp(x - p[2]))
    t = range(-10, stop = 10, length = 40)
    u = model(t, [1.0, 2.0]) + 0.01 * randn(rng, length(t))
    p0 = [0.5, 0.5]
    A = Curvefit(u, t, model, p0)
    @test_throws DataInterpolations.IntegralNotFoundError integral(A, 0.0, 1.0)
    @test_throws DataInterpolations.IntegralNotFoundError integral(A, 5.0)
end

@testset "BSplineInterpolation" begin
    t = Float64[0, 62.25, 109.66, 162.66, 205.8, 252.3]
    u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
    test_integral(
        BSplineInterpolation;
        args = [u, t, 2, :Uniform],
        name = "BSpline Interpolation (d=2, Uniform)",
        test_cache_parameters = false
    )
    test_integral(
        BSplineInterpolation;
        args = [u, t, 3, :Average],
        name = "BSpline Interpolation (d=3, Average)",
        test_cache_parameters = false
    )
    u2 = vcat(u', u')
    test_integral(
        BSplineInterpolation;
        args = [u2, t, 2, :Uniform],
        name = "BSpline Interpolation (d=2, Uniform): Matrix",
        test_cache_parameters = false
    )
    u_vov = [[u_, u_] for u_ in u]
    test_integral(
        BSplineInterpolation;
        args = [u_vov, t, 2, :Uniform],
        name = "BSpline Interpolation (d=2, Uniform): Vector{Vector}",
        test_cache_parameters = false
    )
end

@testset "BSplineApprox" begin
    t = Float64[0, 62.25, 109.66, 162.66, 205.8, 252.3]
    u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
    test_integral(
        BSplineApprox;
        args = [u, t, 2, 4, :Uniform],
        name = "BSpline Approx (d=2, Uniform)",
        test_cache_parameters = false
    )
    test_integral(
        BSplineApprox;
        args = [u, t, 3, 4, :Average],
        name = "BSpline Approx (d=3, Average)",
        test_cache_parameters = false
    )
    u2 = vcat(u', u')
    test_integral(
        BSplineApprox;
        args = [u2, t, 2, 4, :Uniform],
        name = "BSpline Approx (d=2, Uniform): Matrix",
        test_cache_parameters = false
    )
    u_vov = [[u_, u_] for u_ in u]
    test_integral(
        BSplineApprox;
        args = [u_vov, t, 2, 4, :Uniform],
        name = "BSpline Approx (d=2, Uniform): Vector{Vector}",
        test_cache_parameters = false
    )
end

@testset "SmoothArcLengthInterpolation" begin
    u = [0.3 -1.5 3.1; -0.2 0.2 -1.5; 10.4 -37.2 -5.8]
    test_integral(
        SmoothArcLengthInterpolation;
        args = [u], kwargs = Pair[:m => 5],
        name = "Smooth Arc Length Interpolation"
    )
    u_vov = [u[:, i] for i in 1:size(u, 2)]
    test_integral(
        SmoothArcLengthInterpolation;
        args = [u_vov], kwargs = Pair[:m => 5],
        name = "Smooth Arc Length Interpolation (Vector of Vectors)"
    )
end

# issue #385
@testset "Integrals with unitful numbers" begin
    u = rand(5)u"m"
    A = ConstantInterpolation(u, (1:5)u"s")
    @test @inferred(integral(A, 4u"s")) ≈ sum(u[1:3]) * u"s"
end

@testset "cumulative_integral" begin
    A = ConstantInterpolation(["A", "B", "C"], [0.0, 0.25, 0.75])
    for cache_parameter in (true, false)
        @test @inferred(DataInterpolations.cumulative_integral(A, cache_parameter)) ===
            nothing
    end

    A = ConstantInterpolation([3.1, 2.5, 4.7], [0.0, 0.25, 0.75])
    @test @inferred(DataInterpolations.cumulative_integral(A, false)) == Float64[]
    @test @inferred(DataInterpolations.cumulative_integral(A, true)) == [0.775, 2.025]
end
