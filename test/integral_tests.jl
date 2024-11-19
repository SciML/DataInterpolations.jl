using DataInterpolations, Test
using QuadGK
using DataInterpolations: integral
using Optim, ForwardDiff
using RegularizationTools
using StableRNGs

function test_integral(method; args = [], kwargs = [], name::String)
    func = method(args...; kwargs..., extrapolation_left = ExtrapolationType.extension,
        extrapolation_right = ExtrapolationType.extension)
    (; t) = func
    t1 = minimum(t)
    t2 = maximum(t)
    @testset "$name" begin
        # integral(A, t1, t2)
        qint, err = quadgk(func, t1, t2; atol = 1e-12, rtol = 1e-12)
        aint = integral(func, t1, t2)
        @test isapprox(qint, aint, atol = 1e-6, rtol = 1e-8)

        # integral(A, t)
        qint, err = quadgk(func, t1, (t1 + t2) / 2; atol = 1e-12, rtol = 1e-12)
        aint = integral(func, (t1 + t2) / 2)
        @test isapprox(qint, aint, atol = 1e-6, rtol = 1e-8)

        # integral(A, t1, t), integral(A, t, t2), integral(A, t, t)
        ts = range(t1, t2; length = 100)
        for t in ts
            qint, err = quadgk(func, t1, t; atol = 1e-12, rtol = 1e-12)
            aint1 = integral(func, t1, t)
            @test isapprox(qint, aint1, atol = 1e-5, rtol = 1e-8)
            aint2 = integral(func, t)
            @test aint1 == aint2

            qint, err = quadgk(func, t, t2; atol = 1e-12, rtol = 1e-12)
            aint = integral(func, t, t2)
            @test isapprox(qint, aint, atol = 1e-5, rtol = 1e-8)

            aint = integral(func, t, t)
            @test aint == 0.0
        end

        # integrals with extrapolation
        qint, err = quadgk(func, t1 - 5.0, (t1 + t2) / 2; atol = 1e-12, rtol = 1e-12)
        aint = integral(func, t1 - 5.0, (t1 + t2) / 2)
        @test isapprox(qint, aint, atol = 1e-6, rtol = 1e-8)

        qint, err = quadgk(func, (t1 + t2) / 2, t2 + 5.0; atol = 1e-12, rtol = 1e-12)
        aint = integral(func, (t1 + t2) / 2, t2 + 5.0)
        @test isapprox(qint, aint, atol = 1e-6, rtol = 1e-8)
    end
    func = method(args...; kwargs...)
    @test_throws DataInterpolations.DownExtrapolationError integral(func, t[1] - 1.0)
    @test_throws DataInterpolations.UpExtrapolationError integral(func, t[end] + 1.0)
    @test_throws DataInterpolations.DownExtrapolationError integral(func, t[1] - 1.0, t[2])
    @test_throws DataInterpolations.UpExtrapolationError integral(func, t[1], t[end] + 1.0)
end

@testset "LinearInterpolation" begin
    u = 2.0collect(1:10)
    t = 1.0collect(1:10)
    test_integral(
        LinearInterpolation; args = [u, t], name = "Linear Interpolation (Vector)")
    u = round.(rand(100), digits = 5)
    t = 1.0collect(1:100)
    test_integral(LinearInterpolation; args = [u, t],
        name = "Linear Interpolation (Vector) with random points")
end

@testset "QuadraticInterpolation" begin
    u = [1.0, 4.0, 9.0, 16.0]
    t = [1.0, 2.0, 3.0, 4.0]
    test_integral(
        QuadraticInterpolation; args = [u, t], name = "Quadratic Interpolation (Vector)")
    u = [3.0, 0.0, 3.0, 0.0]
    t = [1.0, 2.0, 3.0, 4.0]
    test_integral(QuadraticInterpolation;
        args = [u, t, :Backward],
        name = "Quadratic Interpolation (Vector)")
    u = round.(rand(100), digits = 5)
    t = 1.0collect(1:10)
    test_integral(QuadraticInterpolation; args = [u, t],
        name = "Quadratic Interpolation (Vector) with random points")
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
    u = round.(rand(100), digits = 5)
    t = 1.0collect(1:100)
    test_integral(
        QuadraticSpline; args = [u, t], name = "Quadratic Spline (Vector) with random points")
end

@testset "CubicSpline" begin
    u = [0.0, 1.0, 3.0]
    t = [-1.0, 0.0, 1.0]
    test_integral(CubicSpline; args = [u, t], name = "Cubic Spline (Vector)")
    u = round.(rand(100), digits = 5)
    t = 1.0collect(1:100)
    test_integral(
        CubicSpline; args = [u, t], name = "Cubic Spline (Vector) with random points")
end

@testset "AkimaInterpolation" begin
    u = [0.0, 2.0, 1.0, 3.0, 2.0, 6.0, 5.5, 5.5, 2.7, 5.1, 3.0]
    t = collect(0.0:10.0)
    test_integral(AkimaInterpolation; args = [u, t], name = "Akima Interpolation (Vector)")
    u = round.(rand(100), digits = 5)
    t = 1.0collect(1:100)
    test_integral(
        AkimaInterpolation; args = [u, t], name = "Akima Interpolation (Vector) with random points")
end

@testset "CubicHermiteSpline" begin
    du = [-0.047, -0.058, 0.054, 0.012, -0.068, 0.0]
    u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
    t = [0.0, 62.25, 109.66, 162.66, 205.8, 252.3]
    test_integral(CubicHermiteSpline; args = [du, u, t],
        name = "Cubic Hermite Spline (Vector)")

    u = round.(rand(100), digits = 5)
    t = 1.0collect(1:100)
    du = diff(u) ./ diff(t)
    push!(du, 0)
    test_integral(CubicHermiteSpline; args = [du, u, t],
        name = "Cubic Hermite Spline (Vector) with random points")
end

@testset "QuinticHermiteSpline" begin
    ddu = [0.0, -0.00033, 0.0051, -0.0067, 0.0029, 0.0]
    du = [-0.047, -0.058, 0.054, 0.012, -0.068, 0.0]
    u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
    t = [0.0, 62.25, 109.66, 162.66, 205.8, 252.3]
    test_integral(QuinticHermiteSpline; args = [ddu, du, u, t],
        name = "Quintic Hermite Spline (Vector)")

    u = round.(rand(100), digits = 5)
    t = 1.0collect(1:100)
    du = diff(u) ./ diff(t)
    push!(du, 0)
    ddu = diff(du) ./ diff(t)
    push!(ddu, 0)
    test_integral(QuinticHermiteSpline; args = [ddu, du, u, t],
        name = "Quintic Hermite Spline (Vector) with random points")
end

@testset "RegularizationSmooth" begin
    npts = 50
    xmin = 0.0
    xspan = 3 / 2 * π
    x = collect(range(xmin, xmin + xspan, length = npts))
    rng = StableRNG(655)
    x = x + xspan / npts * (rand(rng, npts) .- 0.5)
    # select a subset randomly
    idx = unique(rand(rng, collect(eachindex(x)), 20))
    t = x[unique(idx)]
    npts = length(t)
    ut = sin.(t)
    stdev = 1e-1 * maximum(ut)
    u = ut + stdev * randn(rng, npts)
    # data must be ordered if t̂ is not provided
    idx = sortperm(t)
    tₒ = t[idx]
    uₒ = u[idx]
    test_integral(RegularizationSmooth;
        args = [uₒ, tₒ],
        kwargs = [:alg => :fixed],
        name = "RegularizationSmooth")
end

@testset "Curvefit" begin
    rng = StableRNG(12345)
    model(x, p) = @. p[1] / (1 + exp(x - p[2]))
    t = range(-10, stop = 10, length = 40)
    u = model(t, [1.0, 2.0]) + 0.01 * randn(rng, length(t))
    p0 = [0.5, 0.5]
    A = Curvefit(u, t, model, p0, LBFGS())
    @test_throws DataInterpolations.IntegralNotFoundError integral(A, 0.0, 1.0)
    @test_throws DataInterpolations.IntegralNotFoundError integral(A, 5.0)
end

@testset "BSplineInterpolation" begin
    t = [0, 62.25, 109.66, 162.66, 205.8, 252.3]
    u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
    A = BSplineInterpolation(u, t, 2, :Uniform, :Uniform)
    @test_throws DataInterpolations.IntegralNotFoundError integral(A, 1.0, 100.0)
    @test_throws DataInterpolations.IntegralNotFoundError integral(A, 50.0)
end

@testset "BSplineApprox" begin
    t = [0, 62.25, 109.66, 162.66, 205.8, 252.3]
    u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
    A = BSplineApprox(u, t, 2, 4, :Uniform, :Uniform)
    @test_throws DataInterpolations.IntegralNotFoundError integral(A, 1.0, 100.0)
    @test_throws DataInterpolations.IntegralNotFoundError integral(A, 50.0)
end
