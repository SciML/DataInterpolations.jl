using DataInterpolations, Test
using QuadGK
using DataInterpolations: integral
using Optim, ForwardDiff
using RegularizationTools
using StableRNGs

function test_integral(method, u, t; args = [], kwargs = [], name::String)
    func = method(u, t, args...; kwargs..., extrapolate = true)
    t1 = minimum(t)
    t2 = maximum(t)
    @testset "$name" begin
        # integral(A, t1, t2)
        qint, err = quadgk(func, t1, t2; atol = 1e-12, rtol = 1e-12)
        aint = integral(func, t1, t2)
        @test isapprox(qint, aint, atol = 1e-8)

        # integral(A, t)
        qint, err = quadgk(func, t1, (t1 + t2) / 2; atol = 1e-12, rtol = 1e-12)
        aint = integral(func, (t1 + t2) / 2)
        @test isapprox(qint, aint, atol = 1e-8)

        # integrals with extrapolation
        qint, err = quadgk(func, t1 - 5.0, (t1 + t2) / 2; atol = 1e-12, rtol = 1e-12)
        aint = integral(func, t1 - 5.0, (t1 + t2) / 2)
        @test isapprox(qint, aint, atol = 1e-8)

        qint, err = quadgk(func, (t1 + t2) / 2, t2 + 5.0; atol = 1e-12, rtol = 1e-12)
        aint = integral(func, (t1 + t2) / 2, t2 + 5.0)
        @test isapprox(qint, aint, atol = 1e-8)
    end
    func = method(u, t, args...; kwargs...)
    @test_throws DataInterpolations.ExtrapolationError integral(func, t[1] - 1.0)
    @test_throws DataInterpolations.ExtrapolationError integral(func, t[end] + 1.0)
    @test_throws DataInterpolations.ExtrapolationError integral(func, t[1] - 1.0, t[2])
    @test_throws DataInterpolations.ExtrapolationError integral(func, t[1], t[end] + 1.0)
end

@testset "LinearInterpolation" begin
    u = 2.0collect(1:10)
    t = 1.0collect(1:10)
    test_integral(LinearInterpolation, u, t; name = "Linear Interpolation (Vector)")
end

@testset "QuadraticInterpolation" begin
    u = [1.0, 4.0, 9.0, 16.0]
    t = [1.0, 2.0, 3.0, 4.0]
    test_integral(QuadraticInterpolation, u, t; name = "Quadratic Interpolation (Vector)")
    u = [3.0, 0.0, 3.0, 0.0]
    t = [1.0, 2.0, 3.0, 4.0]
    test_integral(QuadraticInterpolation,
        u,
        t;
        args = [:Backward],
        name = "Quadratic Interpolation (Vector)")
end

@testset "LagrangeInterpolation" begin
    u = [1.0, 4.0, 9.0]
    t = [1.0, 2.0, 3.0]
    A = LagrangeInterpolation(u, t)
    @test_throws DataInterpolations.IntegralNotFoundError integral(A, 1.0, 2.0)
    @test_throws DataInterpolations.IntegralNotFoundError integral(A, 5.0)
end

@testset "QuadraticSpline" begin
    u = [0.0, 1.0, 3.0]
    t = [-1.0, 0.0, 1.0]
    test_integral(QuadraticSpline, u, t; name = "Quadratic Spline (Vector)")
end

@testset "CubicSpline" begin
    u = [0.0, 1.0, 3.0]
    t = [-1.0, 0.0, 1.0]
    test_integral(CubicSpline, u, t; name = "Cubic Spline Interpolation (Vector)")
end

@testset "AkimaInterpolation" begin
    u = [0.0, 2.0, 1.0, 3.0, 2.0, 6.0, 5.5, 5.5, 2.7, 5.1, 3.0]
    t = collect(0.0:10.0)
    test_integral(AkimaInterpolation, u, t; name = "Akima Interpolation (Vector)")
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
    test_integral(RegularizationSmooth,
        uₒ,
        tₒ;
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
