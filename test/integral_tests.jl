using DataInterpolations, Test
using QuadGK
using DataInterpolations: integral
using Optim, ForwardDiff
using RegularizationTools

function test_integral(func, tspan, name::String)
    t1 = minimum(tspan)
    t2 = maximum(tspan)
    @testset "$name" begin
        qint, err = quadgk(func, t1, t2)
        aint = integral(func, t1, t2)
        @test isapprox(qint, aint, atol = 1e-8)
    end
end

@testset "LinearInterpolation" begin
    u = 2.0collect(1:10)
    t = 1.0collect(1:10)
    A = LinearInterpolation(u, t)
    test_integral(A, t, "Linear Interpolation (Vector)")
end

@testset "QuadraticInterpolation" begin
    u = [1.0, 4.0, 9.0, 16.0]
    t = [1.0, 2.0, 3.0, 4.0]
    A = QuadraticInterpolation(u, t)
    test_integral(A, t, "Quadratic Interpolation (Vector)")
    u = [3.0, 0.0, 3.0, 0.0]
    t = [1.0, 2.0, 3.0, 4.0]
    A_f = QuadraticInterpolation(u, t)
    A_b = QuadraticInterpolation(u, t, :Backward)
    @test integral(A_f, 1.0, 2.0) ≈ 1.0
    @test integral(A_f, 2.0, 3.0) ≈ 2.0
    @test integral(A_f, 3.0, 4.0) ≈ 2.0
    @test integral(A_f, 1.0, 3.0) ≈ 3.0
    @test integral(A_f, 2.0, 4.0) ≈ 4.0
    @test integral(A_f, 1.0, 4.0) ≈ 5.0
    @test integral(A_b, 1.0, 2.0) ≈ 1.0
    @test integral(A_b, 2.0, 3.0) ≈ 1.0
    @test integral(A_b, 3.0, 4.0) ≈ 2.0
    @test integral(A_b, 1.0, 3.0) ≈ 2.0
    @test integral(A_b, 2.0, 4.0) ≈ 3.0
    @test integral(A_b, 1.0, 4.0) ≈ 4.0
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
    A = QuadraticSpline(u, t)
    test_integral(A, t, "Quadratic Spline (Vector)")
end

@testset "CubicSpline" begin
    u = [0.0, 1.0, 3.0]
    t = [-1.0, 0.0, 1.0]
    A = CubicSpline(u, t)
    test_integral(A, t, "Cubic Spline Interpolation (Vector)")
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
    A = RegularizationSmooth(uₒ, tₒ; alg = :fixed)
    test_integral(A, tₒ, "RegularizationSmooth")
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
