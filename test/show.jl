using DataInterpolations
using CurveFit, StableRNGs
using RegularizationTools

t = [1.0, 2.0, 3.0, 4.0, 5.0]
x = [1.0, 2.0, 3.0, 4.0, 5.0]

@testset "Generic Cases" begin
    function test_show_line(A)
        @testset "$(nameof(typeof(A)))" begin
            @test startswith(
                sprint(io -> show(io, MIME"text/plain"(), A)),
                "$(nameof(typeof(A))) with $(length(A.t)) points\n"
            )
        end
    end
    methods = [
        LinearInterpolation(x, t),
        AkimaInterpolation(x, t),
        QuadraticSpline(x, t),
        CubicSpline(x, t),
    ]
    test_show_line.(methods)
end

@testset "Specific Cases" begin
    @testset "QuadraticInterpolation" begin
        A = QuadraticInterpolation(x, t)
        @test startswith(
            sprint(io -> show(io, MIME"text/plain"(), A)),
            "QuadraticInterpolation with 5 points, Forward mode\n"
        )
    end
    @testset "LagrangeInterpolation" begin
        A = LagrangeInterpolation(x, t)
        @test startswith(
            sprint(io -> show(io, MIME"text/plain"(), A)),
            "LagrangeInterpolation with 5 points, with order 4\n"
        )
    end
    @testset "ConstantInterpolation" begin
        A = ConstantInterpolation(x, t)
        @test startswith(
            sprint(io -> show(io, MIME"text/plain"(), A)),
            "ConstantInterpolation with 5 points, in left direction\n"
        )
    end
    @testset "BSplineInterpolation" begin
        A = BSplineInterpolation(x, t, 3, :Uniform, :Uniform)
        @test startswith(
            sprint(io -> show(io, MIME"text/plain"(), A)),
            "BSplineInterpolation with 5 points, with degree 3\n"
        )
    end
    @testset "BSplineApprox" begin
        A = BSplineApprox(x, t, 2, 4, :Uniform, :Uniform)
        @test startswith(
            sprint(io -> show(io, MIME"text/plain"(), A)),
            "BSplineApprox with 5 points, with degree 2, number of control points 4\n"
        )
    end
end

@testset "CurveFit" begin
    rng = StableRNG(12345)
    model(x, p) = @. p[1] / (1 + exp(x - p[2]))
    t = range(-10, stop = 10, length = 40)
    u = model(t, [1.0, 2.0]) + 0.01 * randn(rng, length(t))
    p0 = [0.5, 0.5]
    A = Curvefit(u, t, model, p0)
    @test startswith(
        sprint(io -> show(io, MIME"text/plain"(), A)),
        "Curvefit with 40 points.\n"
    )
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
    stdev = 1.0e-1 * maximum(ut)
    u = ut + stdev * randn(rng, npts)
    # data must be ordered if t̂ is not provided
    idx = sortperm(t)
    tₒ = t[idx]
    uₒ = u[idx]
    A = RegularizationSmooth(uₒ, tₒ; alg = :fixed)
    @test startswith(
        sprint(io -> show(io, MIME"text/plain"(), A)),
        "RegularizationSmooth with 15 points, with regularization coefficient 1.0\n"
    )
end
