using DataInterpolations, Test
using FindFirstFunctions: searchsortedfirstcorrelated
using FiniteDifferences
using DataInterpolations: derivative
using Symbolics
using StableRNGs
using RegularizationTools
using Optim
using ForwardDiff

function test_derivatives(method; args = [], kwargs = [], name::String)
    func = method(args...; kwargs..., extrapolate = true)
    (; t) = func
    trange = collect(range(minimum(t) - 5.0, maximum(t) + 5.0, step = 0.1))
    trange_exclude = filter(x -> !in(x, t), trange)
    @testset "$name" begin
        # Rest of the points
        for _t in trange_exclude
            cdiff = central_fdm(5, 1; geom = true)(func, _t)
            adiff = derivative(func, _t)
            @test isapprox(cdiff, adiff, atol = 1e-8)
            adiff2 = derivative(func, _t, 2)
            cdiff2 = central_fdm(5, 1; geom = true)(t -> derivative(func, t), _t)
            @test isapprox(cdiff2, adiff2, atol = 1e-8)
        end

        # Interpolation time points
        for _t in t[2:(end - 1)]
            if func isa BSplineInterpolation || func isa BSplineApprox ||
               func isa CubicHermiteSpline
                fdiff = forward_fdm(5, 1; geom = true)(func, _t)
                fdiff2 = forward_fdm(5, 1; geom = true)(t -> derivative(func, t), _t)
            else
                fdiff = backward_fdm(5, 1; geom = true)(func, _t)
                fdiff2 = backward_fdm(5, 1; geom = true)(t -> derivative(func, t), _t)
            end
            adiff = derivative(func, _t)
            adiff2 = derivative(func, _t, 2)
            @test isapprox(fdiff, adiff, atol = 1e-8)
            @test isapprox(fdiff2, adiff2, atol = 1e-8)
            # Cached index
            if hasproperty(func, :idx_prev)
                @test abs(func.idx_prev[] -
                          searchsortedfirstcorrelated(func.t, _t, func.idx_prev[])) <= 1
            end
        end

        # t = t0
        fdiff = forward_fdm(5, 1; geom = true)(func, t[1])
        adiff = derivative(func, t[1])
        @test isapprox(fdiff, adiff, atol = 1e-8)
        if !(func isa BSplineInterpolation || func isa BSplineApprox)
            fdiff2 = forward_fdm(5, 1; geom = true)(t -> derivative(func, t), t[1])
            adiff2 = derivative(func, t[1], 2)
            @test isapprox(fdiff2, adiff2, atol = 1e-8)
        end

        # t = tend
        fdiff = backward_fdm(5, 1; geom = true)(func, t[end])
        adiff = derivative(func, t[end])
        @test isapprox(fdiff, adiff, atol = 1e-8)
        if !(func isa BSplineInterpolation || func isa BSplineApprox)
            fdiff2 = backward_fdm(5, 1; geom = true)(t -> derivative(func, t), t[end])
            adiff2 = derivative(func, t[end], 2)
            @test isapprox(fdiff2, adiff2, atol = 1e-8)
        end
    end
    @test_throws DataInterpolations.DerivativeNotFoundError derivative(
        func, t[1], 3)
    func = method(args...)
    @test_throws DataInterpolations.ExtrapolationError derivative(func, t[1] - 1.0)
    @test_throws DataInterpolations.ExtrapolationError derivative(func, t[end] + 1.0)
    @test_throws DataInterpolations.DerivativeNotFoundError derivative(
        func, t[1], 3)
end

@testset "Linear Interpolation" begin
    u = vcat(collect(1:5), 2 * collect(6:10))
    t = 1.0collect(1:10)
    test_derivatives(
        LinearInterpolation; args = [u, t], name = "Linear Interpolation (Vector)")
    u = vcat(2.0collect(1:10)', 3.0collect(1:10)')
    test_derivatives(
        LinearInterpolation; args = [u, t], name = "Linear Interpolation (Matrix)")
end

@testset "Quadratic Interpolation" begin
    u = [1.0, 4.0, 9.0, 16.0]
    t = [1.0, 2.0, 3.0, 4.0]
    test_derivatives(QuadraticInterpolation, args = [u, t],
        name = "Quadratic Interpolation (Vector)")
    test_derivatives(QuadraticInterpolation;
        args = [u, t, :Backward],
        name = "Quadratic Interpolation (Vector), backward")
    u = [1.0 4.0 9.0 16.0; 1.0 4.0 9.0 16.0]
    test_derivatives(QuadraticInterpolation;
        args = [u, t],
        name = "Quadratic Interpolation (Matrix)")
end

@testset "Lagrange Interpolation" begin
    u = [1.0, 4.0, 9.0]
    t = [1.0, 2.0, 3.0]
    test_derivatives(
        LagrangeInterpolation; args = [u, t], name = "Lagrange Interpolation (Vector)")
    u = [1.0 4.0 9.0; 1.0 2.0 3.0]
    test_derivatives(
        LagrangeInterpolation; args = [u, t], name = "Lagrange Interpolation (Matrix)")
    u = [[1.0, 4.0, 9.0], [3.0, 7.0, 4.0], [5.0, 4.0, 1.0]]
    test_derivatives(LagrangeInterpolation; args = [u, t],
        name = "Lagrange Interpolation (Vector of Vectors)")
    u = [[3.0 1.0 4.0; 1.0 5.0 9.0], [2.0 6.0 5.0; 3.0 5.0 8.0], [9.0 7.0 9.0; 3.0 2.0 3.0]]
    test_derivatives(LagrangeInterpolation; args = [u, t],
        name = "Lagrange Interpolation (Vector of Matrices)")
end

@testset "Akima Interpolation" begin
    u = [0.0, 2.0, 1.0, 3.0, 2.0, 6.0, 5.5, 5.5, 2.7, 5.1, 3.0]
    t = collect(0.0:10.0)
    test_derivatives(AkimaInterpolation; args = [u, t], name = "Akima Interpolation")
    @testset "Akima smooth derivative at end points" begin
        A = AkimaInterpolation(u, t)
        @test derivative(A, t[1]) ≈ derivative(A, nextfloat(t[1]))
        @test derivative(A, t[end]) ≈ derivative(A, prevfloat(t[end]))
    end
end

@testset "Quadratic Spline" begin
    u = [0.0, 1.0, 3.0]
    t = [-1.0, 0.0, 1.0]
    test_derivatives(
        QuadraticSpline; args = [u, t], name = "Quadratic Interpolation (Vector)")
    u = [[1.0, 2.0, 9.0], [3.0, 7.0, 5.0], [5.0, 4.0, 1.0]]
    test_derivatives(QuadraticSpline; args = [u, t],
        name = "Quadratic Interpolation (Vector of Vectors)")
    u = [[1.0 4.0 9.0; 5.0 9.0 2.0], [3.0 7.0 4.0; 6.0 5.0 3.0], [5.0 4.0 1.0; 2.0 3.0 8.0]]
    test_derivatives(QuadraticSpline; args = [u, t],
        name = "Quadratic Interpolation (Vector of Matrices)")
end

@testset "Cubic Spline" begin
    u = [0.0, 1.0, 3.0]
    t = [-1.0, 0.0, 1.0]
    test_derivatives(
        CubicSpline; args = [u, t], name = "Cubic Spline Interpolation (Vector)")
    u = [[1.0, 2.0, 9.0], [3.0, 7.0, 5.0], [5.0, 4.0, 1.0]]
    test_derivatives(CubicSpline; args = [u, t],
        name = "Cubic Spline Interpolation (Vector of Vectors)")
    u = [[1.0 4.0 9.0; 5.0 9.0 2.0], [3.0 7.0 4.0; 6.0 5.0 3.0], [5.0 4.0 1.0; 2.0 3.0 8.0]]
    test_derivatives(CubicSpline; args = [u, t],
        name = "Cubic Spline Interpolation (Vector of Matrices)")
end

@testset "BSplines" begin
    t = [0, 62.25, 109.66, 162.66, 205.8, 252.3]
    u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
    test_derivatives(BSplineInterpolation;
        args = [u, t, 2,
            :Uniform,
            :Uniform],
        name = "BSpline Interpolation (Uniform, Uniform)")
    test_derivatives(BSplineInterpolation;
        args = [u, t, 2,
            :ArcLen,
            :Average],
        name = "BSpline Interpolation (Arclen, Average)")
    test_derivatives(BSplineApprox;
        args = [u, t,
            3,
            4,
            :Uniform,
            :Uniform],
        name = "BSpline Approx (Uniform, Uniform)")
end

@testset "Cubic Hermite Spline" begin
    du = [-0.047, -0.058, 0.054, 0.012, -0.068, 0.0]
    u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
    t = [0.0, 62.25, 109.66, 162.66, 205.8, 252.3]
    test_derivatives(CubicHermiteSpline; args = [du, u, t],
        name = "Cubic Hermite Spline")
    A = CubicHermiteSpline(du, u, t; extrapolate = true)
    @test derivative.(Ref(A), t) ≈ du
    @test derivative(A, 100.0)≈0.0105409 rtol=1e-5
    @test derivative(A, 300.0)≈-0.0806717 rtol=1e-5
end

@testset "Quintic Hermite Spline" begin
    ddu = [0.0, -0.00033, 0.0051, -0.0067, 0.0029, 0.0]
    du = [-0.047, -0.058, 0.054, 0.012, -0.068, 0.0]
    u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
    t = [0.0, 62.25, 109.66, 162.66, 205.8, 252.3]
    test_derivatives(QuinticHermiteSpline; args = [ddu, du, u, t],
        name = "Quintic Hermite Spline")
    A = QuinticHermiteSpline(ddu, du, u, t; extrapolate = true)
    @test derivative.(Ref(A), t) ≈ du
    @test derivative.(Ref(A), t, 2) ≈ ddu
    @test derivative(A, 100.0)≈0.0103916 rtol=1e-5
    @test derivative(A, 300.0)≈0.0331361 rtol=1e-5
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
    test_derivatives(RegularizationSmooth; args = [uₒ, tₒ],
        kwargs = [:alg => :fixed],
        name = "RegularizationSmooth")
end

@testset "Curvefit" begin
    rng = StableRNG(12345)
    model(x, p) = @. p[1] / (1 + exp(x - p[2]))
    t = range(-10, stop = 10, length = 40)
    u = model(t, [1.0, 2.0]) + 0.01 * randn(rng, length(t))
    p0 = [0.5, 0.5]
    test_derivatives(Curvefit; args = [u, t, model, p0, LBFGS()], name = "Curvefit")
end

@testset "Symbolic derivatives" begin
    u = [0.0, 1.5, 0.0]
    t = [0.0, 0.5, 1.0]
    A = QuadraticSpline(u, t)
    @variables τ, ω(τ)
    D = Symbolics.Differential(τ)
    D2 = Symbolics.Differential(τ)^2
    expr = A(ω)
    @test isequal(Symbolics.derivative(expr, τ), D(ω) * DataInterpolations.derivative(A, ω))

    derivexpr1 = expand_derivatives(substitute(D(A(ω)), Dict(ω => 0.5τ)))
    derivexpr2 = expand_derivatives(substitute(D2(A(ω)), Dict(ω => 0.5τ)))
    symfunc1 = Symbolics.build_function(derivexpr1, τ; expression = Val{false})
    symfunc2 = Symbolics.build_function(derivexpr2, τ; expression = Val{false})
    @test symfunc1(0.5) == 0.5 * 3
    @test symfunc2(0.5) == 0.5 * 6

    u = [0.0, 1.5, 0.0]
    t = [0.0, 0.5, 1.0]
    @variables τ
    D = Symbolics.Differential(τ)
    D2 = Symbolics.Differential(τ)^2
    D3 = Symbolics.Differential(τ)^3
    f = LinearInterpolation(u, t)
    df = expand_derivatives(D(f(τ)))
    df2 = expand_derivatives(D2(f(τ)))
    df3 = expand_derivatives(D3(f(τ)))
    symfunc1 = Symbolics.build_function(df, τ; expression = Val{false})
    symfunc2 = Symbolics.build_function(df2, τ; expression = Val{false})
    symfunc3 = Symbolics.build_function(df3, τ; expression = Val{false})
    ts = 0.0:0.1:1.0
    @test all(map(ti -> symfunc1(ti) == derivative(f, ti), ts))
    @test all(map(ti -> symfunc2(ti) == derivative(f, ti, 2), ts))
    @test_throws DataInterpolations.DerivativeNotFoundError symfunc3(ts[1])
end

@testset "Jacobian tests" begin
    u = rand(5)
    t = 0:4
    interp = LinearInterpolation(u, t, extrapolate = true)
    grad1 = ForwardDiff.derivative(interp, 2.4)

    myvec = rand(20) .* 4.0
    interp(myvec)

    grad = ForwardDiff.jacobian(interp, myvec)
end
