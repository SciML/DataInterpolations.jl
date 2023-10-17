using DataInterpolations, Test
using FiniteDifferences
using DataInterpolations: derivative
# using Symbolics

function test_derivatives(method, u, t, args...; name::String)
    func = method(u, t, args...; extrapolate = true)
    trange = range(minimum(t) - 5.0, maximum(t) + 5.0, length = 32)
    @testset "$name" begin
        for t in trange
            cdiff = central_fdm(5, 1; geom = true)(_t -> func(_t), t)
            adiff = derivative(func, t)
            @test isapprox(cdiff, adiff, atol = 1e-8)
        end
    end
    func = method(u, t, args...)
    @test_throws DataInterpolations.ExtrapolationError derivative(func, t[1] - 1.0)
    @test_throws DataInterpolations.ExtrapolationError derivative(func, t[end] + 1.0)
end

@testset "Linear Interpolation" begin
    u = 2.0collect(1:10)
    t = 1.0collect(1:10)
    test_derivatives(LinearInterpolation, u, t; name = "Linear Interpolation (Vector)")
    u = vcat(2.0collect(1:10)', 3.0collect(1:10)')
    test_derivatives(LinearInterpolation, u, t; name = "Linear Interpolation (Matrix)")
end

@testset "Quadratic Interpolation" begin
    u = [1.0, 4.0, 9.0, 16.0]
    t = [1.0, 2.0, 3.0, 4.0]
    test_derivatives(QuadraticInterpolation,
        u,
        t;
        name = "Quadratic Interpolation (Vector)")
    test_derivatives(QuadraticInterpolation,
        u,
        t,
        :Backward;
        name = "Quadratic Interpolation (Vector), backward")
    u = [1.0 4.0 9.0 16.0; 1.0 4.0 9.0 16.0]
    test_derivatives(QuadraticInterpolation,
        u,
        t;
        name = "Quadratic Interpolation (Matrix)")
end

@testset "Lagrange Interpolation" begin
    u = [1.0, 4.0, 9.0]
    t = [1.0, 2.0, 3.0]
    test_derivatives(LagrangeInterpolation, u, t; name = "Lagrange Interpolation (Vector)")
    u = [1.0 4.0 9.0; 1.0 2.0 3.0]
    test_derivatives(LagrangeInterpolation, u, t; name = "Lagrange Interpolation (Matrix)")
    u = [[1.0, 4.0, 9.0], [3.0, 7.0, 4.0], [5.0, 4.0, 1.0]]
    test_derivatives(LagrangeInterpolation,
        u,
        t;
        name = "Lagrange Interpolation (Vector of Vectors)")
    u = [[3.0 1.0 4.0; 1.0 5.0 9.0], [2.0 6.0 5.0; 3.0 5.0 8.0], [9.0 7.0 9.0; 3.0 2.0 3.0]]
    test_derivatives(LagrangeInterpolation,
        u,
        t;
        name = "Lagrange Interpolation (Vector of Matrices)")
end

@testset "Akima Interpolation" begin
    u = [0.0, 2.0, 1.0, 3.0, 2.0, 6.0, 5.5, 5.5, 2.7, 5.1, 3.0]
    t = collect(0.0:10.0)
    test_derivatives(AkimaInterpolation, u, t; name = "Akima Interpolation")
    @testset "Akima smooth derivative at end points" begin
        A = AkimaInterpolation(u, t)
        @test derivative(A, t[1]) ≈ derivative(A, nextfloat(t[1]))
        @test derivative(A, t[end]) ≈ derivative(A, prevfloat(t[end]))
    end
end

@testset "Quadratic Spline" begin
    u = [0.0, 1.0, 3.0]
    t = [-1.0, 0.0, 1.0]
    test_derivatives(QuadraticSpline, u, t; name = "Quadratic Interpolation (Vector)")
    u = [[1.0, 2.0, 9.0], [3.0, 7.0, 5.0], [5.0, 4.0, 1.0]]
    test_derivatives(QuadraticSpline,
        u,
        t;
        name = "Quadratic Interpolation (Vector of Vectors)")
    u = [[1.0 4.0 9.0; 5.0 9.0 2.0], [3.0 7.0 4.0; 6.0 5.0 3.0], [5.0 4.0 1.0; 2.0 3.0 8.0]]
    test_derivatives(QuadraticSpline,
        u,
        t;
        name = "Quadratic Interpolation (Vector of Matrices)")
end

@testset "Cubic Spline" begin
    u = [0.0, 1.0, 3.0]
    t = [-1.0, 0.0, 1.0]
    test_derivatives(CubicSpline, u, t; name = "Cubic Spline Interpolation (Vector)")
    u = [[1.0, 2.0, 9.0], [3.0, 7.0, 5.0], [5.0, 4.0, 1.0]]
    test_derivatives(CubicSpline,
        u,
        t;
        name = "Cubic Spline Interpolation (Vector of Vectors)")
    u = [[1.0 4.0 9.0; 5.0 9.0 2.0], [3.0 7.0 4.0; 6.0 5.0 3.0], [5.0 4.0 1.0; 2.0 3.0 8.0]]
    test_derivatives(CubicSpline,
        u,
        t;
        name = "Cubic Spline Interpolation (Vector of Matrices)")
end

@testset "BSplines" begin
    t = [0, 62.25, 109.66, 162.66, 205.8, 252.3]
    u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
    test_derivatives(BSplineInterpolation,
        u,
        t,
        2,
        :Uniform,
        :Uniform;
        name = "BSpline Interpolation (Uniform, Uniform)")
    test_derivatives(BSplineInterpolation,
        u,
        t,
        2,
        :ArcLen,
        :Average;
        name = "BSpline Interpolation (Arclen, Average)")
    test_derivatives(BSplineApprox,
        u,
        t,
        3,
        4,
        :Uniform,
        :Uniform;
        name = "BSpline Approx (Uniform, Uniform)")
end

@testset "Symbolic derivatives" begin
    u = [0.0, 1.5, 0.0]
    t = [0.0, 0.5, 1.0]
    A = QuadraticSpline(u, t)
    @variables τ, ω(τ)
    D = Symbolics.Differential(τ)
    expr = A(ω)
    @test isequal(Symbolics.derivative(expr, τ), D(ω) * DataInterpolations.derivative(A, ω))

    derivexpr = expand_derivatives(substitute(D(A(ω)), Dict(ω => 0.5τ)))
    symfunc = Symbolics.build_function(derivexpr, τ; expression = Val{false})
    @test symfunc(0.5) == 0.5 * 3

    u = [0.0, 1.5, 0.0]
    t = [0.0, 0.5, 1.0]
    @variables τ
    D = Symbolics.Differential(τ)
    f = LinearInterpolation(u, t)
    df = expand_derivatives(D(f(τ)))
    symfunc = Symbolics.build_function(df, τ; expression = Val{false})
    ts = 0.0:0.1:1.0
    @test all(map(ti -> symfunc(ti) == derivative(f, ti), ts))
end
