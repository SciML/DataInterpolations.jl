using DataInterpolations, Test
using FiniteDifferences
using DataInterpolations: derivative

function test_derivatives(func, tspan, name::String)
    trange = range(minimum(tspan) - 5.0, maximum(tspan) + 5.0, length = 32)
    @testset "$name" begin
        for t in trange
            cdiff = central_fdm(5, 1; geom = true)(_t -> func(_t), t)
            adiff = derivative(func, t)
            @test isapprox(cdiff, adiff, atol = 1e-8)
        end
    end
end

# Linear Interpolation
u = 2.0collect(1:10)
t = 1.0collect(1:10)
A = LinearInterpolation(u, t)

test_derivatives(A, t, "Linear Interpolation (Vector)")

u = vcat(2.0collect(1:10)', 3.0collect(1:10)')
A = LinearInterpolation(u, t)

test_derivatives(A, t, "Linear Interpolation (Matrix)")

# Quadratic Interpolation
u = [1.0, 4.0, 9.0, 16.0]
t = [1.0, 2.0, 3.0, 4.0]
A = QuadraticInterpolation(u, t)

test_derivatives(A, t, "Quadratic Interpolation (Vector)")

Ab = QuadraticInterpolation(u, t, :Backward)

test_derivatives(Ab, t, "Quadratic Interpolation (Vector), backward")

u = [1.0 4.0 9.0 16.0; 1.0 4.0 9.0 16.0]
A = QuadraticInterpolation(u, t)

test_derivatives(A, t, "Quadratic Interpolation (Matrix)")

@testset "Backward Quadratic Interpolation" begin
    u = [0.5, 0.0, 0.5, 0.0]
    t = [1.0, 2.0, 3.0, 4.0]
    A_f = QuadraticInterpolation(u, t)
    A_b = QuadraticInterpolation(u, t, :Backward)
    @test derivative(A_f, 1.5) ≈ -0.5
    @test derivative(A_b, 1.5) ≈ -0.5
    @test derivative(A_f, 2.25) ≈ 0.75
    @test derivative(A_b, 2.25) ≈ 0.25
    @test derivative(A_f, 2.75) ≈ 0.25
    @test derivative(A_b, 2.75) ≈ 0.75
    @test derivative(A_f, 3.5) ≈ -0.5
    @test derivative(A_b, 3.5) ≈ -0.5
end

# Lagrange Interpolation
u = [1.0, 4.0, 9.0]
t = [1.0, 2.0, 3.0]
A = LagrangeInterpolation(u, t)

test_derivatives(A, t, "Lagrange Interpolation (Vector)")

# Lagrange Interpolation
u = [1.0 4.0 9.0; 1.0 2.0 3.0]
t = [1.0, 2.0, 3.0]
A = LagrangeInterpolation(u, t)

test_derivatives(A, t, "Lagrange Interpolation (Matrix)")

# Lagrange Interpolation
u = [[1.0, 4.0, 9.0], [3.0, 7.0, 4.0], [5.0, 4.0, 1.0]]
t = [1.0, 2.0, 3.0]
A = LagrangeInterpolation(u, t)

test_derivatives(A, t, "Lagrange Interpolation (Vector of Vectors)")

# Lagrange Interpolation
u = [[3.0 1.0 4.0; 1.0 5.0 9.0], [2.0 6.0 5.0; 3.0 5.0 8.0], [9.0 7.0 9.0; 3.0 2.0 3.0]]
t = [1.0, 2.0, 3.0]
A = LagrangeInterpolation(u, t)

test_derivatives(A, t, "Lagrange Interpolation (Vector of Matrices)")

# Akima Interpolation
u = [0.0, 2.0, 1.0, 3.0, 2.0, 6.0, 5.5, 5.5, 2.7, 5.1, 3.0]
t = collect(0.0:10.0)
A = AkimaInterpolation(u, t)

test_derivatives(A, t, "Akima Interpolation")

@testset "Akima smooth derivative at end points" begin
    @test derivative(A, t[1]) ≈ derivative(A, nextfloat(t[1]))
    @test derivative(A, t[end]) ≈ derivative(A, prevfloat(t[end]))
end

# QuadraticSpline Interpolation
u = [0.0, 1.0, 3.0]
t = [-1.0, 0.0, 1.0]

A = QuadraticSpline(u, t)

test_derivatives(A, t, "Quadratic Interpolation (Vector)")

# QuadraticSpline Interpolation
u = [[1.0, 2.0, 9.0], [3.0, 7.0, 5.0], [5.0, 4.0, 1.0]]
t = [-1.0, 0.0, 1.0]

A = QuadraticSpline(u, t)

test_derivatives(A, t, "Quadratic Interpolation (Vector of Vectors)")

# QuadraticSpline Interpolation
u = [[1.0 4.0 9.0; 5.0 9.0 2.0], [3.0 7.0 4.0; 6.0 5.0 3.0], [5.0 4.0 1.0; 2.0 3.0 8.0]]
t = [-1.0, 0.0, 1.0]

A = QuadraticSpline(u, t)

test_derivatives(A, t, "Quadratic Interpolation (Vector of Matrices)")

# CubicSpline Interpolation
u = [0.0, 1.0, 3.0]
t = [-1.0, 0.0, 1.0]

A = CubicSpline(u, t)

test_derivatives(A, t, "Cubic Spline Interpolation (Vector)")

# CubicSpline Interpolation
u = [[1.0, 2.0, 9.0], [3.0, 7.0, 5.0], [5.0, 4.0, 1.0]]
t = [-1.0, 0.0, 1.0]

A = CubicSpline(u, t)

test_derivatives(A, t, "Cubic Spline Interpolation (Vector of Vectors)")

# CubicSpline Interpolation
u = [[1.0 4.0 9.0; 5.0 9.0 2.0], [3.0 7.0 4.0; 6.0 5.0 3.0], [5.0 4.0 1.0; 2.0 3.0 8.0]]
t = [-1.0, 0.0, 1.0]

A = CubicSpline(u, t)

test_derivatives(A, t, "Cubic Spline Interpolation (Vector of Matrices)")

# BSpline Interpolation and Approximation
t = [0, 62.25, 109.66, 162.66, 205.8, 252.3]
u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]

A = BSplineInterpolation(u, t, 2, :Uniform, :Uniform)

test_derivatives(A, t, "BSpline Interpolation (Uniform, Uniform)")

A = BSplineInterpolation(u, t, 2, :ArcLen, :Average)

test_derivatives(A, t, "BSpline Interpolation (Arclen, Average)")

A = BSplineApprox(u, t, 3, 4, :Uniform, :Uniform)

test_derivatives(A, t, "BSpline Approx (Uniform, Uniform)")

using Symbolics

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
