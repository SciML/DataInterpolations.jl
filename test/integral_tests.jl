using DataInterpolations, Test
using QuadGK
using DataInterpolations: integral

function test_integral(func, tspan, name::String)
    t1 = minimum(tspan)
    t2 = maximum(tspan)
    @testset "$name" begin
        qint, err = quadgk(func, t1, t2)
        aint = integral(func, t1, t2)
        @test isapprox(qint, aint, atol = 1e-8)
    end
end

# Linear Interpolation
u = 2.0collect(1:10)
t = 1.0collect(1:10)
A = LinearInterpolation(u, t)

test_integral(A, t, "Linear Interpolation (Vector)")

# Quadratic Interpolation
u = [1.0, 4.0, 9.0, 16.0]
t = [1.0, 2.0, 3.0, 4.0]
A = QuadraticInterpolation(u, t)

test_integral(A, t, "Quadratic Interpolation (Vector)")

# Quadratic forward/backward Interpolation
@testset "QuadraticInterpolation - forward/backward modes" begin
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

# QuadraticSpline Interpolation
u = [0.0, 1.0, 3.0]
t = [-1.0, 0.0, 1.0]

A = QuadraticSpline(u, t)

test_integral(A, t, "Quadratic Spline (Vector)")

# CubicSpline Interpolation
u = [0.0, 1.0, 3.0]
t = [-1.0, 0.0, 1.0]

A = CubicSpline(u, t)

test_integral(A, t, "Cubic Spline Interpolation (Vector)")
