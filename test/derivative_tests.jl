using DataInterpolations, Test
using Random
using FiniteDifferences
using DataInterpolations: derivative

function test_derivatives(func, tspan, name::String)
    trange = range(minimum(tspan), maximum(tspan), length=32)[2:end - 1]
    @testset "$name" begin
        for t in trange
            cdiff = central_fdm(5, 1)(func, t)
            adiff = derivative(func, t)
            @test isapprox(cdiff, adiff)
        end
    end
end

# Linear Interpolation
u = 2.0collect(1:10)
t = 1.0collect(1:10)
A = LinearInterpolation(u,t)

test_derivatives(A, t, "Linear Interpolation (Vector)")

u = vcat(2.0collect(1:10)', 3.0collect(1:10)')
A = LinearInterpolation(u,t)

test_derivatives(A, t, "Linear Interpolation (Matrix)")

# Quadratic Interpolation
u = [1.0, 4.0, 9.0, 16.0]
t = [1.0, 2.0, 3.0, 4.0]
A = QuadraticInterpolation(u,t)

test_derivatives(A, t, "Quadratic Interpolation (Vector)")

u = [1.0 4.0 9.0 16.0; 1.0 4.0 9.0 16.0]
A = QuadraticInterpolation(u,t)

test_derivatives(A, t, "Quadratic Interpolation (Matrix)")

# Lagrange Interpolation
u = [1.0, 4.0, 9.0]
t = [1.0, 2.0, 3.0]
A = LagrangeInterpolation(u,t)

test_derivatives(A, t, "Lagrange Interpolation (Vector)")

# Akima Interpolation
u = [0.0, 2.0, 1.0, 3.0, 2.0, 6.0, 5.5, 5.5, 2.7, 5.1, 3.0]
t = collect(0.0:10.0)
A = AkimaInterpolation(u, t)

test_derivatives(A, t, "Akima Interpolation")

# QuadraticSpline Interpolation
u = [0.0, 1.0, 3.0]
t = [-1.0, 0.0, 1.0]

A = QuadraticSpline(u,t)

#             Solution ->
#             f(x) = (x+1)^2 for x -> [-1.0, 0.0]
#             f(x) = 1+2x    for x -> [0.0, 1.0]

test_derivatives(A, t, "Quadratic Interpolation")


# CubicSpline Interpolation
u = [0.0, 1.0, 3.0]
t = [-1.0, 0.0, 1.0]

A = CubicSpline(u,t)

#             Solution ->
#             f(x) = 1 + 1.5x + x^2 + 0.5x^3 for x -> [-1.0, 0.0]
#             f(x) = 1 + 1.5x + x^2 - 0.5x^3   for x -> [0.0, 1.0]

test_derivatives(A, t, "Cubic Spline Interpolation")

# BSpline Interpolation and Approximation
t = [0,62.25,109.66,162.66,205.8,252.3]
u = [14.7,11.51,10.41,14.95,12.24,11.22]

A = BSplineInterpolation(u,t,2,:Uniform,:Uniform)

# test_derivatives(A, t, "BSpline Interpolation (Uniform, Uniform)")

A = BSplineInterpolation(u,t,2,:ArcLen,:Average)

# test_derivatives(A, t, "BSpline Interpolation (Arclen, Average)")

A = BSplineApprox(u,t,2,4,:Uniform,:Uniform)

# test_derivatives(A, t, "BSpline Approx (Uniform, Uniform)")