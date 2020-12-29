using DataInterpolations, Test
using Random
using FiniteDifferences
using DataInterpolations: derivative

function test_derivatives(func, tspan, name::String)
  trange = range(minimum(tspan), maximum(tspan), length=32)[2:end-1]
  @testset "$name" begin
    for t in trange
      # Linearly spaced points might lead to evaluations outside
      # trange
      cdiff = central_fdm(5, 1; geom=true)(_t -> func(_t), t)
      adiff = derivative(func, t)
      @test isapprox(cdiff, adiff, atol=1e-8)
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

# Lagrange Interpolation
u = [1.0 4.0 9.0; 1.0 2.0 3.0]
t = [1.0, 2.0, 3.0]
A = LagrangeInterpolation(u,t)

test_derivatives(A, t, "Lagrange Interpolation (Matrix)")

# Akima Interpolation
u = [0.0, 2.0, 1.0, 3.0, 2.0, 6.0, 5.5, 5.5, 2.7, 5.1, 3.0]
t = collect(0.0:10.0)
A = AkimaInterpolation(u, t)

test_derivatives(A, t, "Akima Interpolation")

# QuadraticSpline Interpolation
u = [0.0, 1.0, 3.0]
t = [-1.0, 0.0, 1.0]

A = QuadraticSpline(u,t)

test_derivatives(A, t, "Quadratic Interpolation")

# CubicSpline Interpolation
u = [0.0, 1.0, 3.0]
t = [-1.0, 0.0, 1.0]

A = CubicSpline(u,t)

test_derivatives(A, t, "Cubic Spline Interpolation")
