using DataInterpolations, Test

# Linear Interpolation
u = 2.0collect(1:10)
t = 1.0collect(1:10)
A = LinearInterpolation(u,t)

@test A(1) == u[1]
@test A(5) == u[5]
@test A(5.5) == 11.0

u = vcat(2.0collect(1:10)',3.0collect(1:10)')
A = LinearInterpolation(u,t)

@test A(1) == u[:,1]
@test A(5) == u[:,5]
@test A(5.5) == [11.0,16.5]

# Quadratic Interpolation
u = [1.0, 4.0, 9.0, 16.0]
t = [1.0, 2.0, 3.0, 4.0]
A = QuadraticInterpolation(u,t)

@test A(2.0) == 4.0
@test A(1.5) == 2.25
@test A(3.5) == 12.25
@test A(2.5) == 6.25

u = [1.0 4.0 9.0 16.0; 1.0 4.0 9.0 16.0]
A = QuadraticInterpolation(u,t)

@test A(2.0) == [4.0,4.0]
@test A(1.5) == [2.25,2.25]
@test A(3.5) == [12.25,12.25]
@test A(2.5) == [6.25,6.25]

# Lagrange Interpolation
u = [1.0, 4.0, 9.0, 16.0]
t = [1.0, 2.0, 3.0, 4.0]
A = LagrangeInterpolation(u,t,2)

@test A(2.0) == 4.0
@test A(1.5) == 2.25
@test A(3.5) == 12.25

u = [1.0, 8.0, 27.0]
t = [1.0, 2.0, 3.0]
A = LagrangeInterpolation(u,t,2)

@test A(2.0) == 8.0
@test A(1.5) == 3.0
@test A(2.7) ≈ 20.04

u = [1.0 4.0 9.0 16.0; 1.0 4.0 9.0 16.0]
A = LagrangeInterpolation(u,t,2)

@test A(2.0) == [4.0,4.0]
@test A(1.5) == [2.25,2.25]
@test A(3.5) ≈ [12.25,12.25]

# QuadraticSpline Interpolation
u = [0.0, 1.0, 3.0]
t = [-1.0, 0.0, 1.0]

A = QuadraticSpline(u,t)

#             Solution ->
#             f(x) = (x+1)^2 for x -> [-1.0, 0.0]
#             f(x) = 1+2x    for x -> [0.0, 1.0]

@test A(-0.5) == 0.25
@test A(0.7) == 2.4
@test A(-1.0) == 0.0
@test A(0.0) == 1.0
@test A(1.0) == 3.0


# CubicSpline Interpolation
u = [0.0, 1.0, 3.0]
t = [-1.0, 0.0, 1.0]

A = CubicSpline(u,t)

#             Solution ->
#             f(x) = 1 + 1.5x + x^2 + 0.5x^3 for x -> [-1.0, 0.0]
#             f(x) = 1 + 1.5x + x^2 - 0.5x^3   for x -> [0.0, 1.0]

@test A(-0.5) == 0.4375
@test A(0.5) == 1.9375
@test A(-0.7) == 0.2685
@test A(0.3) == 1.5265
@test A(-1.0) == 0.0
@test A(0.0) == 1.0
@test A(1.0) == 3.0
