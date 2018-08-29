using DataInterpolations, Test

# Linear Interpolation
u = 2.0collect(1:10)
t = 1.0collect(1:10)
A = LinearInterpolation(u,t)

@test A(5) == u[5]
@test A(5.5) == 11.0

u = vcat(2.0collect(1:10)',3.0collect(1:10)')
A = LinearInterpolation(u,t)
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

@test A(1.5) == 2.25
@test A(3.5) == 12.25

u = [1.0, 8.0, 27.0]
t = [1.0, 2.0, 3.0]
A = LagrangeInterpolation(u,t,2)

@test A(1.5) == 3.0
@test A(2.7) ≈ 20.04

u = [1.0 4.0 9.0 16.0; 1.0 4.0 9.0 16.0]
A = LagrangeInterpolation(u,t,2)

@test A(1.5) == [2.25,2.25]
@test A(3.5) ≈ [12.25,12.25]
