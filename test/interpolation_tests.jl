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
