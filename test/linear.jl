using DataInterpolations, Test

u = 2.0collect(1:10)
t = 1.0collect(1:10)
A = LinearInterpolation(u,t)

@test A(5) == u[5]
@test A(5.5) == 11.0

u = vcat(2.0collect(1:10)',3.0collect(1:10)')
A = LinearInterpolation(u,t)
@test A(5) == u[:,5]
@test A(5.5) == [11.0,16.5]
