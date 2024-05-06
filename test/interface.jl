using DataInterpolations
u = 2.0collect(1:10)
t = 1.0collect(1:10)
A = LinearInterpolation(u, t)

for i in 1:10
    @test u[i] == A.u[i]
end

for i in 1:10
    @test t[i] == A.t[i]
end

using Symbolics
u = 2.0collect(1:10)
t = 1.0collect(1:10)
A = LinearInterpolation(u, t)

@variables t x(t)
substitute(A(t), Dict(t => x))
