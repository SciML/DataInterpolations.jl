using DataInterpolations
u = 2.0collect(1:10)
t = 1.0collect(1:10)
A = LinearInterpolation(u, t)

for i in 1:10
    @test u[i] == A[i]
end

for i in 11:20
    @test t[i - 10] == A[i]
end

A = LinearInterpolation{false}(u, t, true)
for i in 1:10
    @test u[i] == A[i]
end

using Symbolics
u = 2.0collect(1:10)
t = 1.0collect(1:10)
A = LinearInterpolation(u, t)

@variables t x(t)
substitute(A(t), Dict(t => x))
