u = 2.0collect(1:10)
t = 1.0collect(1:10)
A = LinearInterpolation(u,t)

@test length(A) == 20
for i in 1:10
  @test u[i] == A[i]
end

for i in 11:20
  @test t[i-10] == A[i]
end

A = LinearInterpolation{false}(u,t)
@test length(A) == 10
for i in 1:10
  @test u[i] == A[i]
end
