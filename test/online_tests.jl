using DataInterpolations, Test

t = [1,2,3]
u = [0, 1, 0]

for di in [LinearInterpolation, QuadraticInterpolation, ConstantInterpolation]
    li = LinearInterpolation(copy(u), copy(t))
    append!(li, u, t)
    @test li == LinearInterpolation(vcat(u, u), vcat(t, t))
end