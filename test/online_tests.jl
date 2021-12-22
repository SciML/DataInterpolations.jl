using DataInterpolations, Test

t = [1,2,3]
u = [0, 1, 0]

for di in [LinearInterpolation, QuadraticInterpolation, ConstantInterpolation]
    li = di(copy(u), copy(t))
    append!(li, u, t)
    @test li == di(vcat(u, u), vcat(t, t))

    li = di(copy(u), copy(t))
    push!(li, 4, 1)
    @test li == di(vcat(u, 1), vcat(t, 4))
end