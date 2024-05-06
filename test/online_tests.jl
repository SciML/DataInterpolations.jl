using DataInterpolations

t = [1, 2, 3]
u = [0, 1, 0]

for di in [LinearInterpolation, QuadraticInterpolation, ConstantInterpolation]
    li = di(copy(u), copy(t))
    append!(li, u, t)
    li2 = di(vcat(u, u), vcat(t, t))
    @test li.u == li2.u
    @test li.t == li2.t

    li = di(copy(u), copy(t))
    push!(li, 1, 4)
    li2 = di(vcat(u, 1), vcat(t, 4))
    @test li.u == li2.u
    @test li.t == li2.t
end
