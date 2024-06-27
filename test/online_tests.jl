using DataInterpolations

t1 = Float64[1, 2, 3]
u1 = Float64[0, 1, 0]

t2 = Float64[4, 5, 6]
u2 = Float64[1, 2, 1]

ts = 1.0:0.5:6.0

for method in [LinearInterpolation, QuadraticInterpolation, ConstantInterpolation]
    func1 = method(u1, t1)
    append!(func1, u2, t2)
    func2 = method(vcat(u1, u2), vcat(t1, t2))
    @test func1.u == func2.u
    @test func1.t == func2.t
    for name in propertynames(func1.p)
       @test getfield(func1.p, name) == getfield(func2.p, name) 
    end
    @test func1(ts) == func2(ts)

    func1 = method(u1, t1)
    push!(func1, 1.0, 4.0)
    func2 = method(vcat(u, 1.0), vcat(t, 4.0))
    @test func1.u == func2.u
    @test func1.t == func2.t
    for name in propertynames(func1.p)
        @test getfield(func1.p, name) == getfield(func2.p, name) 
     end
    @test func1(ts) == func2(ts)
end
