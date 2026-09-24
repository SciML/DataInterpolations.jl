using DataInterpolations

t1 = [1.0, 2.0, 3.0]
u1 = [0.0, 1.0, 0.0]

t2 = [4.0, 5.0, 6.0]
u2 = [1.0, 2.0, 1.0]

ts_append = 1.0:0.5:6.0
ts_push = 1.0:0.5:4.0

@testset "$method" for method in [
        LinearInterpolation, QuadraticInterpolation, ConstantInterpolation,
    ]
    func1 = method(copy(u1), copy(t1); cache_parameters = true)
    append!(func1, u2, t2)
    func2 = method(vcat(u1, u2), vcat(t1, t2); cache_parameters = true)
    @test func1.u == func2.u
    @test func1.t == func2.t
    for name in propertynames(func1.p)
        @test getfield(func1.p, name) == getfield(func2.p, name)
    end
    @test func1(ts_append) == func2(ts_append)
    @test func1.I == func2.I

    func1 = method(copy(u1), copy(t1); cache_parameters = true)
    push!(func1, 1.0, 4.0)
    func2 = method(vcat(u1, 1.0), vcat(t1, 4.0); cache_parameters = true)
    @test func1.u == func2.u
    @test func1.t == func2.t
    for name in propertynames(func1.p)
        @test getfield(func1.p, name) == getfield(func2.p, name)
    end
    @test func1(ts_push) == func2(ts_push)
    @test func1.I == func2.I
end
