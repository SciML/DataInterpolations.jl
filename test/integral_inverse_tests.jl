using DataInterpolations
using DataInterpolations: integral, invert_integral

function test_integral_inverses(method; args = [], kwargs = [])
    A = method(args...; kwargs...)
    @test hasfield(typeof(A), :I)
    A_intinv = invert_integral(A)
    @test A_intinv isa DataInterpolations.AbstractIntegralInverseInterpolation
    ts = range(first(A.t), last(A.t), length = 100)
    Is = integral.(Ref(A), ts)
    ts_ = A_intinv.(Is)
    @test ts ≈ ts_
end

@testset "Linear Interpolation" begin
    t = collect(1:5)
    u = [1.0, 1.0, 2.0, 4.0, 3.0]
    test_integral_inverses(LinearInterpolation; args = [u, t])

    u = [1.0, -1.0, 2.0, 4.0, 3.0]
    A = LinearInterpolation(u, t)
    @test_throws DataInterpolations.IntegralNotInvertibleError invert_integral(A)
end
