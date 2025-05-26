using DataInterpolations
using DataInterpolations: integral, derivative, invert_integral

function test_integral_inverses(method; args = [], kwargs = [])
    A = method(args...; kwargs..., extrapolation = ExtrapolationType.Extension)
    @test hasfield(typeof(A), :I)
    A_intinv = invert_integral(A)
    @test A_intinv isa DataInterpolations.AbstractIntegralInverseInterpolation
    ts = range(first(A.t), last(A.t), length = 100)
    Is = integral.(Ref(A), ts)
    ts_ = A_intinv.(Is)
    @test ts ≈ ts_

    for I in Is
        cdiff = forward_fdm(5, 1; geom = true)(A_intinv, I)
        adiff = derivative(A_intinv, I)
        @test cdiff ≈ adiff
    end

    @test @inferred(A(ts[37])) == A(ts[37])
end

function test_integral_inverse_extrapolation()
    # Linear function with constant extrapolation
    t = collect(0:4)
    u = [0.0, 2.0, 3.0, 4.0]
    A = LinearInterpolation(u, t, extrapolation = ExtrapolationType.Constant)

    A_intinv = invert_integral(A, extrapolation_left = ExtrapolationType.Extension,
        extrapolation_right = ExtrapolationType.Extension)

    # for a linear function, the integral is quadratic
    # but the constant extrapolation part is linear.
    A_5 = 1 / 2 * 4.0^2 + 5.0
    @test A_int_inv(A_5) ≈ 5.0
end

@testset "Linear Interpolation" begin
    t = collect(1:5)
    u = [1.0, 1.0, 2.0, 4.0, 3.0]
    test_integral_inverses(LinearInterpolation; args = [u, t])

    u = [1.0, -1.0, 2.0, 4.0, 3.0]
    A = LinearInterpolation(u, t)
    @test_throws DataInterpolations.IntegralNotInvertibleError invert_integral(A)

    test_integral_inverse_extrapolation()
end

@testset "Constant Interpolation" begin
    t = collect(1:5)
    u = [1.0, 1.0, 2.0, 4.0, 3.0]
    test_integral_inverses(ConstantInterpolation; args = [u, t])
    test_integral_inverses(ConstantInterpolation; args = [u, t], kwargs = [:dir => :right])

    u = [1.0, -1.0, 2.0, 4.0, 3.0]
    A = ConstantInterpolation(u, t)
    @test_throws DataInterpolations.IntegralNotInvertibleError invert_integral(A)
end

t = collect(1:5)
u = [1.0, 1.0, 2.0, 4.0, 3.0]
A = QuadraticInterpolation(u, t)
@test_throws DataInterpolations.IntegralInverseNotFoundError invert_integral(A)
