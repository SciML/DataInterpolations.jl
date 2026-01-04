using DataInterpolations, Test
using ForwardDiff
using QuadGK
using Unitful

function test_extrapolation(method, u, t)
    @testset "Extrapolation errors" begin
        A = method(u, t)
        @test A.extrapolation_right == ExtrapolationType.None
        @test A.extrapolation_left == ExtrapolationType.None
        for (error_type, t_eval) in zip(
                (
                    DataInterpolations.LeftExtrapolationError,
                    DataInterpolations.RightExtrapolationError,
                ),
                (first(t) - 1, last(t) + 1)
            )
            @test_throws error_type A(t_eval)
            @test_throws error_type DataInterpolations.derivative(
                A, t_eval
            )
            @test_throws error_type DataInterpolations.derivative(
                A, t_eval, 2
            )
            @test_throws error_type DataInterpolations.integral(
                A, t_eval
            )
        end
    end

    for extrapolation_type in instances(ExtrapolationType.T)
        (extrapolation_type == ExtrapolationType.None) && continue
        @testset "extrapolation type $extrapolation_type" begin
            A = method(u, t; extrapolation = extrapolation_type)

            t_eval = first(t) - 1.5
            @test DataInterpolations.derivative(A, t_eval) ≈
                ForwardDiff.derivative(A, t_eval)

            t_eval = last(t) + 1.5
            @test DataInterpolations.derivative(A, t_eval) ≈
                ForwardDiff.derivative(A, t_eval)

            T = last(A.t) - first(A.t)
            t1 = first(t) - 2.5T
            t2 = last(t) + 3.5T
            @test DataInterpolations.integral(A, t1, t2) ≈
                quadgk(A, t1, t2; atol = 1.0e-12, rtol = 1.0e-12)[1]
        end
    end
    return
end

@testset "Constant Interpolation with Unitful" begin
    t_un = [1.0, 2.0]u"s"
    u_un = [1.0, 2.0]u"m"

    for extrapolation_type in [ExtrapolationType.Constant, ExtrapolationType.Linear]
        # Left extrapolation
        A = ConstantInterpolation(u_un, t_un; extrapolation_left = extrapolation_type)
        t_eval = 0.0u"s"
        @test @inferred(A(t_eval)) == 1.0u"m"
        @test @inferred(A([t_eval])) == [1.0u"m"]
        @test A([t_eval]) isa Vector{typeof(1.0u"m")}

        # Right extrapolation
        A = ConstantInterpolation(u_un, t_un; extrapolation_right = extrapolation_type)
        t_eval = 3.0u"s"
        @test @inferred(A(t_eval)) == 2.0u"m"
        @test @inferred(A([t_eval])) == [2.0u"m"]
        @test A([t_eval]) isa Vector{typeof(2.0u"m")}
    end
end

@testset "Linear Interpolation with Unitful" begin
    t_un = [1.0, 2.0]u"s"
    u_un = [1.0, 2.0]u"m"

    # Left constant extrapolation
    A = LinearInterpolation(u_un, t_un; extrapolation_left = ExtrapolationType.Constant)
    t_eval = 0.0u"s"
    @test @inferred(A(t_eval)) == 1.0u"m"
    @test @inferred(A([t_eval])) == [1.0u"m"]
    @test A([t_eval]) isa Vector{typeof(1.0u"m")}

    # Right constant extrapolation
    A = LinearInterpolation(u_un, t_un; extrapolation_right = ExtrapolationType.Constant)
    t_eval = 3.0u"s"
    @test @inferred(A(t_eval)) == 2.0u"m"
    @test @inferred(A([t_eval])) == [2.0u"m"]
    @test A([t_eval]) isa Vector{typeof(2.0u"m")}

    # Left linear extrapolation
    A = LinearInterpolation(u_un, t_un; extrapolation_left = ExtrapolationType.Linear)
    t_eval = 0.0u"s"
    @test @inferred(A(t_eval)) == 0.0u"m"
    @test @inferred(A([t_eval])) == [0.0u"m"]
    @test A([t_eval]) isa Vector{typeof(0.0u"m")}

    # Right constant extrapolation
    A = LinearInterpolation(u_un, t_un; extrapolation_right = ExtrapolationType.Linear)
    t_eval = 3.0u"s"
    @test @inferred(A(t_eval)) == 3.0u"m"
    @test @inferred(A([t_eval])) == [3.0u"m"]
    @test A([t_eval]) isa Vector{typeof(3.0u"m")}
end

@testset "Linear Interpolation" begin
    u = [1.0, 2.0]
    t = [1.0, 2.0]

    test_extrapolation(LinearInterpolation, u, t)

    for extrapolation_type in [ExtrapolationType.Linear, ExtrapolationType.Extension]
        # Left extrapolation
        A = LinearInterpolation(u, t; extrapolation_left = extrapolation_type)
        t_eval = 0.0
        @test A(t_eval) == 0.0
        @test DataInterpolations.derivative(A, t_eval) == 1.0
        @test DataInterpolations.derivative(A, t_eval, 2) == 0.0
        @test DataInterpolations.integral(A, t_eval) == -0.5
        t_eval = 3.0

        # Right extrapolation
        A = LinearInterpolation(u, t; extrapolation_right = extrapolation_type)
        t_eval = 3.0
        @test A(t_eval) == 3.0
        @test DataInterpolations.derivative(A, t_eval) == 1.0
        @test DataInterpolations.derivative(A, t_eval, 2) == 0.0
        @test DataInterpolations.integral(A, t_eval) == 4.0
        t_eval = 0.0
    end
end

@testset "Quadratic Interpolation" begin
    u = [1.0, 3.0, 2.0]
    t = 1:3

    test_extrapolation(QuadraticInterpolation, u, t)

    # Linear left extrapolation
    A = QuadraticInterpolation(u, t; extrapolation_left = ExtrapolationType.Linear)
    t_eval = 0.0
    @test A(t_eval) ≈ -2.5
    @test DataInterpolations.derivative(A, t_eval) == 3.5
    @test DataInterpolations.derivative(A, t_eval, 2) == 0.0
    @test DataInterpolations.integral(A, t_eval) ≈ 0.75

    # Linear right extrapolation
    A = QuadraticInterpolation(u, t; extrapolation_right = ExtrapolationType.Linear)
    t_eval = 4.0
    @test A(t_eval) ≈ -0.5
    @test DataInterpolations.derivative(A, t_eval) == -2.5
    @test DataInterpolations.derivative(A, t_eval, 2) == 0.0
    @test DataInterpolations.integral(A, t[end], t_eval) ≈ 0.75

    # Extension left extrapolation
    f = t -> (-3t^2 + 13t - 8) / 2
    df = t -> (-6t + 13) / 2
    A = QuadraticInterpolation(u, t; extrapolation_left = ExtrapolationType.Extension)
    t_eval = 0.0
    @test A(t_eval) ≈ -4.0
    @test DataInterpolations.derivative(A, t_eval) == df(t_eval)
    @test DataInterpolations.derivative(A, t_eval, 2) == -3
    @test DataInterpolations.integral(A, t_eval) ≈ 1.25

    # Extension right extrapolation
    A = QuadraticInterpolation(u, t; extrapolation_right = ExtrapolationType.Extension)
    t_eval = 4.0
    @test A(t_eval) ≈ -2.0
    @test DataInterpolations.derivative(A, t_eval) == df(t_eval)
    @test DataInterpolations.derivative(A, t_eval, 2) == -3
    @test DataInterpolations.integral(A, t_eval) ≈ 5.25
end
