using DataInterpolations

function test_extrapolation_errors(method, u, t)
    A = method(u, t)
    @test A.extrapolation_right == ExtrapolationType.none
    @test A.extrapolation_left == ExtrapolationType.none
    for (error_type, t_eval) in zip(
        (DataInterpolations.DownExtrapolationError,
            DataInterpolations.UpExtrapolationError),
        (first(t) - 1, last(t) + 1))
        @test_throws error_type A(t_eval)
        @test_throws error_type DataInterpolations.derivative(
            A, t_eval)
        @test_throws error_type DataInterpolations.derivative(
            A, t_eval, 2)
        @test_throws error_type DataInterpolations.integral(
            A, t_eval)
    end
end

function test_constant_extrapolation(method, u, t)
    A = method(u, t; extrapolation_left = ExtrapolationType.constant,
        extrapolation_right = ExtrapolationType.constant)
    t_lower = first(t) - 1
    t_upper = last(t) + 1
    @test A(t_lower) == first(u)
    @test A(t_upper) == last(u)
    @test DataInterpolations.derivative(A, t_lower) == 0
    @test DataInterpolations.derivative(A, t_upper) == 0
    @test DataInterpolations.integral(A, t_lower, first(t)) ≈
          first(u) * (first(t) - t_lower)
    @test DataInterpolations.integral(A, last(t), t_upper) ≈ last(u) * (t_upper - last(t))
end

@testset "Linear Interpolation" begin
    u = [1.0, 2.0]
    t = [1.0, 2.0]

    test_extrapolation_errors(LinearInterpolation, u, t)
    test_constant_extrapolation(LinearInterpolation, u, t)

    for extrapolation_type in [ExtrapolationType.linear, ExtrapolationType.extension]
        # Down extrapolation
        A = LinearInterpolation(u, t; extrapolation_left = extrapolation_type)
        t_eval = 0.0
        @test A(t_eval) == 0.0
        @test DataInterpolations.derivative(A, t_eval) == 1.0
        @test DataInterpolations.derivative(A, t_eval, 2) == 0.0
        @test DataInterpolations.integral(A, t_eval) == -0.5
        t_eval = 3.0

        # Up extrapolation
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

    test_extrapolation_errors(QuadraticInterpolation, u, t)
    test_constant_extrapolation(LinearInterpolation, u, t)

    # Linear down extrapolation
    A = QuadraticInterpolation(u, t; extrapolation_left = ExtrapolationType.linear)
    t_eval = 0.0
    @test A(t_eval) ≈ -2.5
    @test DataInterpolations.derivative(A, t_eval) == 3.5
    @test DataInterpolations.derivative(A, t_eval, 2) == 0.0
    @test DataInterpolations.integral(A, t_eval) ≈ 0.75

    # Linear up extrapolation
    A = QuadraticInterpolation(u, t; extrapolation_right = ExtrapolationType.linear)
    t_eval = 4.0
    @test A(t_eval) ≈ -0.5
    @test DataInterpolations.derivative(A, t_eval) == -2.5
    @test DataInterpolations.derivative(A, t_eval, 2) == 0.0
    @test DataInterpolations.integral(A, t[end], t_eval) ≈ 0.75

    # Extension down extrapolation
    f = t -> (-3t^2 + 13t - 8) / 2
    df = t -> (-6t + 13) / 2
    A = QuadraticInterpolation(u, t; extrapolation_left = ExtrapolationType.extension)
    t_eval = 0.0
    @test A(t_eval) ≈ -4.0
    @test DataInterpolations.derivative(A, t_eval) == df(t_eval)
    @test DataInterpolations.derivative(A, t_eval, 2) == -3
    @test DataInterpolations.integral(A, t_eval) ≈ 1.25

    # Extension up extrapolation
    A = QuadraticInterpolation(u, t; extrapolation_right = ExtrapolationType.extension)
    t_eval = 4.0
    @test A(t_eval) ≈ -2.0
    @test DataInterpolations.derivative(A, t_eval) == df(t_eval)
    @test DataInterpolations.derivative(A, t_eval, 2) == -3
    @test DataInterpolations.integral(A, t_eval) ≈ 5.25
end
