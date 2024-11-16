using DataInterpolations

@testset "Linear Interpolation" begin
    u = [1.0, 2.0]
    t = [1.0, 2.0]

    A = LinearInterpolation(u, t; extrapolation_down = ExtrapolationType.constant)
    t_eval = 0.0
    @test A(t_eval) == 1.0
    @test DataInterpolations.derivative(A, t_eval) == 0.0
    @test DataInterpolations.derivative(A, t_eval, 2) == 0.0
    @test DataInterpolations.integral(A, t_eval) == -1.0
    t_eval = 3.0
    @test_throws DataInterpolations.UpExtrapolationError A(t_eval)
    @test_throws DataInterpolations.UpExtrapolationError DataInterpolations.derivative(A, t_eval)
    @test_throws DataInterpolations.UpExtrapolationError DataInterpolations.derivative(A, t_eval, 2)
    @test_throws DataInterpolations.UpExtrapolationError DataInterpolations.integral(A, t_eval)

    A = LinearInterpolation(u, t; extrapolation_up = ExtrapolationType.constant)
    t_eval = 3.0
    @test A(t_eval) == 2.0
    @test DataInterpolations.derivative(A, t_eval) == 0.0
    @test DataInterpolations.derivative(A, t_eval, 2) == 0.0
    @test DataInterpolations.integral(A, t_eval) == 3.5
    t_eval = 0.0
    @test_throws DataInterpolations.DownExtrapolationError A(t_eval)
    @test_throws DataInterpolations.DownExtrapolationError DataInterpolations.derivative(A, t_eval)
    @test_throws DataInterpolations.DownExtrapolationError DataInterpolations.derivative(A, t_eval, 2)
    @test_throws DataInterpolations.DownExtrapolationError DataInterpolations.integral(A, t_eval)

    for extrapolation_type in [ExtrapolationType.linear, ExtrapolationType.extension]
        A = LinearInterpolation(u, t; extrapolation_down = extrapolation_type)
        t_eval = 0.0
        @test A(t_eval) == 0.0
        @test DataInterpolations.derivative(A, t_eval) == 1.0
        @test DataInterpolations.derivative(A, t_eval, 2) == 0.0
        @test DataInterpolations.integral(A, t_eval) == -0.5
        t_eval = 3.0
        @test_throws DataInterpolations.UpExtrapolationError A(t_eval)
        @test_throws DataInterpolations.UpExtrapolationError DataInterpolations.derivative(A, t_eval)
        @test_throws DataInterpolations.UpExtrapolationError DataInterpolations.derivative(A, t_eval, 2)
        @test_throws DataInterpolations.UpExtrapolationError DataInterpolations.integral(A, t_eval)

        A = LinearInterpolation(u, t; extrapolation_up = extrapolation_type)
        t_eval = 3.0
        @test A(t_eval) == 3.0
        @test DataInterpolations.derivative(A, t_eval) == 1.0
        @test DataInterpolations.derivative(A, t_eval, 2) == 0.0
        @test DataInterpolations.integral(A, t_eval) == 4.0
        t_eval = 0.0
        @test_throws DataInterpolations.DownExtrapolationError A(t_eval)
        @test_throws DataInterpolations.DownExtrapolationError DataInterpolations.derivative(A, t_eval)
        @test_throws DataInterpolations.DownExtrapolationError DataInterpolations.derivative(A, t_eval, 2)
        @test_throws DataInterpolations.DownExtrapolationError DataInterpolations.integral(A, t_eval)
    end
end

@testset "Quadratic Interpolation" begin
    u = [1.0, 3.0, 2.0]
    t = 1:3
    f = t -> (-3t^2 + 13t - 8)/2

    A = QuadraticInterpolation(u, t; extrapolation_down = ExtrapolationType.constant)
    t_eval = 0.0
    @test A(t_eval) ≈ 1.0 
    @test DataInterpolations.derivative(A, t_eval) == 0.0
    @test DataInterpolations.derivative(A, t_eval, 2) == 0.0
    @test DataInterpolations.integral(A, t_eval) ≈ -1.0
    t_eval = 4.0
    @test_throws DataInterpolations.UpExtrapolationError DataInterpolations.derivative(A, t_eval)
    @test_throws DataInterpolations.UpExtrapolationError DataInterpolations.derivative(A, t_eval, 2)
    @test_throws DataInterpolations.UpExtrapolationError DataInterpolations.integral(A, t_eval)

    A = QuadraticInterpolation(u, t; extrapolation_up = ExtrapolationType.constant)
    t_eval = 4.0
    @test A(t_eval) ≈ 2.0
    @test DataInterpolations.derivative(A, t_eval) == 0.0
    @test DataInterpolations.derivative(A, t_eval, 2) == 0.0
    @test DataInterpolations.integral(A, t[end], t_eval) ≈ 2.0
    t_eval = 0.0
    @test_throws DataInterpolations.DownExtrapolationError DataInterpolations.derivative(A, t_eval)
    @test_throws DataInterpolations.DownExtrapolationError DataInterpolations.derivative(A, t_eval, 2)
    @test_throws DataInterpolations.DownExtrapolationError DataInterpolations.integral(A, t_eval)
end