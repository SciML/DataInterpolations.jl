using DataInterpolations
using Symbolics

@testset "Interface" begin
    u = 2.0collect(1:10)
    t = 1.0collect(1:10)
    A = LinearInterpolation(u, t)

    for i in 1:10
        @test u[i] == A.u[i]
    end

    for i in 1:10
        @test t[i] == A.t[i]
    end
end

@testset "Symbolics" begin
    u = 2.0collect(1:10)
    t = 1.0collect(1:10)
    A = LinearInterpolation(u, t; extrapolation_up = ExtrapolationType.extension,
        extrapolation_down = ExtrapolationType.extension)
    B = LinearInterpolation(u .^ 2, t; extrapolation_up = ExtrapolationType.extension,
        extrapolation_down = ExtrapolationType.extension)
    @variables t x(t)
    substitute(A(t), Dict(t => x))
    t_val = 2.7
    @test substitute(A(t), Dict(t => t_val)) == A(t_val)
    @test substitute(B(A(t)), Dict(t => t_val)) == B(A(t_val))
    @test substitute(A(B(A(t))), Dict(t => t_val)) == A(B(A(t_val)))
end

@testset "Type Inference" begin
    u = 2.0collect(1:10)
    t = 1.0collect(1:10)
    methods = [
        ConstantInterpolation, LinearInterpolation,
        QuadraticInterpolation, LagrangeInterpolation,
        QuadraticSpline, CubicSpline, AkimaInterpolation
    ]
    @testset "$method" for method in methods
        @inferred method(u, t)
    end
    @testset "BSplineInterpolation" begin
        @inferred BSplineInterpolation(u, t, 3, :Uniform, :Uniform)
        @inferred BSplineInterpolation(u, t, 3, :ArcLen, :Average)
    end
    @testset "BSplineApprox" begin
        @inferred BSplineApprox(u, t, 3, 5, :Uniform, :Uniform)
        @inferred BSplineApprox(u, t, 3, 5, :ArcLen, :Average)
    end
    du = ones(10)
    ddu = zeros(10)
    @testset "Hermite Splines" begin
        @inferred CubicHermiteSpline(du, u, t)
        @inferred PCHIPInterpolation(u, t)
        @inferred QuinticHermiteSpline(ddu, du, u, t)
    end
end
