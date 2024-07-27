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
    A = LinearInterpolation(u, t)
    @variables t x(t)
    substitute(A(t), Dict(t => x))
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
