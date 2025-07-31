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
    A = LinearInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
    B = LinearInterpolation(u .^ 2, t; extrapolation = ExtrapolationType.Extension)
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

@testset "Output Type" begin
    # Test consistency between eltype(u) and type of the output
    u = Float32[-0.676367f0, 0.8449812f0, 1.2366607f0, -0.13347931f0, 1.9928657f0,
    -0.63596356f0, 0.76009744f0, -0.30632544f0, 0.34649512f0, -0.3846099f0]
    t = 0.1f0:0.1f0:1.0f0
    for extrapolation_flag in instances(ExtrapolationType.T)
        (extrapolation_flag == ExtrapolationType.None) && continue
        aki = AkimaInterpolation(u, t; extrapolation = extrapolation_flag)
        @test eltype(aki.([-2.0f0, 0.5f0, 3.0f0])) == Float32
    end
end
