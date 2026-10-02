using DataInterpolations
using StaticArrays: SVector, SA

function test_cached_integration(method, args...)
    A_c = method(args...; cache_parameters = true)
    A_nc = method(args...; cache_parameters = false)
    return @test DataInterpolations.integral(A_c, last(A_c.t)) ≈
        DataInterpolations.integral(A_nc, last(A_nc.t))
end

@testset "Linear Interpolation" begin
    u = [1.0, 5.0, 3.0, 4.0, 4.0]
    t = collect(1:5)
    A = LinearInterpolation(u, t; cache_parameters = true)
    @test A.p.slope ≈ [4.0, -2.0, 1.0, 0.0]
    test_cached_integration(LinearInterpolation, u, t)
end

@testset "Smoothed constant Interpolation" begin
    u = [1.0, 5.0, 3.0, 4.0, 4.0]
    t = collect(1:5)
    A = SmoothedConstantInterpolation(u, t; cache_parameters = true)
    @test A.p.d ≈ [0.5, 0.5, 0.5, 0.5, 0.5]
    @test A.p.c ≈ [0.0, 2.0, -1.0, 0.5, 0.0]
end

@testset "Quadratic Interpolation" begin
    u = [1.0, 5.0, 3.0, 4.0, 4.0]
    t = collect(1:5)
    A = QuadraticInterpolation(u, t; cache_parameters = true)
    @test A.p.α ≈ [-3.0, 1.5, -0.5, -0.5]
    @test A.p.β ≈ [7.0, -3.5, 1.5, 0.5]
    test_cached_integration(QuadraticInterpolation, u, t)
end

@testset "Quadratic Spline" begin
    u = [1.0, 5.0, 3.0, 4.0, 4.0]
    t = collect(1:5)
    A = QuadraticSpline(u, t; cache_parameters = true)
    @test A.p.α ≈ [-9.5, 3.5, -0.5, -0.5]
    @test A.p.β ≈ [13.5, -5.5, 1.5, 0.5]
    test_cached_integration(QuadraticSpline, u, t)
end

@testset "Cubic Spline" begin
    u = [1, 5, 3, 4, 4]
    t = collect(1:5)
    A = CubicSpline(u, t; cache_parameters = true)
    @test A.p.c₁ ≈ [6.839285714285714, 1.642857142857143, 4.589285714285714, 4.0]
    @test A.p.c₂ ≈ [1.0, 6.839285714285714, 1.642857142857143, 4.589285714285714]
    test_cached_integration(CubicSpline, u, t)
end

@testset "Cubic Hermite Spline" begin
    du = [5.0, 3.0, 6.0, 8.0, 1.0]
    u = [1.0, 5.0, 3.0, 4.0, 4.0]
    t = collect(1:5)
    A = CubicHermiteSpline(du, u, t; cache_parameters = true)
    @test A.p.c₁ ≈ [-1.0, -5.0, -5.0, -8.0]
    @test A.p.c₂ ≈ [0.0, 13.0, 12.0, 9.0]
    test_cached_integration(CubicHermiteSpline, du, u, t)
end

@testset "Quintic Hermite Spline" begin
    ddu = [0.0, 3.0, 6.0, 4.0, 5.0]
    du = [5.0, 3.0, 6.0, 8.0, 1.0]
    u = [1.0, 5.0, 3.0, 4.0, 4.0]
    t = collect(1:5)
    A = QuinticHermiteSpline(ddu, du, u, t; cache_parameters = true)
    @test A.p.c₁ ≈ [-1.0, -6.5, -8.0, -10.0]
    @test A.p.c₂ ≈ [1.0, 19.5, 20.0, 19.0]
    @test A.p.c₃ ≈ [1.5, -37.5, -37.0, -26.5]
    test_cached_integration(QuinticHermiteSpline, ddu, du, u, t)
end

@testset "StaticArray parameter and integral caches" begin
    u = SA[1.0, 2.0, 4.0]
    t = SA[0.0, 1.0, 2.0]
    u_v = [1.0, 2.0, 4.0]
    t_v = [0.0, 1.0, 2.0]

    for cp in (false, true)
        A = LinearInterpolation(u, t; cache_parameters = cp)
        n = cp ? 2 : 0
        @test A.p.slope isa SVector{n, Float64}
        @test A.I isa SVector{n, Float64}
        A_v = LinearInterpolation(u_v, t_v; cache_parameters = cp)
        @test A(1.5) ≈ A_v(1.5)
        @test DataInterpolations.integral(A, 0.0, 2.0) ≈
            DataInterpolations.integral(A_v, 0.0, 2.0)
        @test A_v.p.slope isa Vector{Float64}
        @test A_v.I isa Vector{Float64}
    end

    A = LinearInterpolation(u, t; cache_parameters = true)
    @test A.p.slope ≈ [1.0, 2.0]

    u4 = SA[1.0, 5.0, 3.0, 4.0]
    t4 = SA[1.0, 2.0, 3.0, 4.0]
    for (method, args) in (
            (QuadraticInterpolation, (u4, t4)),
            (QuadraticSpline, (u4, t4)),
            (CubicSpline, (u4, t4)),
            (SmoothedConstantInterpolation, (u4, t4)),
            (CubicHermiteSpline, (SA[5.0, 3.0, 6.0, 8.0], u4, t4)),
            (QuinticHermiteSpline, (SA[0.0, 3.0, 6.0, 4.0], SA[5.0, 3.0, 6.0, 8.0], u4, t4)),
        )
        for cp in (false, true)
            A = method(args...; cache_parameters = cp)
            fields = filter(n -> isdefined(A.p, n), (:slope, :α, :β, :c₁, :c₂, :c₃, :d, :c))
            for f in fields
                @test getfield(A.p, f) isa SVector
            end
            @test A.I isa SVector
            @test A(2.5) ≈ method(collect.(args)...; cache_parameters = cp)(2.5)
        end
    end
end
