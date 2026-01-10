using DataInterpolations
using AllocCheck: @check_allocs
using StaticArrays: SVector

@check_allocs(test_allocs(itp, x) = itp(x)) # Reuse function definition to save on compilation time

@testset "Allocation-free interpolation tests" begin
    @testset "LinearInterpolation" begin
        t = 1.0collect(1:10)
        x = 1:10
        y = 2:4
        u_ = x' .* y
        u = [u_[:, i] for i in 1:size(u_, 2)]
        u_s = [convert(SVector{length(y), eltype(u_)}, i) for i in u]
        A_s = LinearInterpolation(u_s, t; extrapolation = ExtrapolationType.Extension)
        @test_nowarn test_allocs(A_s, 0)
    end

    @testset "QuadraticInterpolation" begin
        t = [1.0, 2.0, 3.0, 4.0]
        u_ = [1.0, 4.0, 9.0, 16.0]' .* ones(5)
        u = [u_[:, i] for i in 1:size(u_, 2)]
        u_s = [convert(SVector{length(u[1])}, i) for i in u]
        A_s = QuadraticInterpolation(u_s, t; extrapolation = ExtrapolationType.Extension)
        @test_nowarn test_allocs(A_s, 0)
    end

    @testset "ConstantInterpolation" begin
        t = [1.0, 2.0, 3.0, 4.0]
        u = [[1.0, 2.0], [0.0, 1.0], [1.0, 2.0], [0.0, 1.0]]
        u_s = [convert(SVector{length(first(u))}, i) for i in u]
        A_s = ConstantInterpolation(u_s, t; extrapolation = ExtrapolationType.Extension)
        @test_nowarn test_allocs(A_s, 0)
    end

    @testset "SmoothedConstantInterpolation" begin
        u = [0.0, 2.0, 1.0, 3.0]
        t = [1.2, 2.5, 5.7, 8.7]
        d_max = 0.5
        u_ = u' .* ones(5)
        uv = [u_[:, i] for i in 1:size(u_, 2)]
        u_s = [convert(SVector{length(uv[1])}, i) for i in uv]
        A_s = SmoothedConstantInterpolation(u_s, t; d_max)
        @test_nowarn test_allocs(A_s, 1.9)
    end

    @testset "QuadraticSpline" begin
        t = [-1.0, 0.0, 1.0]
        u_ = [0.0, 1.0, 3.0]' .* ones(4)
        u = [u_[:, i] for i in 1:size(u_, 2)]
        u_s = [convert(SVector{length(u[1])}, i) for i in u]
        A_s = QuadraticSpline(u_s, t; extrapolation = ExtrapolationType.Extension)
        @test_nowarn test_allocs(A_s, 0.7)
    end

    @testset "CubicSpline" begin
        t = [-1.0, 0.0, 1.0]
        u_ = [0.0, 1.0, 3.0]' .* ones(4)
        u = [u_[:, i] for i in 1:size(u_, 2)]
        u_s = [convert(SVector{length(first(u))}, i) for i in u]
        A_s = CubicSpline(u_s, t; extrapolation = ExtrapolationType.Extension)
        @test_nowarn test_allocs(A_s, 0)
    end

    @testset "CubicHermiteSpline" begin
        du = [-0.047, -0.058, 0.054, 0.012, -0.068, 0.0]
        u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
        t = [0.0, 62.25, 109.66, 162.66, 205.8, 252.3]
        u2 = [[u[i], u[i] + 1] for i in eachindex(u)]
        du2 = [[du[i], du[i]] for i in eachindex(du)]
        u2_s = [convert(SVector{length(u2[1])}, i) for i in u2]
        du2_s = [convert(SVector{length(du2[1])}, i) for i in du2]
        A2_s = CubicHermiteSpline(
            du2_s, u2_s, t; extrapolation = ExtrapolationType.Extension
        )
        @test_nowarn test_allocs(A2_s, 0.7)
    end

    @testset "QuinticHermiteSpline" begin
        ddu = [0.0, -0.00033, 0.0051, -0.0067, 0.0029, 0.0]
        du = [-0.047, -0.058, 0.054, 0.012, -0.068, 0.0]
        u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
        t = [0.0, 62.25, 109.66, 162.66, 205.8, 252.3]
        u2 = [[u[i], u[i] + 1] for i in eachindex(u)]
        du2 = [[du[i], du[i]] for i in eachindex(du)]
        ddu2 = [[ddu[i], ddu[i]] for i in eachindex(ddu)]
        u2_s = [convert(SVector{length(u2[1])}, i) for i in u2]
        du2_s = [convert(SVector{length(du2[1])}, i) for i in du2]
        ddu2_s = [convert(SVector{length(du2[1])}, i) for i in ddu2]
        A2_s = QuinticHermiteSpline(
            ddu2_s, du2_s, u2_s, t; extrapolation = ExtrapolationType.Extension
        )
        @test_nowarn test_allocs(A2_s, 0.7)
    end
end
