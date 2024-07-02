using DataInterpolations

@testset "Linear Interpolation" begin
    u = [1.0, 5.0, 3.0, 4.0, 4.0]
    t = collect(1:5)
    A = LinearInterpolation(u, t)
    @test A.p.slope ≈ [4.0, -2.0, 1.0, 0.0]
end

@testset "Quadratic Interpolation" begin
    u = [1.0, 5.0, 3.0, 4.0, 4.0]
    t = collect(1:5)
    A = QuadraticInterpolation(u, t)
    @test A.p.l₀ ≈ [0.5, 2.5, 1.5]
    @test A.p.l₁ ≈ [-5.0, -3.0, -4.0]
    @test A.p.l₂ ≈ [1.5, 2.0, 2.0]
end

@testset "Quadratic Spline" begin
    u = [1.0, 5.0, 3.0, 4.0, 4.0]
    t = collect(1:5)
    A = QuadraticSpline(u, t)
    @test A.p.σ ≈ [4.0, -10.0, 13.0, -14.0]
end

@testset "Cubic Spline" begin
    u = [1, 5, 3, 4, 4]
    t = collect(1:5)
    A = CubicSpline(u, t)
    @test A.p.c₁ ≈ [6.839285714285714, 1.642857142857143, 4.589285714285714, 4.0]
    @test A.p.c₂ ≈ [1.0, 6.839285714285714, 1.642857142857143, 4.589285714285714]
end
