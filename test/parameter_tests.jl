using DataInterpolations

@testset "Cubic Hermite Spline" begin
    du = [5.0, 3.0, 6.0, 8.0, 1.0]
    u = [1.0, 5.0, 3.0, 4.0, 4.0]
    t = collect(1:5)
    A = CubicHermiteSpline(du, u, t)
    @test A.p.c₁ ≈ [-1.0, -5.0, -5.0, -8.0]
    @test A.p.c₂ ≈ [0.0, 13.0, 12.0, 9.0]
end

@testset "Quintic Hermite Spline" begin
    ddu = [0.0, 3.0, 6.0, 4.0, 5.0]
    du = [5.0, 3.0, 6.0, 8.0, 1.0]
    u = [1.0, 5.0, 3.0, 4.0, 4.0]
    t = collect(1:5)
    A = QuinticHermiteSpline(ddu, du, u, t)
    @test A.p.c₁ ≈ [-1.0, -6.5, -8.0, -10.0]
    @test A.p.c₂ ≈ [1.0, 19.5, 20.0, 19.0]
    @test A.p.c₃ ≈ [1.5, -37.5, -37.0, -26.5]
end
