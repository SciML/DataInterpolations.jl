using DataInterpolations, Test
using StableRNGs
using Optim
import ForwardDiff

@testset "Linear Interpolation" begin
    for t in (1.0:10.0, 1.0collect(1:10))
        u = 2.0collect(1:10)
        #t = 1.0collect(1:10)
        A = LinearInterpolation(u, t; extrapolate = true)

        for (_t, _u) in zip(t, u)
            @test A(_t) == _u
        end
        @test A(0) == 0.0
        @test A(5.5) == 11.0
        @test A(11) == 22

        u = vcat(2.0collect(1:10)', 3.0collect(1:10)')
        A = LinearInterpolation(u, t; extrapolate = true)

        for (_t, _u) in zip(t, eachcol(u))
            @test A(_t) == _u
        end
        @test A(0) == [0.0, 0.0]
        @test A(5.5) == [11.0, 16.5]
        @test A(11) == [22, 33]

        x = 1:10
        y = 2:4
        u_ = x' .* y
        u = [u_[:, i] for i in 1:size(u_, 2)]
        A = LinearInterpolation(u, t; extrapolate = true)
        @test A(0) == [0.0, 0.0, 0.0]
        @test A(5.5) == [11.0, 16.5, 22.0]
        @test A(11) == [22.0, 33.0, 44.0]
    end

    x = 1:10
    y = 2:4
    u_ = x' .* y
    u = [u_[:, i:(i + 1)] for i in 1:2:10]
    t = 1.0collect(2:2:10)
    A = LinearInterpolation(u, t; extrapolate = true)

    @test A(0) == [-2.0 0.0; -3.0 0.0; -4.0 0.0]
    @test A(3) == [4.0 6.0; 6.0 9.0; 8.0 12.0]
    @test A(5) == [8.0 10.0; 12.0 15.0; 16.0 20.0]

    # with NaNs (#113)
    u = [NaN, 1.0, 2.0, 3.0]
    t = 1:4
    A = LinearInterpolation(u, t; extrapolate = true)
    @test isnan(A(1.0))
    @test A(2.0) == 1.0
    @test A(2.5) == 1.5
    @test A(3.0) == 2.0
    @test A(4.0) == 3.0

    u = [0.0, NaN, 2.0, 3.0]
    A = LinearInterpolation(u, t; extrapolate = true)
    @test A(1.0) == 0.0
    @test isnan(A(2.0))
    @test isnan(A(2.5))
    @test A(3.0) == 2.0
    @test A(4.0) == 3.0

    u = [0.0, NaN, 2.0, 3.0]
    A = LinearInterpolation(u, t; extrapolate = true)
    @test A(1.0) == 0.0
    @test isnan(A(2.0))
    @test isnan(A(2.5))
    @test A(3.0) == 2.0
    @test A(4.0) == 3.0

    u = [0.0, 1.0, NaN, 3.0]
    A = LinearInterpolation(u, t; extrapolate = true)
    @test A(1.0) == 0.0
    @test A(2.0) == 1.0
    @test isnan(A(2.5))
    @test isnan(A(3.0))
    @test A(4.0) == 3.0

    u = [0.0, 1.0, 2.0, NaN]
    A = LinearInterpolation(u, t; extrapolate = true)
    @test A(1.0) == 0.0
    @test A(2.0) == 1.0
    @test A(3.0) == 2.0
    @test isnan(A(3.5))
    @test isnan(A(4.0))

    # Test type stability
    u = Float32.(1:5)
    t = Float32.(1:5)
    A1 = LinearInterpolation(u, t; extrapolate = true)
    u = 1:5
    t = 1:5
    A2 = LinearInterpolation(u, t; extrapolate = true)
    u = [1 // i for i in 1:5]
    t = (1:5)
    A3 = LinearInterpolation(u, t; extrapolate = true)
    u = [1 // i for i in 1:5]
    t = [1 // (6 - i) for i in 1:5]
    A4 = LinearInterpolation(u, t; extrapolate = true)

    F32 = Float32(1)
    F64 = Float64(1)
    I32 = Int32(1)
    I64 = Int64(1)
    R32 = Int32(1) // Int32(1)
    R64 = 1 // 1
    for A in Any[A1, A2, A3, A4]
        @test @inferred(A(F32)) === A(F32)
        @test @inferred(A(F64)) === A(F64)
        @test @inferred(A(I32)) === A(I32)
        @test @inferred(A(I64)) === A(I64)
        @test @inferred(A(R32)) === A(R32)
        @test @inferred(A(R64)) === A(R64)
    end

    # Nan time value:
    t = 0.0:3  # Floats
    u = [0, -2, -1, -2]
    A = LinearInterpolation(u, t; extrapolate = true)
    dA = t -> ForwardDiff.derivative(A, t)
    @test isnan(dA(NaN))

    t = 0:3  # Integers
    u = [0, -2, -1, -2]
    A = LinearInterpolation(u, t; extrapolate = true)
    dA = t -> ForwardDiff.derivative(A, t)
    @test isnan(dA(NaN))

    # Test derivative at point gives derivative to the right (except last is to left):
    ts = t[begin:(end - 1)]
    @test dA.(ts) == dA.(ts .+ 0.5)
    # Test last derivative is to the left:
    @test dA(last(t)) == dA(last(t) - 0.5)

    # Test array-valued interpolation
    u = collect.(2.0collect(1:10))
    t = 1.0collect(1:10)
    A = LinearInterpolation(u, t; extrapolate = true)
    @test A(0) == fill(0.0)
    @test A(5.5) == fill(11.0)
    @test A(11) == fill(22)

    # Test extrapolation
    u = 2.0collect(1:10)
    t = 1.0collect(1:10)
    A = LinearInterpolation(u, t; extrapolate = true)
    @test A(-1.0) == -2.0
    @test A(11.0) == 22.0
    A = LinearInterpolation(u, t)
    @test_throws DataInterpolations.ExtrapolationError A(-1.0)
    @test_throws DataInterpolations.ExtrapolationError A(11.0)
end

@testset "Quadratic Interpolation" begin
    u = [1.0, 4.0, 9.0, 16.0]
    t = [1.0, 2.0, 3.0, 4.0]
    A = QuadraticInterpolation(u, t; extrapolate = true)

    for (_t, _u) in zip(t, u)
        @test A(_t) == _u
    end
    @test A(0.0) == 0.0
    @test A(1.5) == 2.25
    @test A(2.5) == 6.25
    @test A(3.5) == 12.25
    @test A(5.0) == 25

    # backward-looking interpolation
    u = [1.0, 4.0, 9.0, 16.0]
    t = [1.0, 2.0, 3.0, 4.0]
    A = QuadraticInterpolation(u, t, :Backward; extrapolate = true)

    for (_t, _u) in zip(t, u)
        @test A(_t) == _u
    end
    @test A(0.0) == 0.0
    @test A(1.5) == 2.25
    @test A(2.5) == 6.25
    @test A(3.5) == 12.25
    @test A(5.0) == 25

    # Test both forward and backward-looking quadratic interpolation
    u = [1.0, 4.5, 6.0, 2.0]
    t = [1.0, 2.0, 3.0, 4.0]
    A_f = QuadraticInterpolation(u, t, :Forward)
    A_b = QuadraticInterpolation(u, t, :Backward)

    for (_t, _u) in zip(t, u)
        @test A_f(_t) == _u
        @test A_b(_t) == _u
    end
    l₀, l₁, l₂ = 0.375, 0.75, -0.125
    # In the first subinterval they're the same (no other option)
    @test A_f(1.5) ≈ l₀ * u[1] + l₁ * u[2] + l₂ * u[3]
    @test A_b(1.5) ≈ l₀ * u[1] + l₁ * u[2] + l₂ * u[3]
    # In the second subinterval they should be different
    @test A_f(2.5) ≈ l₀ * u[2] + l₁ * u[3] + l₂ * u[4]
    @test A_b(2.5) ≈ l₂ * u[1] + l₁ * u[2] + l₀ * u[3]
    # In the last subinterval they should be the same again
    @test A_f(3.5) ≈ l₂ * u[2] + l₁ * u[3] + l₀ * u[4]
    @test A_b(3.5) ≈ l₂ * u[2] + l₁ * u[3] + l₀ * u[4]

    # Matrix interpolation test
    u = [1.0 4.0 9.0 16.0; 1.0 4.0 9.0 16.0]
    A = QuadraticInterpolation(u, t; extrapolate = true)

    for (_t, _u) in zip(t, eachcol(u))
        @test A(_t) == _u
    end
    @test A(0.0) == [0.0, 0.0]
    @test A(1.5) == [2.25, 2.25]
    @test A(2.5) == [6.25, 6.25]
    @test A(3.5) == [12.25, 12.25]
    @test A(5.0) == [25.0, 25.0]

    u_ = [1.0, 4.0, 9.0, 16.0]' .* ones(5)
    u = [u_[:, i] for i in 1:size(u_, 2)]
    A = QuadraticInterpolation(u, t; extrapolate = true)
    @test A(0) == zeros(5)
    @test A(1.5) == 2.25 * ones(5)
    @test A(2.5) == 6.25 * ones(5)
    @test A(3.5) == 12.25 * ones(5)
    @test A(5.0) == 25.0 * ones(5)

    u = [repeat(u[i], 1, 3) for i in 1:4]
    A = QuadraticInterpolation(u, t; extrapolate = true)
    @test A(0) == zeros(5, 3)
    @test A(1.5) == 2.25 * ones(5, 3)
    @test A(2.5) == 6.25 * ones(5, 3)
    @test A(3.5) == 12.25 * ones(5, 3)
    @test A(5.0) == 25.0 * ones(5, 3)

    # Test extrapolation
    u = [1.0, 4.5, 6.0, 2.0]
    t = [1.0, 2.0, 3.0, 4.0]
    A = QuadraticInterpolation(u, t; extrapolate = true)
    @test A(0.0) == -4.5
    @test A(5.0) == -7.5
    A = QuadraticInterpolation(u, t)
    @test_throws DataInterpolations.ExtrapolationError A(0.0)
    @test_throws DataInterpolations.ExtrapolationError A(5.0)
end

@testset "Lagrange Interpolation" begin
    u = [1.0, 4.0, 9.0]
    t = [1.0, 2.0, 3.0]
    A = LagrangeInterpolation(u, t)

    @test A(2.0) == 4.0
    @test A(1.5) == 2.25

    u = [1.0, 8.0, 27.0, 64.0]
    t = [1.0, 2.0, 3.0, 4.0]
    A = LagrangeInterpolation(u, t)

    @test A(2.0) == 8.0
    @test A(1.5) ≈ 3.375
    @test A(3.5) ≈ 42.875

    u = [1.0 4.0 9.0 16.0; 1.0 4.0 9.0 16.0]
    A = LagrangeInterpolation(u, t)

    @test A(2.0) == [4.0, 4.0]
    @test A(1.5) ≈ [2.25, 2.25]
    @test A(3.5) ≈ [12.25, 12.25]

    u_ = [1.0, 4.0, 9.0]' .* ones(4)
    u = [u_[:, i] for i in 1:size(u_, 2)]
    t = [1.0, 2.0, 3.0]
    A = LagrangeInterpolation(u, t)

    @test A(2.0) == 4.0 * ones(4)
    @test A(1.5) == 2.25 * ones(4)

    u_ = [1.0, 8.0, 27.0, 64.0]' .* ones(4)
    u = [u_[:, i] for i in 1:size(u_, 2)]
    t = [1.0, 2.0, 3.0, 4.0]
    A = LagrangeInterpolation(u, t)

    @test A(2.0) == 8.0 * ones(4)
    @test A(1.5) ≈ 3.375 * ones(4)
    @test A(3.5) ≈ 42.875 * ones(4)

    u = [repeat(u[i], 1, 3) for i in 1:4]
    A = LagrangeInterpolation(u, t)

    @test A(2.0) == 8.0 * ones(4, 3)
    @test A(1.5) ≈ 3.375 * ones(4, 3)
    @test A(3.5) ≈ 42.875 * ones(4, 3)

    # Test extrapolation
    u = [1.0, 4.0, 9.0]
    t = [1.0, 2.0, 3.0]
    A = LagrangeInterpolation(u, t; extrapolate = true)
    @test A(0.0) == 0.0
    @test A(4.0) == 16.0
    A = LagrangeInterpolation(u, t)
    @test_throws DataInterpolations.ExtrapolationError A(-1.0)
    @test_throws DataInterpolations.ExtrapolationError A(4.0)
end

@testset "Akima Interpolation" begin
    u = [0.0, 2.0, 1.0, 3.0, 2.0, 6.0, 5.5, 5.5, 2.7, 5.1, 3.0]
    t = collect(0.0:10.0)
    A = AkimaInterpolation(u, t)

    @test A(0.0) ≈ 0.0
    @test A(0.5) ≈ 1.375
    @test A(1.0) ≈ 2.0
    @test A(1.5) ≈ 1.5
    @test A(2.5) ≈ 1.953125
    @test A(3.5) ≈ 2.484375
    @test A(4.5) ≈ 4.1363636363636366866103344
    @test A(5.1) ≈ 5.9803623910336236590978842
    @test A(6.5) ≈ 5.5067291516462386624652936
    @test A(7.2) ≈ 5.2031367459745245795943447
    @test A(8.6) ≈ 4.1796554159017080820603951
    @test A(9.9) ≈ 3.4110386597938129327189927
    @test A(10.0) ≈ 3.0

    # Test extrapolation
    A = AkimaInterpolation(u, t; extrapolate = true)
    @test A(-1.0) == 0.0
    @test A(11.0) == 3.0
    A = AkimaInterpolation(u, t)
    @test_throws DataInterpolations.ExtrapolationError A(-1.0)
    @test_throws DataInterpolations.ExtrapolationError A(11.0)
end

@testset "ConstantInterpolation" begin
    t = [1.0, 2.0, 3.0, 4.0]

    @testset "Vector case" for u in [[1.0, 2.0, 0.0, 1.0], ["B", "C", "A", "B"]]
        A = ConstantInterpolation(u, t, dir = :right; extrapolate = true)
        @test A(0.5) == u[1]
        @test A(1.0) == u[1]
        @test A(1.5) == u[2]
        @test A(2.0) == u[2]
        @test A(2.5) == u[3]
        @test A(3.0) == u[3]
        @test A(3.5) == u[1]
        @test A(4.0) == u[1]
        @test A(4.5) == u[1]

        A = ConstantInterpolation(u, t; extrapolate = true) # dir=:left is default
        @test A(0.5) == u[1]
        @test A(1.0) == u[1]
        @test A(1.5) == u[1]
        @test A(2.0) == u[2]
        @test A(2.5) == u[2]
        @test A(3.0) == u[3]
        @test A(3.5) == u[3]
        @test A(4.0) == u[1]
        @test A(4.5) == u[1]
    end

    @testset "Matrix case" for u in [
        [1.0 2.0 0.0 1.0; 1.0 2.0 0.0 1.0],
        ["B" "C" "A" "B"; "B" "C" "A" "B"],
    ]
        A = ConstantInterpolation(u, t, dir = :right; extrapolate = true)
        @test A(0.5) == u[:, 1]
        @test A(1.0) == u[:, 1]
        @test A(1.5) == u[:, 2]
        @test A(2.0) == u[:, 2]
        @test A(2.5) == u[:, 3]
        @test A(3.0) == u[:, 3]
        @test A(3.5) == u[:, 1]
        @test A(4.0) == u[:, 1]
        @test A(4.5) == u[:, 1]

        A = ConstantInterpolation(u, t; extrapolate = true) # dir=:left is default
        @test A(0.5) == u[:, 1]
        @test A(1.0) == u[:, 1]
        @test A(1.5) == u[:, 1]
        @test A(2.0) == u[:, 2]
        @test A(2.5) == u[:, 2]
        @test A(3.0) == u[:, 3]
        @test A(3.5) == u[:, 3]
        @test A(4.0) == u[:, 1]
        @test A(4.5) == u[:, 1]
    end

    @testset "Vector of Vectors case" for u in [
        [[1.0, 2.0], [0.0, 1.0], [1.0, 2.0], [0.0, 1.0]],
        [["B", "C"], ["A", "B"], ["B", "C"], ["A", "B"]]]
        A = ConstantInterpolation(u, t, dir = :right; extrapolate = true)
        @test A(0.5) == u[1]
        @test A(1.0) == u[1]
        @test A(1.5) == u[2]
        @test A(2.0) == u[2]
        @test A(2.5) == u[3]
        @test A(3.0) == u[3]
        @test A(3.5) == u[4]
        @test A(4.0) == u[4]
        @test A(4.5) == u[4]

        A = ConstantInterpolation(u, t; extrapolate = true) # dir=:left is default
        @test A(0.5) == u[1]
        @test A(1.0) == u[1]
        @test A(1.5) == u[1]
        @test A(2.0) == u[2]
        @test A(2.5) == u[2]
        @test A(3.0) == u[3]
        @test A(3.5) == u[3]
        @test A(4.0) == u[4]
        @test A(4.5) == u[4]
    end

    @testset "Vector of Matrices case" for u in [
        [[1.0 2.0; 1.0 2.0], [0.0 1.0; 0.0 1.0], [1.0 2.0; 1.0 2.0], [0.0 1.0; 0.0 1.0]],
        [["B" "C"; "B" "C"], ["A" "B"; "A" "B"], ["B" "C"; "B" "C"], ["A" "B"; "A" "B"]]]
        A = ConstantInterpolation(u, t, dir = :right; extrapolate = true)
        @test A(0.5) == u[1]
        @test A(1.0) == u[1]
        @test A(1.5) == u[2]
        @test A(2.0) == u[2]
        @test A(2.5) == u[3]
        @test A(3.0) == u[3]
        @test A(3.5) == u[4]
        @test A(4.0) == u[4]
        @test A(4.5) == u[4]

        A = ConstantInterpolation(u, t; extrapolate = true) # dir=:left is default
        @test A(0.5) == u[1]
        @test A(1.0) == u[1]
        @test A(1.5) == u[1]
        @test A(2.0) == u[2]
        @test A(2.5) == u[2]
        @test A(3.0) == u[3]
        @test A(3.5) == u[3]
        @test A(4.0) == u[4]
        @test A(4.5) == u[4]
    end

    # Test extrapolation
    u = [1.0, 2.0, 0.0, 1.0]
    A = ConstantInterpolation(u, t; extrapolate = true)
    @test A(-1.0) == 1.0
    @test A(11.0) == 1.0
    A = ConstantInterpolation(u, t)
    @test_throws DataInterpolations.ExtrapolationError A(-1.0)
    @test_throws DataInterpolations.ExtrapolationError A(11.0)
end

@testset "QuadraticSpline Interpolation" begin
    u = [0.0, 1.0, 3.0]
    t = [-1.0, 0.0, 1.0]

    A = QuadraticSpline(u, t; extrapolate = true)

    # Solution
    P₁ = x -> (x + 1)^2 # for x ∈ [-1, 0]
    P₂ = x -> 2 * x + 1   # for x ∈ [ 0, 1]

    for (_t, _u) in zip(t, u)
        @test A(_t) == _u
    end
    @test A(-2.0) == P₁(-2.0)
    @test A(-0.5) == P₁(-0.5)
    @test A(0.7) == P₂(0.7)
    @test A(2.0) == P₂(2.0)

    u_ = [0.0, 1.0, 3.0]' .* ones(4)
    u = [u_[:, i] for i in 1:size(u_, 2)]
    A = QuadraticSpline(u, t; extrapolate = true)
    @test A(-2.0) == P₁(-2.0) * ones(4)
    @test A(-0.5) == P₁(-0.5) * ones(4)
    @test A(0.7) == P₂(0.7) * ones(4)
    @test A(2.0) == P₂(2.0) * ones(4)

    u = [repeat(u[i], 1, 3) for i in 1:3]
    A = QuadraticSpline(u, t; extrapolate = true)
    @test A(-2.0) == P₁(-2.0) * ones(4, 3)
    @test A(-0.5) == P₁(-0.5) * ones(4, 3)
    @test A(0.7) == P₂(0.7) * ones(4, 3)
    @test A(2.0) == P₂(2.0) * ones(4, 3)

    # Test extrapolation
    u = [0.0, 1.0, 3.0]
    t = [-1.0, 0.0, 1.0]
    A = QuadraticSpline(u, t; extrapolate = true)
    @test A(-2.0) == 1.0
    @test A(2.0) == 5.0
    A = QuadraticSpline(u, t)
    @test_throws DataInterpolations.ExtrapolationError A(-2.0)
    @test_throws DataInterpolations.ExtrapolationError A(2.0)
end

@testset "CubicSpline Interpolation" begin
    u = [0.0, 1.0, 3.0]
    t = [-1.0, 0.0, 1.0]

    A = CubicSpline(u, t; extrapolate = true)

    # Solution
    P₁ = x -> 1 + 1.5x + x^2 + 0.5x^3 # for x ∈ [-1.0, 0.0]
    P₂ = x -> 1 + 1.5x + x^2 - 0.5x^3 # for x ∈ [0.0, 1.0]

    for (_t, _u) in zip(t, u)
        @test A(_t) == _u
    end
    for x in (-1.5, -0.5, -0.7)
        @test A(x) ≈ P₁(x)
    end
    for x in (0.3, 0.5, 1.5)
        @test A(x) ≈ P₂(x)
    end

    u_ = [0.0, 1.0, 3.0]' .* ones(4)
    u = [u_[:, i] for i in 1:size(u_, 2)]
    A = CubicSpline(u, t; extrapolate = true)
    for x in (-1.5, -0.5, -0.7)
        @test A(x) ≈ P₁(x) * ones(4)
    end
    for x in (0.3, 0.5, 1.5)
        @test A(x) ≈ P₂(x) * ones(4)
    end

    u = [repeat(u[i], 1, 3) for i in 1:3]
    A = CubicSpline(u, t; extrapolate = true)
    for x in (-1.5, -0.5, -0.7)
        @test A(x) ≈ P₁(x) * ones(4, 3)
    end
    for x in (0.3, 0.5, 1.5)
        @test A(x) ≈ P₂(x) * ones(4, 3)
    end

    # Test extrapolation
    u = [0.0, 1.0, 3.0]
    t = [-1.0, 0.0, 1.0]
    A = CubicSpline(u, t; extrapolate = true)
    @test A(-2.0) ≈ -2.0
    @test A(2.0) ≈ 4.0
    A = CubicSpline(u, t)
    @test_throws DataInterpolations.ExtrapolationError A(-2.0)
    @test_throws DataInterpolations.ExtrapolationError A(2.0)
end

@testset "BSplines" begin
    # BSpline Interpolation and Approximation
    t = [0, 62.25, 109.66, 162.66, 205.8, 252.3]
    u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]

    @testset "BSplineInterpolation" begin
        A = BSplineInterpolation(u, t, 2, :Uniform, :Uniform)

        @test [A(25.0), A(80.0)] == [13.454197730061425, 10.305633616059845]
        @test [A(190.0), A(225.0)] == [14.07428439395079, 11.057784141519251]
        @test [A(t[1]), A(t[end])] == [u[1], u[end]]

        # Test extrapolation
        A = BSplineInterpolation(u, t, 2, :Uniform, :Uniform; extrapolate = true)
        @test A(-1.0) == u[1]
        @test A(300.0) == u[end]
        A = BSplineInterpolation(u, t, 2, :Uniform, :Uniform)
        @test_throws DataInterpolations.ExtrapolationError A(-1.0)
        @test_throws DataInterpolations.ExtrapolationError A(300.0)

        A = BSplineInterpolation(u, t, 2, :ArcLen, :Average)

        @test [A(25.0), A(80.0)] == [13.363814458968486, 10.685201117692609]
        @test [A(190.0), A(225.0)] == [13.437481084762863, 11.367034741256463]
        @test [A(t[1]), A(t[end])] == [u[1], u[end]]

        # Test extrapolation
        A = BSplineInterpolation(u, t, 2, :ArcLen, :Average; extrapolate = true)
        @test A(-1.0) == u[1]
        @test A(300.0) == u[end]
        A = BSplineInterpolation(u, t, 2, :ArcLen, :Average)
        @test_throws DataInterpolations.ExtrapolationError A(-1.0)
        @test_throws DataInterpolations.ExtrapolationError A(300.0)
    end

    @testset "BSplineApprox" begin
        A = BSplineApprox(u, t, 2, 4, :Uniform, :Uniform)

        @test [A(25.0), A(80.0)] ≈ [12.979802931218234, 10.914310609953178]
        @test [A(190.0), A(225.0)] ≈ [13.851245975109263, 12.963685868886575]
        @test [A(t[1]), A(t[end])] ≈ [u[1], u[end]]

        # Test extrapolation
        A = BSplineApprox(u, t, 2, 4, :Uniform, :Uniform; extrapolate = true)
        @test A(-1.0) == u[1]
        @test A(300.0) == u[end]
        A = BSplineApprox(u, t, 2, 4, :Uniform, :Uniform)
        @test_throws DataInterpolations.ExtrapolationError A(-1.0)
        @test_throws DataInterpolations.ExtrapolationError A(300.0)
    end
end

# Curvefit Interpolation
rng = StableRNG(12345)
model(x, p) = @. p[1] / (1 + exp(x - p[2]))
t = range(-10, stop = 10, length = 40)
u = model(t, [1.0, 2.0]) + 0.01 * randn(rng, length(t))
p0 = [0.5, 0.5]

A = Curvefit(u, t, model, p0, LBFGS())

ts = [-7.0, -2.0, 0.0, 2.5, 5.0]
vs = [
    1.0013468217936277,
    0.9836755196317837,
    0.8833959853995836,
    0.3810348276782708,
    0.048062978598861855,
]
us = A.(ts)

@test vs ≈ us

# missing values handling tests
u = [1.0, 4.0, 9.0, 16.0, 25.0, missing, missing]
t = [1.0, 2.0, 3.0, 4.0, missing, 6.0, missing]
A = QuadraticInterpolation(u, t)

@test A(2.0) == 4.0
@test A(1.5) == 2.25
@test A(3.5) == 12.25
@test A(2.5) == 6.25

u = copy(hcat(u, u)')
A = QuadraticInterpolation(u, t)

@test A(2.0) == [4.0, 4.0]
@test A(1.5) == [2.25, 2.25]
@test A(3.5) == [12.25, 12.25]
@test A(2.5) == [6.25, 6.25]

# ForwardDiff compatibility with respect to coefficients

function square(INTERPOLATION_TYPE, c)  # elaborate way to write f(x) = x²
    xs = -4.0:2.0:4.0
    ys = [c^2 + x for x in xs]
    itp = INTERPOLATION_TYPE(ys, xs)
    return itp(0.0)
end

# generate versions of this function with different interpolators
f_quadratic_spline = c -> square(QuadraticSpline, c)
f_cubic_spline = c -> square(CubicSpline, c)

@test ForwardDiff.derivative(f_quadratic_spline, 2.0) ≈ 4.0
@test ForwardDiff.derivative(f_quadratic_spline, 4.0) ≈ 8.0
@test ForwardDiff.derivative(f_cubic_spline, 2.0) ≈ 4.0
@test ForwardDiff.derivative(f_cubic_spline, 4.0) ≈ 8.0
