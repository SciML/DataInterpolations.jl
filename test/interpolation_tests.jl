using DataInterpolations
using FindFirstFunctions: searchsortedfirstcorrelated
using StableRNGs
using Optim, ForwardDiff
using BenchmarkTools
using Unitful

function test_interpolation_type(T)
    @test T <: DataInterpolations.AbstractInterpolation
    @test hasfield(T, :u)
    @test hasfield(T, :t)
    @test hasfield(T, :extrapolation_right)
    @test hasfield(T, :extrapolation_left)
    @test hasfield(T, :iguesser)
    @test !isempty(methods(DataInterpolations._interpolate, (T, Any, Number)))
    @test !isempty(methods(DataInterpolations._integral, (T, Number, Number, Number)))
    @test !isempty(methods(DataInterpolations._derivative, (T, Any, Number)))
end

function test_cached_index(A)
    for t in range(first(A.t), last(A.t); length = 2 * length(A.t) - 1)
        A(t)
        idx = searchsortedfirstcorrelated(A.t, t, A.iguesser)
        @test abs(A.iguesser.idx_prev[] -
                  searchsortedfirstcorrelated(A.t, t, A.iguesser)) <= 2
    end
end

@testset "Linear Interpolation" begin
    test_interpolation_type(LinearInterpolation)

    for t in (1.0:10.0, 1.0collect(1:10))
        u = 2.0collect(1:10)
        #t = 1.0collect(1:10)
        A = LinearInterpolation(u, t; extrapolation = ExtrapolationType.Extension)

        for (_t, _u) in zip(t, u)
            @test A(_t) == _u
        end
        @test A(0) == 0.0
        @test A(5.5) == 11.0
        @test A(11) == 22

        u = vcat(2.0collect(1:10)', 3.0collect(1:10)')
        A = LinearInterpolation(u, t; extrapolation = ExtrapolationType.Extension)

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
        A = LinearInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
        @test A(0) == [0.0, 0.0, 0.0]
        @test A(5.5) == [11.0, 16.5, 22.0]
        @test A(11) == [22.0, 33.0, 44.0]
    end

    x = 1:10
    y = 2:4
    u_ = x' .* y
    u = [u_[:, i:(i + 1)] for i in 1:2:10]
    t = 1.0collect(2:2:10)
    A = LinearInterpolation(u, t; extrapolation = ExtrapolationType.Extension)

    @test A(0) == [-2.0 0.0; -3.0 0.0; -4.0 0.0]
    @test A(3) == [4.0 6.0; 6.0 9.0; 8.0 12.0]
    @test A(5) == [8.0 10.0; 12.0 15.0; 16.0 20.0]
    test_cached_index(A)

    # with NaNs (#113)
    u = [NaN, 1.0, 2.0, 3.0]
    t = 1:4
    A = LinearInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
    @test isnan(A(1.0))
    @test A(2.0) == 1.0
    @test A(2.5) == 1.5
    @test A(3.0) == 2.0
    @test A(4.0) == 3.0

    u = [0.0, NaN, 2.0, 3.0]
    A = LinearInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
    @test A(1.0) == 0.0
    @test isnan(A(2.0))
    @test isnan(A(2.5))
    @test A(3.0) == 2.0
    @test A(4.0) == 3.0

    u = [0.0, 1.0, NaN, 3.0]
    A = LinearInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
    @test A(1.0) == 0.0
    @test A(2.0) == 1.0
    @test isnan(A(2.5))
    @test isnan(A(3.0))
    @test A(4.0) == 3.0

    u = [0.0, 1.0, 2.0, NaN]
    A = LinearInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
    @test A(1.0) == 0.0
    @test A(2.0) == 1.0
    @test A(3.0) == 2.0
    @test isnan(A(3.5))
    @test isnan(A(4.0))

    u = [0.0, 1.0, 2.0, NaN]
    A = LinearInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
    @test A(1.0) == 0.0
    @test A(2.0) == 1.0
    @test A(3.0) == 2.0
    @test isnan(A(3.5))
    @test isnan(A(4.0))

    # Test type stability
    u = Float32.(1:5)
    t = Float32.(1:5)
    A1 = LinearInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
    u = 1:5
    t = 1:5
    A2 = LinearInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
    u = [1 // i for i in 1:5]
    t = (1:5)
    A3 = LinearInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
    u = [1 // i for i in 1:5]
    t = [1 // (6 - i) for i in 1:5]
    A4 = LinearInterpolation(u, t; extrapolation = ExtrapolationType.Extension)

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

    # NaN time value for Unitful arrays: issue #365
    t = (0:3)u"s" # Unitful quantities  
    u = [0, -2, -1, -2]u"m"
    A = LinearInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
    @test isnan(A(NaN * u"s"))

    # Nan time value:
    t = 0.0:3  # Floats
    u = [0, -2, -1, -2]
    A = LinearInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
    dA = t -> ForwardDiff.derivative(A, t)
    @test isnan(dA(NaN))

    t = 0:3  # Integers
    u = [0, -2, -1, -2]
    A = LinearInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
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
    A = LinearInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
    @test A(0) == fill(0.0)
    @test A(5.5) == fill(11.0)
    @test A(11) == fill(22)

    # Test constant -Inf interpolation
    u = [-Inf, -Inf]
    t = [0.0, 1.0]
    A = LinearInterpolation(u, t)
    @test A(0.0) == -Inf
    @test A(0.5) == -Inf

    # Test extrapolation
    u = 2.0collect(1:10)
    t = 1.0collect(1:10)
    A = LinearInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
    @test A(-1.0) == -2.0
    @test A(11.0) == 22.0
    A = LinearInterpolation(u, t)
    @test_throws DataInterpolations.LeftExtrapolationError A(-1.0)
    @test_throws DataInterpolations.RightExtrapolationError A(11.0)
    @test_throws DataInterpolations.LeftExtrapolationError A([-1.0, 11.0])
end

@testset "Quadratic Interpolation" begin
    test_interpolation_type(QuadraticInterpolation)

    u = [1.0, 4.0, 9.0, 16.0]
    t = [1.0, 2.0, 3.0, 4.0]
    A = QuadraticInterpolation(u, t; extrapolation = ExtrapolationType.Extension)

    for (_t, _u) in zip(t, u)
        @test A(_t) == _u
    end
    @test A(0.0) == 0.0
    @test A(1.5) == 2.25
    @test A(2.5) == 6.25
    @test A(3.5) == 12.25
    @test A(5.0) == 25
    test_cached_index(A)

    # backward-looking interpolation
    u = [1.0, 4.0, 9.0, 16.0]
    t = [1.0, 2.0, 3.0, 4.0]
    A = QuadraticInterpolation(
        u, t, :Backward; extrapolation = ExtrapolationType.Extension)

    for (_t, _u) in zip(t, u)
        @test A(_t) == _u
    end
    @test A(0.0) == 0.0
    @test A(1.5) == 2.25
    @test A(2.5) == 6.25
    @test A(3.5) == 12.25
    @test A(5.0) == 25
    test_cached_index(A)

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

    test_cached_index(A_f)
    test_cached_index(A_b)

    # Matrix interpolation test
    u = [1.0 4.0 9.0 16.0; 1.0 4.0 9.0 16.0]
    A = QuadraticInterpolation(u, t; extrapolation = ExtrapolationType.Extension)

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
    A = QuadraticInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
    @test A(0) == zeros(5)
    @test A(1.5) == 2.25 * ones(5)
    @test A(2.5) == 6.25 * ones(5)
    @test A(3.5) == 12.25 * ones(5)
    @test A(5.0) == 25.0 * ones(5)

    u = [repeat(u[i], 1, 3) for i in 1:4]
    A = QuadraticInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
    @test A(0) == zeros(5, 3)
    @test A(1.5) == 2.25 * ones(5, 3)
    @test A(2.5) == 6.25 * ones(5, 3)
    @test A(3.5) == 12.25 * ones(5, 3)
    @test A(5.0) == 25.0 * ones(5, 3)

    # Test extrapolation
    u = [1.0, 4.5, 6.0, 2.0]
    t = [1.0, 2.0, 3.0, 4.0]
    A = QuadraticInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
    @test A(0.0) == -4.5
    @test A(5.0) == -7.5
    A = QuadraticInterpolation(u, t)
    @test_throws DataInterpolations.LeftExtrapolationError A(0.0)
    @test_throws DataInterpolations.RightExtrapolationError A(5.0)
end

@testset "Lagrange Interpolation" begin
    test_interpolation_type(LagrangeInterpolation)

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
    A = LagrangeInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
    @test A(0.0) == 0.0
    @test A(4.0) == 16.0
    A = LagrangeInterpolation(u, t)
    @test_throws DataInterpolations.LeftExtrapolationError A(-1.0)
    @test_throws DataInterpolations.RightExtrapolationError A(4.0)
end

@testset "Akima Interpolation" begin
    test_interpolation_type(AkimaInterpolation)

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
    test_cached_index(A)

    # Test extrapolation
    A = AkimaInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
    @test A(-1.0) ≈ -5.0
    @test A(11.0) ≈ -3.924742268041234
    A = AkimaInterpolation(u, t)
    @test_throws DataInterpolations.LeftExtrapolationError A(-1.0)
    @test_throws DataInterpolations.RightExtrapolationError A(11.0)
end

@testset "ConstantInterpolation" begin
    test_interpolation_type(ConstantInterpolation)

    t = [1.0, 2.0, 3.0, 4.0]

    @testset "Vector case" for u in [[1.0, 2.0, 0.0, 1.0], ["B", "C", "A", "B"]]
        A = ConstantInterpolation(
            u, t, dir = :right; extrapolation = ExtrapolationType.Extension)
        @test A(0.5) == u[1]
        @test A(1.0) == u[1]
        @test A(1.5) == u[2]
        @test A(2.0) == u[2]
        @test A(2.5) == u[3]
        @test A(3.0) == u[3]
        @test A(3.5) == u[1]
        @test A(4.0) == u[1]
        @test A(4.5) == u[1]
        test_cached_index(A)

        A = ConstantInterpolation(u, t; extrapolation = ExtrapolationType.Extension) # dir=:left is default
        @test A(0.5) == u[1]
        @test A(1.0) == u[1]
        @test A(1.5) == u[1]
        @test A(2.0) == u[2]
        @test A(2.5) == u[2]
        @test A(3.0) == u[3]
        @test A(3.5) == u[3]
        @test A(4.0) == u[1]
        @test A(4.5) == u[1]
        test_cached_index(A)
    end

    @testset "Matrix case" for u in [
        [1.0 2.0 0.0 1.0; 1.0 2.0 0.0 1.0],
        ["B" "C" "A" "B"; "B" "C" "A" "B"]
    ]
        A = ConstantInterpolation(
            u, t, dir = :right; extrapolation = ExtrapolationType.Extension)
        @test A(0.5) == u[:, 1]
        @test A(1.0) == u[:, 1]
        @test A(1.5) == u[:, 2]
        @test A(2.0) == u[:, 2]
        @test A(2.5) == u[:, 3]
        @test A(3.0) == u[:, 3]
        @test A(3.5) == u[:, 1]
        @test A(4.0) == u[:, 1]
        @test A(4.5) == u[:, 1]
        test_cached_index(A)

        A = ConstantInterpolation(u, t; extrapolation = ExtrapolationType.Extension) # dir=:left is default
        @test A(0.5) == u[:, 1]
        @test A(1.0) == u[:, 1]
        @test A(1.5) == u[:, 1]
        @test A(2.0) == u[:, 2]
        @test A(2.5) == u[:, 2]
        @test A(3.0) == u[:, 3]
        @test A(3.5) == u[:, 3]
        @test A(4.0) == u[:, 1]
        @test A(4.5) == u[:, 1]
        test_cached_index(A)
    end

    @testset "Vector of Vectors case" for u in [
        [[1.0, 2.0], [0.0, 1.0], [1.0, 2.0], [0.0, 1.0]],
        [["B", "C"], ["A", "B"], ["B", "C"], ["A", "B"]]]
        A = ConstantInterpolation(
            u, t, dir = :right; extrapolation = ExtrapolationType.Extension)
        @test A(0.5) == u[1]
        @test A(1.0) == u[1]
        @test A(1.5) == u[2]
        @test A(2.0) == u[2]
        @test A(2.5) == u[3]
        @test A(3.0) == u[3]
        @test A(3.5) == u[4]
        @test A(4.0) == u[4]
        @test A(4.5) == u[4]
        test_cached_index(A)

        A = ConstantInterpolation(u, t; extrapolation = ExtrapolationType.Extension) # dir=:left is default
        @test A(0.5) == u[1]
        @test A(1.0) == u[1]
        @test A(1.5) == u[1]
        @test A(2.0) == u[2]
        @test A(2.5) == u[2]
        @test A(3.0) == u[3]
        @test A(3.5) == u[3]
        @test A(4.0) == u[4]
        @test A(4.5) == u[4]
        test_cached_index(A)
    end

    @testset "Vector of Matrices case" for u in [
        [[1.0 2.0; 1.0 2.0], [0.0 1.0; 0.0 1.0], [1.0 2.0; 1.0 2.0], [0.0 1.0; 0.0 1.0]],
        [["B" "C"; "B" "C"], ["A" "B"; "A" "B"], ["B" "C"; "B" "C"], ["A" "B"; "A" "B"]]]
        A = ConstantInterpolation(
            u, t, dir = :right; extrapolation = ExtrapolationType.Extension)
        @test A(0.5) == u[1]
        @test A(1.0) == u[1]
        @test A(1.5) == u[2]
        @test A(2.0) == u[2]
        @test A(2.5) == u[3]
        @test A(3.0) == u[3]
        @test A(3.5) == u[4]
        @test A(4.0) == u[4]
        @test A(4.5) == u[4]
        test_cached_index(A)

        A = ConstantInterpolation(u, t; extrapolation = ExtrapolationType.Extension) # dir=:left is default
        @test A(0.5) == u[1]
        @test A(1.0) == u[1]
        @test A(1.5) == u[1]
        @test A(2.0) == u[2]
        @test A(2.5) == u[2]
        @test A(3.0) == u[3]
        @test A(3.5) == u[3]
        @test A(4.0) == u[4]
        @test A(4.5) == u[4]
        test_cached_index(A)
    end

    # Test extrapolation
    u = [1.0, 2.0, 0.0, 1.0]
    A = ConstantInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
    @test A(-1.0) == 1.0
    @test A(11.0) == 1.0
    A = ConstantInterpolation(u, t)
    @test_throws DataInterpolations.LeftExtrapolationError A(-1.0)
    @test_throws DataInterpolations.RightExtrapolationError A(11.0)

    # Test extrapolation with infs with regularly spaced t
    u = [1.67e7, 1.6867e7, 1.7034e7, 1.7201e7, 1.7368e7]
    t = [0.0, 0.1, 0.2, 0.3, 0.4]
    A = ConstantInterpolation(u, t; extrapolation = ExtrapolationType.Extension)
    @test A(Inf) == last(u)
    @test A(-Inf) == first(u)

    # Test extrapolation of integer output
    itp = ConstantInterpolation([2], [0.0]; extrapolation = ExtrapolationType.Constant)
    @test itp(1.0) === 2
    @test itp(-1.0) === 2
end

@testset "QuadraticSpline Interpolation" begin
    test_interpolation_type(QuadraticSpline)

    u = [0.0, 1.0, 3.0]
    t = [-1.0, 0.0, 1.0]

    A = QuadraticSpline(u, t; extrapolation = ExtrapolationType.Extension)

    # Solution
    P₁ = x -> 0.5 * (x + 1) * (x + 2)

    for (_t, _u) in zip(t, u)
        @test A(_t) == _u
    end
    @test A(-2.0) == P₁(-2.0)
    @test A(-0.5) == P₁(-0.5)
    @test A(0.7) == P₁(0.7)
    @test A(2.0) == P₁(2.0)
    test_cached_index(A)

    u_ = [0.0, 1.0, 3.0]' .* ones(4)
    u = [u_[:, i] for i in 1:size(u_, 2)]
    A = QuadraticSpline(u, t; extrapolation = ExtrapolationType.Extension)
    @test A(-2.0) == P₁(-2.0) * ones(4)
    @test A(-0.5) == P₁(-0.5) * ones(4)
    @test A(0.7) == P₁(0.7) * ones(4)
    @test A(2.0) == P₁(2.0) * ones(4)

    u = [repeat(u[i], 1, 3) for i in 1:3]
    A = QuadraticSpline(u, t; extrapolation = ExtrapolationType.Extension)
    @test A(-2.0) == P₁(-2.0) * ones(4, 3)
    @test A(-0.5) == P₁(-0.5) * ones(4, 3)
    @test A(0.7) == P₁(0.7) * ones(4, 3)
    @test A(2.0) == P₁(2.0) * ones(4, 3)

    # Test extrapolation
    u = [0.0, 1.0, 3.0]
    t = [-1.0, 0.0, 1.0]
    A = QuadraticSpline(u, t; extrapolation = ExtrapolationType.Extension)
    @test A(-2.0) == 0.0
    @test A(2.0) == 6.0
    A = QuadraticSpline(u, t)
    @test_throws DataInterpolations.LeftExtrapolationError A(-2.0)
    @test_throws DataInterpolations.RightExtrapolationError A(2.0)
end

@testset "CubicSpline Interpolation" begin
    test_interpolation_type(CubicSpline)

    u = [0.0, 1.0, 3.0]
    t = [-1.0, 0.0, 1.0]

    A = CubicSpline(u, t; extrapolation = ExtrapolationType.Extension)
    test_cached_index(A)

    # Solution
    P₁ = x -> 1 + 1.5x + 0.75 * x^2 + 0.25 * x^3 # for x ∈ [-1.0, 0.0]
    P₂ = x -> 1 + 1.5x + 0.75 * x^2 - 0.25 * x^3 # for x ∈ [0.0, 1.0]

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
    A = CubicSpline(u, t; extrapolation = ExtrapolationType.Extension)
    for x in (-1.5, -0.5, -0.7)
        @test A(x) ≈ P₁(x) * ones(4)
    end
    for x in (0.3, 0.5, 1.5)
        @test A(x) ≈ P₂(x) * ones(4)
    end

    u = [repeat(u[i], 1, 3) for i in 1:3]
    A = CubicSpline(u, t; extrapolation = ExtrapolationType.Extension)
    for x in (-1.5, -0.5, -0.7)
        @test A(x) ≈ P₁(x) * ones(4, 3)
    end
    for x in (0.3, 0.5, 1.5)
        @test A(x) ≈ P₂(x) * ones(4, 3)
    end

    # Test extrapolation
    u = [0.0, 1.0, 3.0]
    t = [-1.0, 0.0, 1.0]
    A = CubicSpline(u, t; extrapolation = ExtrapolationType.Extension)
    @test A(-2.0) ≈ -1.0
    @test A(2.0) ≈ 5.0
    A = CubicSpline(u, t)
    @test_throws DataInterpolations.LeftExtrapolationError A(-2.0)
    @test_throws DataInterpolations.RightExtrapolationError A(2.0)

    @testset "AbstractMatrix" begin
        t = 0.1:0.1:1.0
        u = [sin.(t) cos.(t)]' |> collect
        c = CubicSpline(u, t)
        t_test = 0.1:0.05:1.0
        u_test = reduce(hcat, c.(t_test))
        @test isapprox(u_test[1, :], sin.(t_test), atol = 1e-3)
        @test isapprox(u_test[2, :], cos.(t_test), atol = 1e-3)
    end
    @testset "AbstractArray{T, 3}" begin
        f3d(t) = [sin(t) cos(t);
                  0.0 cos(2t)]
        t = 0.1:0.1:1.0
        u3d = f3d.(t)
        c = CubicSpline(u3d, t)
        t_test = 0.1:0.05:1.0
        u_test = reduce(hcat, c.(t_test))
        f_test = reduce(hcat, f3d.(t_test))
        @test isapprox(u_test, f_test, atol = 1e-2)
    end
end

@testset "BSplines" begin

    # BSpline Interpolation and Approximation
    @testset "BSplineInterpolation" begin
        t = [0, 62.25, 109.66, 162.66, 205.8, 252.3]
        u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
        test_interpolation_type(BSplineInterpolation)
        A = BSplineInterpolation(u, t, 2, :Uniform, :Uniform)

        @test [A(25.0), A(80.0)] == [13.454197730061425, 10.305633616059845]
        @test [A(190.0), A(225.0)] == [14.07428439395079, 11.057784141519251]
        @test [A(t[1]), A(t[end])] == [u[1], u[end]]
        test_cached_index(A)

        # Test extrapolation
        A = BSplineInterpolation(
            u, t, 2, :Uniform, :Uniform; extrapolation = ExtrapolationType.Extension)
        @test A(-1.0) == u[1]
        @test A(300.0) == u[end]
        A = BSplineInterpolation(u, t, 2, :Uniform, :Uniform)
        @test_throws DataInterpolations.LeftExtrapolationError A(-1.0)
        @test_throws DataInterpolations.RightExtrapolationError A(300.0)

        A = BSplineInterpolation(u, t, 2, :ArcLen, :Average)

        @test [A(25.0), A(80.0)] == [13.363814458968486, 10.685201117692609]
        @test [A(190.0), A(225.0)] == [13.437481084762863, 11.367034741256463]
        @test [A(t[1]), A(t[end])] == [u[1], u[end]]

        @test_throws ErrorException("BSplineInterpolation needs at least d + 1, i.e. 4 points.") BSplineInterpolation(
            u[1:3], t[1:3], 3, :Uniform, :Uniform)
        @test_throws ErrorException("BSplineInterpolation needs at least d + 1, i.e. 5 points.") BSplineInterpolation(
            u[1:4], t[1:4], 4, :ArcLen, :Average)
        @test_nowarn BSplineInterpolation(u[1:3], t[1:3], 2, :Uniform, :Uniform)

        # Test extrapolation
        A = BSplineInterpolation(
            u, t, 2, :ArcLen, :Average; extrapolation = ExtrapolationType.Extension)
        @test A(-1.0) == u[1]
        @test A(300.0) == u[end]
        A = BSplineInterpolation(u, t, 2, :ArcLen, :Average)
        @test_throws DataInterpolations.LeftExtrapolationError A(-1.0)
        @test_throws DataInterpolations.RightExtrapolationError A(300.0)

        @testset "AbstractMatrix" begin
            t = 0.1:0.1:1.0
            u2d = [sin.(t) cos.(t)]' |> collect
            A = BSplineInterpolation(u2d, t, 2, :Uniform, :Uniform)
            t_test = 0.1:0.05:1.0
            u_test = reduce(hcat, A.(t_test))
            @test isapprox(u_test[1, :], sin.(t_test), atol = 1e-3)
            @test isapprox(u_test[2, :], cos.(t_test), atol = 1e-3)

            A = BSplineInterpolation(u2d, t, 2, :ArcLen, :Average)
            u_test = reduce(hcat, A.(t_test))
            @test isapprox(u_test[1, :], sin.(t_test), atol = 1e-3)
            @test isapprox(u_test[2, :], cos.(t_test), atol = 1e-3)
        end
        @testset "AbstractArray{T, 3}" begin
            f3d(t) = [sin(t) cos(t);
                      0.0 cos(2t)]
            t = 0.1:0.1:1.0
            u3d = cat(f3d.(t)..., dims = 3)
            A = BSplineInterpolation(u3d, t, 2, :Uniform, :Uniform)
            t_test = 0.1:0.05:1.0
            u_test = reduce(hcat, A.(t_test))
            f_test = reduce(hcat, f3d.(t_test))
            @test isapprox(u_test, f_test, atol = 1e-2)

            A = BSplineInterpolation(u3d, t, 2, :ArcLen, :Average)
            t_test = 0.1:0.05:1.0
            u_test = reduce(hcat, A.(t_test))
            @test isapprox(u_test, f_test, atol = 1e-2)
        end
    end

    @testset "BSplineApprox" begin
        test_interpolation_type(BSplineApprox)
        t = [0, 62.25, 109.66, 162.66, 205.8, 252.3]
        u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
        A = BSplineApprox(u, t, 2, 4, :Uniform, :Uniform)

        @test [A(25.0), A(80.0)] ≈ [12.979802931218234, 10.914310609953178]
        @test [A(190.0), A(225.0)] ≈ [13.851245975109263, 12.963685868886575]
        @test [A(t[1]), A(t[end])] ≈ [u[1], u[end]]
        test_cached_index(A)

        @test_throws ErrorException("BSplineApprox needs at least d + 1, i.e. 3 control points.") BSplineApprox(
            u, t, 2, 2, :Uniform, :Uniform)
        @test_throws ErrorException("BSplineApprox needs at least d + 1, i.e. 4 control points.") BSplineApprox(
            u, t, 3, 3, :ArcLen, :Average)
        @test_nowarn BSplineApprox(u, t, 2, 3, :Uniform, :Uniform)

        # Test extrapolation
        A = BSplineApprox(
            u, t, 2, 4, :Uniform, :Uniform; extrapolation = ExtrapolationType.Extension)
        @test A(-1.0) == u[1]
        @test A(300.0) == u[end]
        A = BSplineApprox(u, t, 2, 4, :Uniform, :Uniform)
        @test_throws DataInterpolations.LeftExtrapolationError A(-1.0)
        @test_throws DataInterpolations.RightExtrapolationError A(300.0)

        @testset "AbstractMatrix" begin
            t = 0.1:0.1:1.0
            u2d = [sin.(t) cos.(t)]' |> collect
            A = BSplineApprox(u2d, t, 2, 5, :Uniform, :Uniform)
            t_test = 0.1:0.05:1.0
            u_test = reduce(hcat, A.(t_test))
            @test isapprox(u_test[1, :], sin.(t_test), atol = 1e-3)
            @test isapprox(u_test[2, :], cos.(t_test), atol = 1e-3)

            A = BSplineApprox(u2d, t, 2, 5, :ArcLen, :Average)
            u_test = reduce(hcat, A.(t_test))
            @test isapprox(u_test[1, :], sin.(t_test), atol = 1e-2)
            @test isapprox(u_test[2, :], cos.(t_test), atol = 1e-2)
        end
        @testset "AbstractArray{T, 3}" begin
            f3d(t) = [sin(t) cos(t);
                      0.0 cos(2t)]
            t = 0.1:0.1:1.0
            u3d = cat(f3d.(t)..., dims = 3)
            A = BSplineApprox(u3d, t, 2, 6, :Uniform, :Uniform)
            t_test = 0.1:0.05:1.0
            u_test = reduce(hcat, A.(t_test))
            f_test = reduce(hcat, f3d.(t_test))
            @test isapprox(u_test, f_test, atol = 1e-2)

            A = BSplineApprox(u3d, t, 2, 7, :ArcLen, :Average)
            t_test = 0.1:0.05:1.0
            u_test = reduce(hcat, A.(t_test))
            @test isapprox(u_test, f_test, atol = 1e-2)
        end
    end
end

@testset "Cubic Hermite Spline" begin
    test_interpolation_type(CubicHermiteSpline)

    du = [-0.047, -0.058, 0.054, 0.012, -0.068, 0.0]
    u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
    t = [0.0, 62.25, 109.66, 162.66, 205.8, 252.3]
    A = CubicHermiteSpline(du, u, t; extrapolation = ExtrapolationType.Extension)
    @test A.(t) ≈ u
    @test A(100.0)≈10.106770 rtol=1e-5
    @test A(300.0)≈9.901542 rtol=1e-5
    test_cached_index(A)
    push!(u, 1.0)
    @test_throws AssertionError CubicHermiteSpline(du, u, t)
end

@testset "PCHIPInterpolation" begin
    u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
    t = [0.0, 62.25, 109.66, 162.66, 205.8, 250.0]
    A = PCHIPInterpolation(u, t)
    @test A isa CubicHermiteSpline
    ts = 0.0:0.1:250.0
    us = A(ts)
    @test all(minimum(u) .<= us)
    @test all(maximum(u) .>= us)
    @test all(A.du[3:4] .== 0.0)
end

@testset "Quintic Hermite Spline" begin
    test_interpolation_type(QuinticHermiteSpline)

    ddu = [0.0, -0.00033, 0.0051, -0.0067, 0.0029, 0.0]
    du = [-0.047, -0.058, 0.054, 0.012, -0.068, 0.0]
    u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
    t = [0.0, 62.25, 109.66, 162.66, 205.8, 252.3]
    A = QuinticHermiteSpline(ddu, du, u, t; extrapolation = ExtrapolationType.Extension)
    @test A.(t) ≈ u
    @test A(100.0)≈10.107996 rtol=1e-5
    @test A(300.0)≈11.364162 rtol=1e-5
    test_cached_index(A)
    push!(u, 1.0)
    @test_throws AssertionError QuinticHermiteSpline(ddu, du, u, t)
end

@testset "Curvefit" begin
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
        0.048062978598861855
    ]
    us = A.(ts)
    @test vs ≈ us

    # Test extrapolation
    A = Curvefit(u, t, model, p0, LBFGS(); extrapolate = true)
    @test A(15.0) == model(15.0, A.pmin)
    A = Curvefit(u, t, model, p0, LBFGS())
    @test_throws DataInterpolations.ExtrapolationError A(15.0)
end

@testset "Type of vector returned" begin
    # Issue https://github.com/SciML/DataInterpolations.jl/issues/253
    t1 = Float32[0.1, 0.2, 0.3, 0.4, 0.5]
    t2 = Float64[0.1, 0.2, 0.3, 0.4, 0.5]
    interps_and_types = [
        (LinearInterpolation(t1, t1), Float32),
        (LinearInterpolation(t1, t2), Float32),
        (LinearInterpolation(t2, t1), Float64),
        (LinearInterpolation(t2, t2), Float64)
    ]
    for i in eachindex(interps_and_types)
        @test eltype(interps_and_types[i][1](t1)) == interps_and_types[i][2]
    end
end

@testset "Plugging vector timepoints" begin
    # Issue https://github.com/SciML/DataInterpolations.jl/issues/267
    t = Float64[1.0, 2.0, 3.0, 4.0, 5.0]
    @testset "utype - Vectors" begin
        interp = LinearInterpolation(rand(5), t)
        @test interp(t) isa Vector{Float64}
    end
    @testset "utype - Vector of Vectors" begin
        interp = LinearInterpolation([rand(2) for _ in 1:5], t)
        @test interp(t) isa Vector{Vector{Float64}}
    end
    @testset "utype - Matrix" begin
        interp = LinearInterpolation(rand(2, 5), t)
        @test interp(t) isa Matrix{Float64}
    end
end

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
