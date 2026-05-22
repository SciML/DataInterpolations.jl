using DataInterpolations
using Reactant
using Test

function test_reactant_interpolation(interp, t_vals; atol = 1.0e-12)
    for t_val in t_vals
        expected = interp(t_val)
        f = Reactant.compile(
            (t) -> interp(t),
            (Reactant.to_rarray(t_val; track_numbers = true),)
        )
        result = Float64(f(Reactant.to_rarray(t_val; track_numbers = true)))
        @test isapprox(result, expected; atol = atol)
    end
    return
end

@testset "LinearInterpolation" begin
    @testset "basic" begin
        interp = LinearInterpolation([0.0, 1.0, 4.0], [0.0, 0.5, 1.0])
        test_reactant_interpolation(interp, [0.0, 0.25, 0.5, 0.75, 1.0])
    end

    @testset "Extension extrapolation" begin
        interp = LinearInterpolation(
            [0.0, 2.0], [0.0, 1.0]; extrapolation = ExtrapolationType.Extension
        )
        test_reactant_interpolation(interp, [-0.5, 0.0, 0.5, 1.0, 1.5])
    end

    @testset "Constant extrapolation" begin
        interp = LinearInterpolation(
            [0.0, 2.0], [0.0, 1.0];
            extrapolation_left = ExtrapolationType.Constant,
            extrapolation_right = ExtrapolationType.Constant
        )
        test_reactant_interpolation(interp, [-0.5, 0.0, 0.5, 1.0, 1.5])
    end

    @testset "many knots" begin
        u = sin.(range(0, 2π; length = 20))
        t = collect(range(0.0, 1.0; length = 20))
        interp = LinearInterpolation(u, t)
        test_reactant_interpolation(interp, range(0.0, 1.0; length = 50))
    end
end

@testset "ConstantInterpolation" begin
    @testset "dir=:left" begin
        interp = ConstantInterpolation([1.0, 3.0, 5.0], [0.0, 0.5, 1.0])
        test_reactant_interpolation(interp, [0.0, 0.25, 0.5, 0.75, 1.0])
    end

    @testset "dir=:right" begin
        interp = ConstantInterpolation(
            [1.0, 3.0, 5.0], [0.0, 0.5, 1.0]; dir = :right
        )
        test_reactant_interpolation(interp, [0.0, 0.25, 0.5, 0.75, 1.0])
    end
end

@testset "QuadraticInterpolation" begin
    interp = QuadraticInterpolation([0.0, 1.0, 0.0], [0.0, 0.5, 1.0])
    test_reactant_interpolation(interp, [0.0, 0.25, 0.5, 0.75, 1.0]; atol = 1.0e-10)
end
