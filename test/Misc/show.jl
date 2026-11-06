using DataInterpolations
using CurveFit, StableRNGs

t = [1.0, 2.0, 3.0, 4.0, 5.0]
x = [1.0, 2.0, 3.0, 4.0, 5.0]

@testset "Generic Cases" begin
    function test_show_line(A)
        @testset "$(nameof(typeof(A)))" begin
            @test startswith(
                sprint(io -> show(io, MIME"text/plain"(), A)),
                "$(nameof(typeof(A))) with $(length(A.t)) points\n"
            )
        end
    end
    methods = [
        LinearInterpolation(x, t),
        AkimaInterpolation(x, t),
        QuadraticSpline(x, t),
        CubicSpline(x, t),
    ]
    test_show_line.(methods)
end

@testset "Specific Cases" begin
    @testset "QuadraticInterpolation" begin
        A = QuadraticInterpolation(x, t)
        @test startswith(
            sprint(io -> show(io, MIME"text/plain"(), A)),
            "QuadraticInterpolation with 5 points, Forward mode\n"
        )
    end
    @testset "LagrangeInterpolation" begin
        A = LagrangeInterpolation(x, t)
        @test startswith(
            sprint(io -> show(io, MIME"text/plain"(), A)),
            "LagrangeInterpolation with 5 points, with order 4\n"
        )
    end
    @testset "ConstantInterpolation" begin
        A = ConstantInterpolation(x, t)
        @test startswith(
            sprint(io -> show(io, MIME"text/plain"(), A)),
            "ConstantInterpolation with 5 points, in left direction\n"
        )
    end
    @testset "BSplineInterpolation" begin
        A = BSplineInterpolation(x, t, 3, :Uniform)
        @test startswith(
            sprint(io -> show(io, MIME"text/plain"(), A)),
            "BSplineInterpolation with 5 points, with degree 3\n"
        )
    end
    @testset "BSplineApprox" begin
        A = BSplineApprox(x, t, 2, 4, :Uniform)
        @test startswith(
            sprint(io -> show(io, MIME"text/plain"(), A)),
            "BSplineApprox with 5 points, with degree 2, number of control points 4\n"
        )
    end
end

@testset "repr round-trip" begin
    # `repr` names FindFirstFunctions, whose search state every interpolation stores.
    m = Module()
    Core.eval(m, :(using DataInterpolations))
    Core.eval(m, :(const FindFirstFunctions = DataInterpolations.FindFirstFunctions))
    function test_repr_roundtrip(A)
        @testset "$(nameof(typeof(A)))" begin
            B = Core.eval(m, Meta.parse(repr(A)))
            @test typeof(B) == typeof(A)
            @test repr(B) == repr(A)
            tm = (first(A.t) + last(A.t)) / 2
            @test B(tm) == A(tm)
        end
    end
    du = ones(5)
    ddu = zeros(5)
    interpolations = [
        LinearInterpolation(x, t),
        LinearInterpolation(x, t; cache_parameters = true),
        QuadraticInterpolation(x, t),
        LagrangeInterpolation(x, t),
        AkimaInterpolation(x, t),
        ConstantInterpolation(x, t),
        SmoothedConstantInterpolation(x, t),
        QuadraticSpline(x, t),
        CubicSpline(x, t),
        BSplineInterpolation(x, t, 3, :Uniform),
        BSplineApprox(x, t, 2, 4, :Uniform),
        CubicHermiteSpline(du, x, t),
        QuinticHermiteSpline(ddu, du, x, t),
        SmoothArcLengthInterpolation([0.3 -0.2 0.0; -1.5 0.2 0.0; 3.1 -1.5 0.0]; m = 5),
        DataInterpolations.invert_integral(LinearInterpolation(x, t)),
        DataInterpolations.invert_integral(ConstantInterpolation(x, t)),
    ]
    test_repr_roundtrip.(interpolations)
end

@testset "Higher-rank u (AbstractArray{T, 3})" begin
    f3d(t_) = [sin(t_) cos(t_); 0.0 cos(2t_)]
    u3d = cat(f3d.(t)...; dims = 3)
    A = LinearInterpolation(u3d, t)
    str = sprint(io -> show(io, MIME"text/plain"(), A))
    @test startswith(str, "LinearInterpolation with 5 points\n")
    @test all(occursin("u$i", str) for i in 1:4)
end

@testset "CurveFit" begin
    rng = StableRNG(12345)
    model(x, p) = @. p[1] / (1 + exp(x - p[2]))
    t = range(-10, stop = 10, length = 40)
    u = model(t, [1.0, 2.0]) + 0.01 * randn(rng, length(t))
    p0 = [0.5, 0.5]
    A = Curvefit(u, t, model, p0)
    @test startswith(
        sprint(io -> show(io, MIME"text/plain"(), A)),
        "Curvefit with 40 points.\n"
    )
    # The model function does not print as a call, so rebuild from the fields directly.
    B = typeof(A)((getfield(A, f) for f in fieldnames(typeof(A)))...)
    @test B(0.5) == A(0.5)
end
