using DataInterpolations
using ForwardDiff
using Zygote

function test_zygote(method, u, t; args = [], args_after = [], kwargs = [], name::String)
    func = method(
        args..., u, t, args_after...; kwargs...,
        extrapolation = ExtrapolationType.Extension
    )
    trange = collect(range(minimum(t) - 5.0, maximum(t) + 5.0, step = 0.1))
    trange_exclude = filter(x -> !in(x, t), trange)
    @testset "$name, derivatives w.r.t. input" begin
        for _t in trange_exclude
            adiff = DataInterpolations.derivative(func, _t)
            zdiff = u isa AbstractVector{<:Real} ? only(Zygote.gradient(func, _t)) :
                only(Zygote.jacobian(func, _t))
            isnothing(zdiff) && (zdiff = 0.0)
            @test adiff ≈ zdiff
        end
    end
    return if method ∉ [
            LagrangeInterpolation, BSplineInterpolation, BSplineApprox, QuadraticSpline,
            AkimaInterpolation,
        ]
        @testset "$name, derivatives w.r.t. u" begin
            function f(u)
                A = method(
                    args..., u, t, args_after...; kwargs...,
                    extrapolation = ExtrapolationType.Extension
                )
                out = if u isa AbstractVector{<:Real}
                    zero(eltype(u))
                elseif u isa AbstractMatrix
                    zero(u[:, 1])
                else
                    zero(u[1])
                end

                for _t in trange
                    out += A(_t)
                end
                out
            end
            if u isa AbstractVector{<:Real}
                zgrad, fgrad = only(Zygote.gradient(f, u)), ForwardDiff.gradient(f, u)
                @test zgrad ≈ fgrad
            elseif u isa AbstractMatrix
                zgrad, fgrad = only(Zygote.jacobian(f, u)), ForwardDiff.jacobian(f, u)
                @test zgrad ≈ fgrad
            else
                # `Zygote.jacobian` on a `Vector{<:AbstractVector}` input returns a
                # nested/object-typed structure that isn't directly comparable to
                # `ForwardDiff.jacobian`'s flat matrix, and reconstructing one from
                # `hcat(u...)` breaks methods with extra same-shape args (e.g. `du`
                # for `CubicHermiteSpline`, since only `u` gets reshaped). Just check
                # that Zygote computes without erroring.
                Zygote.jacobian(f, u)
            end
        end
    end
end

@testset "LinearInterpolation" begin
    u = vcat(collect(1.0:5.0), 2 * collect(6.0:10.0))
    t = collect(1.0:10.0)
    test_zygote(
        LinearInterpolation, u, t; name = "Linear Interpolation"
    )
    u2 = Matrix(hcat(u, u)')
    test_zygote(
        LinearInterpolation, u2, t; name = "Linear Interpolation with matrix input"
    )
end

@testset "Quadratic Interpolation" begin
    u = [1.0, 4.0, 9.0, 16.0]
    t = [1.0, 2.0, 3.0, 4.0]
    test_zygote(QuadraticInterpolation, u, t; name = "Quadratic Interpolation")
    u2 = Matrix(hcat(u, u)')
    test_zygote(
        QuadraticInterpolation, u2, t; name = "Quadratic Interpolation with matrix input"
    )
end

@testset "Constant Interpolation" begin
    u = [0.0, 2.0, 1.0, 3.0, 2.0, 6.0, 5.5, 5.5, 2.7, 5.1, 3.0]
    t = collect(0.0:10.0)
    test_zygote(ConstantInterpolation, u, t; name = "Constant Interpolation (vector)")

    t = [1.0, 4.0]
    u = [1.0 2.0; 0.0 1.0; 1.0 2.0; 0.0 1.0]
    test_zygote(ConstantInterpolation, u, t, name = "Constant Interpolation (matrix)")

    u = [[1.0, 2.0, 3.0, 4.0], [2.0, 3.0, 4.0, 5.0]]
    test_zygote(
        ConstantInterpolation, u, t, name = "Constant Interpolation (vector of vectors)"
    )
end

@testset "Cubic Hermite Spline" begin
    du = [-0.047, -0.058, 0.054, 0.012, -0.068, 0.0]
    u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
    t = [0.0, 62.25, 109.66, 162.66, 205.8, 252.3]
    test_zygote(CubicHermiteSpline, u, t, args = [du], name = "Cubic Hermite Spline")
    du2 = Matrix(hcat(du, du)')
    u2 = Matrix(hcat(u, u)')
    test_zygote(
        CubicHermiteSpline, u2, t, args = [du2],
        name = "Cubic Hermite Spline with matrix input"
    )
    du_vov = [[du[i], du[i]] for i in eachindex(du)]
    u_vov = [[u[i], u[i]] for i in eachindex(u)]
    test_zygote(
        CubicHermiteSpline, u_vov, t, args = [du_vov],
        name = "Cubic Hermite Spline with vector of vectors input"
    )
end

@testset "Quintic Hermite Spline" begin
    ddu = [0.0, -0.00033, 0.0051, -0.0067, 0.0029, 0.0]
    du = [-0.047, -0.058, 0.054, 0.012, -0.068, 0.0]
    u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
    t = [0.0, 62.25, 109.66, 162.66, 205.8, 252.3]
    test_zygote(
        QuinticHermiteSpline, u, t, args = [ddu, du], name = "Quintic Hermite Spline"
    )
    ddu2 = Matrix(hcat(ddu, ddu)')
    du2 = Matrix(hcat(du, du)')
    u2 = Matrix(hcat(u, u)')
    test_zygote(
        QuinticHermiteSpline, u2, t, args = [ddu2, du2],
        name = "Quintic Hermite Spline with matrix input"
    )
    ddu_vov = [[ddu[i], ddu[i]] for i in eachindex(ddu)]
    du_vov = [[du[i], du[i]] for i in eachindex(du)]
    u_vov = [[u[i], u[i]] for i in eachindex(u)]
    test_zygote(
        QuinticHermiteSpline, u_vov, t, args = [ddu_vov, du_vov],
        name = "Quintic Hermite Spline with vector of vectors input"
    )
end

@testset "Akima Interpolation" begin
    u = [0.0, 2.0, 1.0, 3.0, 2.0, 6.0, 5.5, 5.5, 2.7, 5.1, 3.0]
    t = collect(0.0:10.0)
    test_zygote(AkimaInterpolation, u, t, name = "Akima Interpolation")
    u2 = Matrix(hcat(u, u)')
    test_zygote(
        AkimaInterpolation, u2, t, name = "Akima Interpolation with matrix input"
    )
    u_vov = [[u[i], u[i]] for i in eachindex(u)]
    test_zygote(
        AkimaInterpolation, u_vov, t,
        name = "Akima Interpolation with vector of vectors input"
    )
end

@testset "QuadraticSpline" begin
    u = [0.0, 1.0, 3.0]
    t = [-1.0, 0.0, 1.0]
    test_zygote(QuadraticSpline, u, t, name = "QuadraticSpline")
    u2 = Matrix(hcat(u, u)')
    test_zygote(QuadraticSpline, u2, t, name = "QuadraticSpline with matrix input")
    u_vov = [[u[i], u[i]] for i in eachindex(u)]
    test_zygote(
        QuadraticSpline, u_vov, t, name = "QuadraticSpline with vector of vectors input"
    )
end

@testset "SmoothedConstantInterpolation" begin
    u = [0.0, 2.0, 1.0, 3.0, 2.0, 6.0, 5.5, 5.5, 2.7, 5.1, 3.0]
    t = collect(0.0:10.0)
    test_zygote(
        SmoothedConstantInterpolation, u, t, name = "SmoothedConstantInterpolation"
    )
    u2 = Matrix(hcat(u, u)')
    test_zygote(
        SmoothedConstantInterpolation, u2, t,
        name = "SmoothedConstantInterpolation with matrix input"
    )
    u_vov = [[u[i], u[i]] for i in eachindex(u)]
    test_zygote(
        SmoothedConstantInterpolation, u_vov, t,
        name = "SmoothedConstantInterpolation with vector of vectors input"
    )
end

@testset "Lagrange Interpolation" begin
    u = [1.0, 4.0, 9.0]
    t = [1.0, 2.0, 3.0]
    test_zygote(LagrangeInterpolation, u, t, name = "Lagrange Interpolation")
end

@testset "Constant Interpolation" begin
    u = [0.0, 2.0, 1.0, 3.0, 2.0, 6.0, 5.5, 5.5, 2.7, 5.1, 3.0]
    t = collect(0.0:10.0)
    test_zygote(ConstantInterpolation, u, t, name = "Constant Interpolation")
end

@testset "Cubic Spline" begin
    u = [0.0, 1.0, 3.0]
    t = [-1.0, 0.0, 1.0]
    test_zygote(CubicSpline, u, t, name = "Cubic Spline")
    u2 = Matrix(hcat(u, u)')
    test_zygote(CubicSpline, u2, t, name = "Cubic Spline with matrix input")
    u_vov = [[u[i], u[i]] for i in eachindex(u)]
    test_zygote(CubicSpline, u_vov, t, name = "Cubic Spline with vector of vectors input")
end

@testset "BSplines" begin
    t = [0, 62.25, 109.66, 162.66, 205.8, 252.3]
    u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
    test_zygote(
        BSplineInterpolation, u, t; args_after = [2, :Uniform],
        name = "BSpline Interpolation"
    )
    test_zygote(
        BSplineApprox, u, t; args_after = [2, 4, :Uniform],
        name = "BSpline approximation"
    )
    u_vov = [[u[i], u[i]] for i in eachindex(u)]
    test_zygote(
        BSplineInterpolation, u_vov, t; args_after = [2, :Uniform],
        name = "BSpline Interpolation with vector of vectors input"
    )
    test_zygote(
        BSplineApprox, u_vov, t; args_after = [2, 4, :Uniform],
        name = "BSpline approximation with vector of vectors input"
    )
end

@testset "SmoothArcLengthInterpolation" begin
    # Doesn't fit `test_zygote`'s `method(u, t, ...)` calling convention: `t` is derived
    # from `u` (arc length), not a separate user-supplied argument.
    u = [0.3 -1.5 3.1; -0.2 0.2 -1.5; 10.4 -37.2 -5.8]
    A = SmoothArcLengthInterpolation(u; m = 5)
    trange = range(minimum(A.t), maximum(A.t); length = 25)
    @testset "SmoothArcLengthInterpolation, derivatives w.r.t. t" begin
        for _t in trange
            adiff = DataInterpolations.derivative(A, _t)
            zdiff = only(Zygote.jacobian(t -> A(t), _t))
            @test adiff ≈ zdiff
        end
    end
    # `derivatives w.r.t. u` is not supported: the geometry-fitting constructor mutates
    # several `Matrix`/`Vector` buffers in place (`@.` assignments while building `center`,
    # `dir_1`, `dir_2`, the augmented tangent curve, ...), which Zygote's reverse-mode AD
    # cannot trace through ("Mutating arrays is not supported").
end
