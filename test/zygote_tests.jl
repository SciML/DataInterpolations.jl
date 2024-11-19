using DataInterpolations
using ForwardDiff
using Zygote

function test_zygote(method, u, t; args = [], args_after = [], kwargs = [], name::String)
    func = method(args..., u, t, args_after...; kwargs...,
        extrapolation = ExtrapolationType.extension)
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
    if method ∉
       [LagrangeInterpolation, BSplineInterpolation, BSplineApprox, QuadraticSpline]
        @testset "$name, derivatives w.r.t. u" begin
            function f(u)
                A = method(args..., u, t, args_after...; kwargs...,
                    extrapolation = ExtrapolationType.extension)
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
            zgrad, fgrad = if u isa AbstractVector{<:Real}
                Zygote.gradient(f, u), ForwardDiff.gradient(f, u)
            elseif u isa AbstractMatrix
                Zygote.jacobian(f, u), ForwardDiff.jacobian(f, u)
            else
                Zygote.jacobian(f, u), ForwardDiff.jacobian(f, hcat(u...))
            end
        end
    end
end

@testset "LinearInterpolation" begin
    u = vcat(collect(1.0:5.0), 2 * collect(6.0:10.0))
    t = collect(1.0:10.0)
    test_zygote(
        LinearInterpolation, u, t; name = "Linear Interpolation")
end

@testset "Quadratic Interpolation" begin
    u = [1.0, 4.0, 9.0, 16.0]
    t = [1.0, 2.0, 3.0, 4.0]
    test_zygote(QuadraticInterpolation, u, t; name = "Quadratic Interpolation")
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
        ConstantInterpolation, u, t, name = "Constant Interpolation (vector of vectors)")
end

@testset "Cubic Hermite Spline" begin
    du = [-0.047, -0.058, 0.054, 0.012, -0.068, 0.0]
    u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
    t = [0.0, 62.25, 109.66, 162.66, 205.8, 252.3]
    test_zygote(CubicHermiteSpline, u, t, args = [du], name = "Cubic Hermite Spline")
end

@testset "Quintic Hermite Spline" begin
    ddu = [0.0, -0.00033, 0.0051, -0.0067, 0.0029, 0.0]
    du = [-0.047, -0.058, 0.054, 0.012, -0.068, 0.0]
    u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
    t = [0.0, 62.25, 109.66, 162.66, 205.8, 252.3]
    test_zygote(
        QuinticHermiteSpline, u, t, args = [ddu, du], name = "Quintic Hermite Spline")
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
end

@testset "BSplines" begin
    t = [0, 62.25, 109.66, 162.66, 205.8, 252.3]
    u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]
    test_zygote(BSplineInterpolation, u, t; args_after = [2, :Uniform, :Uniform],
        name = "BSpline Interpolation")
    test_zygote(BSplineApprox, u, t; args_after = [2, 4, :Uniform, :Uniform],
        name = "BSpline approximation")
end
