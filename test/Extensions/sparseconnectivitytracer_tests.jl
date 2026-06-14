using SparseConnectivityTracer
using SparseConnectivityTracer: DEFAULT_GRADIENT_TRACER, DEFAULT_HESSIAN_TRACER
using SparseConnectivityTracer: trace_input, Dual, primal
using DataInterpolations
using DataInterpolations: AbstractInterpolation

using LinearAlgebra: I
using Test

myprimal(x) = x
myprimal(d::Dual) = primal(d)

#===========#
# Test data #
#===========#

t = [0.0, 1.0, 2.5, 4.0, 6.0];
t_scalar = 2.0
t_vector = [2.0, 2.5, 3.0]
t_range = 2:5
test_inputs = (t_scalar, t_vector, t_range)

u = sin.(t) # vector
du = cos.(t)
ddu = -sin.(t)

#==================#
# Test definitions #
#==================#

struct InterpolationTest{N, I <: AbstractInterpolation} # N = output dim. of interpolation
    interp::I
    is_der1_zero::Bool
    is_der2_zero::Bool
end
function InterpolationTest(
        interp::I; is_der1_zero = false, is_der2_zero = false
    ) where {T, I <: AbstractInterpolation{T}}
    sz = output_size(interp)
    N = isempty(sz) ? 1 : only(sz)
    return InterpolationTest{N, I}(interp, is_der1_zero, is_der2_zero)
end
testname(t::InterpolationTest{N}) where {N} = "$N-dim $(typeof(t.interp))"

#================#
# Jacobian Tests #
#================#

function test_jacobian(t::InterpolationTest)
    return @testset "Jacobian" begin
        for input in test_inputs
            test_jacobian(t, input)
        end
    end
end
function test_jacobian(t::InterpolationTest{N}, input::Real) where {N}
    N_IN = length(input)
    N_OUT = N * N_IN
    Jref = t.is_der1_zero ? zeros(N, N_IN) : ones(N, N_IN)

    return @testset "input type $(typeof(input)): $N_IN inputs, $N states, $N_OUT outputs" begin
        @testset "Global Jacobian sparsity" begin
            J = jacobian_sparsity(t.interp, input, TracerSparsityDetector())
            @test J ≈ Jref
        end
        @testset "Local Jacobian sparsity" begin
            J = jacobian_sparsity(t.interp, input, TracerLocalSparsityDetector())
            @test J ≈ Jref
        end
    end
end
function test_jacobian(t::InterpolationTest{1}, input::AbstractVector)
    N = 1
    N_IN = length(input)
    N_OUT = N * N_IN
    Jref = t.is_der1_zero ? zeros(N_IN, N_IN) : I(N_IN)

    return @testset "input type $(typeof(input)): $N_IN inputs, $N states, $N_OUT outputs" begin
        @testset "Global Jacobian sparsity" begin
            J = jacobian_sparsity(x -> vec(t.interp(x)), input, TracerSparsityDetector())
            @test J ≈ Jref
        end
        @testset "Local Jacobian sparsity" begin
            J = jacobian_sparsity(
                x -> vec(t.interp(x)), input, TracerLocalSparsityDetector()
            )
            @test J ≈ Jref
        end
    end
end
function test_jacobian(t::InterpolationTest{N}, input::AbstractVector) where {N}
    N_IN = length(input)
    N_OUT = N * N_IN

    # Construct reference Jacobian
    Jref = zeros(Bool, N_OUT, N_IN)
    if !t.is_der1_zero
        for (i, col) in enumerate(eachcol(Jref)) # iterate over outputs
            i0 = 1 + N * (i - 1)
            irange = i0:(i0 + N - 1)
            col[irange] .= true
        end
    end

    return @testset "input type $(typeof(input)): $N_IN inputs, $N states, $N_OUT outputs" begin
        @testset "Global Jacobian sparsity" begin
            J = jacobian_sparsity(x -> vec(t.interp(x)), input, TracerSparsityDetector())
            @test J ≈ Jref
        end
        @testset "Local Jacobian sparsity" begin
            J = jacobian_sparsity(
                x -> vec(t.interp(x)), input, TracerLocalSparsityDetector()
            )
            @test J ≈ Jref
        end
    end
end

#===============#
# Hessian Tests #
#===============#

function test_hessian(t::InterpolationTest)
    return @testset "Hessian" begin
        for input in test_inputs
            test_hessian(t, input)
        end
    end
end
function test_hessian(t::InterpolationTest{1}, input::Real)
    N = 1
    N_IN = length(input)
    N_OUT = N * N_IN
    Href = t.is_der2_zero ? zeros(N_IN, N_IN) : ones(N_IN, N_IN)

    return @testset "input type $(typeof(input)): $N_IN inputs, $N states, $N_OUT outputs" begin
        @testset "Global Hessian sparsity" begin
            H = hessian_sparsity(t.interp, input, TracerSparsityDetector())
            @test H ≈ Href
        end
        @testset "Local Hessian sparsity" begin
            H = hessian_sparsity(t.interp, input, TracerLocalSparsityDetector())
            @test H ≈ Href
        end
    end
end
function test_hessian(t::InterpolationTest{N}, input::Real) where {N} #  N ≠ 1
    N_IN = length(input)
    N_OUT = N * N_IN
    Href = t.is_der2_zero ? zeros(N_IN, N_IN) : ones(N_IN, N_IN)

    return @testset "input type $(typeof(input)): $N_IN inputs, $N states, $N_OUT outputs" begin
        @testset "Global Hessian sparsity" begin
            H = hessian_sparsity(x -> sum(t.interp(x)), input, TracerSparsityDetector())
            @test H ≈ Href
        end
        @testset "Local Hessian sparsity" begin
            H = hessian_sparsity(
                x -> sum(t.interp(x)), input, TracerLocalSparsityDetector()
            )
            @test H ≈ Href
        end
    end
end
function test_hessian(t::InterpolationTest{1}, input::AbstractVector)
    N = 1
    N_IN = length(input)
    N_OUT = N * N_IN
    Href = t.is_der2_zero ? zeros(N_IN, N_IN) : I(N_IN)

    return @testset "input type $(typeof(input)): $N_IN inputs, $N states, $N_OUT outputs" begin
        @testset "Global Hessian sparsity" begin
            H = hessian_sparsity(x -> sum(t.interp(x)), input, TracerSparsityDetector())
            @test H ≈ Href
        end
        @testset "Local Hessian sparsity" begin
            H = hessian_sparsity(
                x -> sum(t.interp(x)), input, TracerLocalSparsityDetector()
            )
            @test H ≈ Href
        end
    end
end
function test_hessian(t::InterpolationTest{N}, input::AbstractVector) where {N} #  N ≠ 1
    N_IN = length(input)
    N_OUT = N * N_IN
    Href = t.is_der2_zero ? zeros(N_IN, N_IN) : I(N_IN)

    return @testset "input type $(typeof(input)): $N_IN inputs, $N states, $N_OUT outputs" begin
        @testset "Global Hessian sparsity" begin
            H = hessian_sparsity(x -> sum(t.interp(x)), input, TracerSparsityDetector())
            @test H ≈ Href
        end
        @testset "Local Hessian sparsity" begin
            H = hessian_sparsity(
                x -> sum(t.interp(x)), input, TracerLocalSparsityDetector()
            )
            @test H ≈ Href
        end
    end
end

function test_output(t::InterpolationTest)
    return @testset "Output sizes and values" begin
        @testset "input type: $(typeof(input))" for input in test_inputs
            out_ref = t.interp(input)
            s_ref = size(out_ref)

            @testset "$T" for T in (DEFAULT_GRADIENT_TRACER, DEFAULT_HESSIAN_TRACER)
                t_tracer = trace_input(T, input)
                out_tracer = t.interp(t_tracer)
                s_tracer = size(out_tracer)
                @test s_tracer == s_ref
            end
            @testset "$T" for T in (
                    Dual{eltype(input), DEFAULT_GRADIENT_TRACER},
                    Dual{eltype(input), DEFAULT_HESSIAN_TRACER},
                )
                t_dual = trace_input(T, input)
                out_dual = t.interp(t_dual)
                s_dual = size(out_dual)
                @test s_dual == s_ref
                @test myprimal.(out_dual) ≈ out_ref
            end
        end
    end
end

#===========#
# Run tests #
#===========#

@testset "1D Interpolations" begin
    @testset "$(testname(t))" for t in (
            InterpolationTest(
                ConstantInterpolation(u, t); is_der1_zero = true, is_der2_zero = true
            ),
            InterpolationTest(LinearInterpolation(u, t); is_der2_zero = true),
            InterpolationTest(QuadraticInterpolation(u, t)),
            InterpolationTest(LagrangeInterpolation(u, t)),
            InterpolationTest(AkimaInterpolation(u, t)),
            InterpolationTest(QuadraticSpline(u, t)),
            InterpolationTest(CubicSpline(u, t)),
            InterpolationTest(BSplineInterpolation(u, t, 3, :ArcLen, :Average)),
            InterpolationTest(BSplineApprox(u, t, 3, 4, :ArcLen, :Average)),
            InterpolationTest(PCHIPInterpolation(u, t)),
            InterpolationTest(CubicHermiteSpline(du, u, t)),
            InterpolationTest(QuinticHermiteSpline(ddu, du, u, t)),
        )
        test_jacobian(t)
        test_hessian(t)
        test_output(t)
        yield()
    end
end

for N in (2, 5)
    local um = rand(N, length(t)) # matrix

    @testset "$(N)D Interpolations" begin
        @testset "$(testname(t))" for t in (
                InterpolationTest(
                    ConstantInterpolation(um, t); is_der1_zero = true, is_der2_zero = true
                ),
                InterpolationTest(LinearInterpolation(um, t); is_der2_zero = true),
                InterpolationTest(QuadraticInterpolation(um, t)),
                InterpolationTest(LagrangeInterpolation(um, t)),
                ## The following interpolations appear to not be supported on N dimensions as of DataInterpolations v6.2.0:
                # InterpolationTest(AkimaInterpolation(um, t)),
                # InterpolationTest(BSplineApprox(um, t, 3, 4, :ArcLen, :Average)),
                # InterpolationTest(QuadraticSpline(um, t)),
                # InterpolationTest(CubicSpline(um, t)),
                # InterpolationTest(BSplineInterpolation(um, t, 3, :ArcLen, :Average)),
                # InterpolationTest(PCHIPInterpolation(um, t)),
            )
            test_jacobian(t)
            test_hessian(t)
            test_output(t)
            yield()
        end
    end
end
