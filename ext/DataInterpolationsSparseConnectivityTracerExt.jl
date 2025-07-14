module DataInterpolationsSparseConnectivityTracerExt

using SparseConnectivityTracer: AbstractTracer, Dual, primal, tracer
using SparseConnectivityTracer: GradientTracer, gradient_tracer_1_to_1
using SparseConnectivityTracer: HessianTracer, hessian_tracer_1_to_1
using FillArrays: Fill # from FillArrays.jl
using DataInterpolations:
    AbstractInterpolation,
    LinearInterpolation,
    QuadraticInterpolation,
    LagrangeInterpolation,
    AkimaInterpolation,
    ConstantInterpolation,
    QuadraticSpline,
    CubicSpline,
    BSplineInterpolation,
    BSplineApprox,
    CubicHermiteSpline,
    # PCHIPInterpolation,
    QuinticHermiteSpline,
    output_size

#===========#
# Utilities #
#===========#

# Limit support to `u` begin an AbstractVector{<:Number} or AbstractMatrix{<:Number},
# to avoid any cases where the output size is dependent on the input value.
# https://github.com/adrhill/SparseConnectivityTracer.jl/pull/234#discussion_r2031038566

function _sct_interpolate(
        ::AbstractInterpolation,
        uType::Type{<:AbstractVector{<:Number}},
        t::GradientTracer,
        is_der_1_zero,
        is_der_2_zero,
    )
    return gradient_tracer_1_to_1(t, is_der_1_zero)
end
function _sct_interpolate(
        ::AbstractInterpolation,
        uType::Type{<:AbstractVector{<:Number}},
        t::HessianTracer,
        is_der_1_zero,
        is_der_2_zero,
    )
    return hessian_tracer_1_to_1(t, is_der_1_zero, is_der_2_zero)
end
function _sct_interpolate(
        interp::AbstractInterpolation,
        uType::Type{<:AbstractMatrix{<:Number}},
        t::GradientTracer,
        is_der_1_zero,
        is_der_2_zero,
    )
    t = gradient_tracer_1_to_1(t, is_der_1_zero)
    N = only(output_size(interp))
    return Fill(t, N)
end
function _sct_interpolate(
        interp::AbstractInterpolation,
        uType::Type{<:AbstractMatrix{<:Number}},
        t::HessianTracer,
        is_der_1_zero,
        is_der_2_zero,
    )
    t = hessian_tracer_1_to_1(t, is_der_1_zero, is_der_2_zero)
    N = only(output_size(interp))
    return Fill(t, N)
end

#===========#
# Overloads #
#===========#

# We assume that with the exception of ConstantInterpolation and LinearInterpolation,
# all interpolations have a non-zero second derivative at some point in the input domain.

for (I, is_der1_zero, is_der2_zero) in (
        (:ConstantInterpolation, true, true),
        (:LinearInterpolation, false, true),
        (:QuadraticInterpolation, false, false),
        (:LagrangeInterpolation, false, false),
        (:AkimaInterpolation, false, false),
        (:QuadraticSpline, false, false),
        (:CubicSpline, false, false),
        (:BSplineInterpolation, false, false),
        (:BSplineApprox, false, false),
        (:CubicHermiteSpline, false, false),
        (:QuinticHermiteSpline, false, false),
    )
    @eval function (interp::$(I){uType})(
            t::AbstractTracer
        ) where {uType <: AbstractArray{<:Number}}
        return _sct_interpolate(interp, uType, t, $is_der1_zero, $is_der2_zero)
    end
end

# Some Interpolations require custom overloads on `Dual` due to mutation of caches.
for I in (
        :LagrangeInterpolation,
        :BSplineInterpolation,
        :BSplineApprox,
        :CubicHermiteSpline,
        :QuinticHermiteSpline,
    )
    @eval function (interp::$(I){uType})(d::Dual) where {uType <: AbstractVector}
        p = interp(primal(d))
        t = interp(tracer(d))
        return Dual(p, t)
    end

    @eval function (interp::$(I){uType})(d::Dual) where {uType <: AbstractMatrix}
        p = interp(primal(d))
        t = interp(tracer(d))
        return Dual.(p, t)
    end
end

end
