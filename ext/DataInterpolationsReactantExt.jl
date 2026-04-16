module DataInterpolationsReactantExt

using DataInterpolations: DataInterpolations,
    AbstractInterpolation,
    LinearInterpolation,
    QuadraticInterpolation,
    ConstantInterpolation,
    ExtrapolationType,
    get_parameters,
    linear_interpolation_parameters
import DataInterpolations: _interpolate
using Reactant: Reactant, TracedRNumber

# =============================================================================
# Branchless interpolation for Reactant traced numbers.
#
# When Reactant traces through an ODE function, the time `t` becomes a
# TracedRNumber. Standard DataInterpolations methods use `if t < ...` which
# fails because TracedRNumber{Bool} can't be used in boolean context.
#
# Strategy: evaluate the interpolation formula for all segments and use
# `ifelse` (which Reactant lowers to stablehlo.select) to pick the right one.
# This is O(n) in knot count but fully traceable.
# =============================================================================

# ===========================================================================
# Helper: get slope for a given segment index (concrete integer)
# ===========================================================================

function _get_slope(A::LinearInterpolation, idx::Int)
    if A.cache_parameters
        return A.p.slope[idx]
    else
        return linear_interpolation_parameters(A.u, A.t, idx)
    end
end

# ===========================================================================
# LinearInterpolation — branchless
# ===========================================================================

function _interpolate(A::LinearInterpolation{<:AbstractVector}, t::TracedRNumber)
    n = length(A.t)

    # Evaluate first segment
    u1 = oftype(t, A.u[1])
    slope1 = oftype(t, _get_slope(A, 1))
    t1 = oftype(t, A.t[1])
    result = u1 + slope1 * (t - t1)

    # Cascade through remaining segments: if t >= A.t[i], use segment i
    for i in 2:(n - 1)
        ui = oftype(t, A.u[i])
        slope_i = oftype(t, _get_slope(A, i))
        ti = oftype(t, A.t[i])
        seg_val = ui + slope_i * (t - ti)
        result = ifelse(t >= ti, seg_val, result)
    end

    # Handle extrapolation
    left_val = _traced_left_extrap_val(A, t)
    right_val = _traced_right_extrap_val(A, t)
    tmin = oftype(t, first(A.t))
    tmax = oftype(t, last(A.t))

    result = ifelse(t < tmin, left_val, result)
    result = ifelse(t > tmax, right_val, result)

    return result
end

function _traced_left_extrap_val(A::LinearInterpolation, t)
    ext = A.extrapolation_left
    if ext == ExtrapolationType.Constant
        return oftype(t, first(A.u))
    else
        # Linear, Extension, or fallback: extend with first segment's slope
        slope = oftype(t, _get_slope(A, 1))
        return oftype(t, first(A.u)) + slope * (t - oftype(t, first(A.t)))
    end
end

function _traced_right_extrap_val(A::LinearInterpolation, t)
    ext = A.extrapolation_right
    n = length(A.t)
    if ext == ExtrapolationType.Constant
        return oftype(t, last(A.u))
    else
        slope = oftype(t, _get_slope(A, n - 1))
        return oftype(t, last(A.u)) + slope * (t - oftype(t, last(A.t)))
    end
end

# ===========================================================================
# ConstantInterpolation — branchless
# ===========================================================================

function _interpolate(A::ConstantInterpolation{<:AbstractVector}, t::TracedRNumber)
    n = length(A.t)

    if A.dir === :left
        # :left — use value at the largest t[i] <= t
        result = oftype(t, A.u[1])
        for i in 2:n
            ti = oftype(t, A.t[i])
            result = ifelse(t >= ti, oftype(t, A.u[i]), result)
        end
    else
        # :right — use value at the smallest t[i] >= t
        result = oftype(t, A.u[n])
        for i in (n - 1):-1:1
            ti = oftype(t, A.t[i])
            result = ifelse(t <= ti, oftype(t, A.u[i]), result)
        end
    end

    # Handle extrapolation
    tmin = oftype(t, first(A.t))
    tmax = oftype(t, last(A.t))

    ext_left = A.extrapolation_left
    left_val = if ext_left == ExtrapolationType.Constant || ext_left == ExtrapolationType.Extension
        oftype(t, first(A.u))
    else
        oftype(t, first(A.u))
    end

    ext_right = A.extrapolation_right
    right_val = if ext_right == ExtrapolationType.Constant || ext_right == ExtrapolationType.Extension
        oftype(t, last(A.u))
    else
        oftype(t, last(A.u))
    end

    result = ifelse(t < tmin, left_val, result)
    result = ifelse(t > tmax, right_val, result)

    return result
end

# ===========================================================================
# QuadraticInterpolation — branchless
# ===========================================================================

function _interpolate(A::QuadraticInterpolation{<:AbstractVector}, t::TracedRNumber)
    n = length(A.t)

    # First segment
    α1, β1 = get_parameters(A, 1)
    Δt = t - oftype(t, A.t[1])
    result = oftype(t, A.u[1]) + Δt * (oftype(t, α1) * Δt + oftype(t, β1))

    # Remaining segments
    for i in 2:(n - 2)
        αi, βi = get_parameters(A, i)
        Δti = t - oftype(t, A.t[i])
        seg_val = oftype(t, A.u[i]) + Δti * (oftype(t, αi) * Δti + oftype(t, βi))
        ti = oftype(t, A.t[i])
        result = ifelse(t >= ti, seg_val, result)
    end

    # Handle extrapolation
    tmin = oftype(t, first(A.t))
    tmax = oftype(t, last(A.t))

    ext_left = A.extrapolation_left
    left_val = if ext_left == ExtrapolationType.Constant
        oftype(t, first(A.u))
    else
        α, β = get_parameters(A, 1)
        Δt_left = t - oftype(t, A.t[1])
        oftype(t, A.u[1]) + Δt_left * (oftype(t, α) * Δt_left + oftype(t, β))
    end

    ext_right = A.extrapolation_right
    last_idx = max(n - 2, 1)
    right_val = if ext_right == ExtrapolationType.Constant
        oftype(t, last(A.u))
    else
        α, β = get_parameters(A, last_idx)
        Δt_right = t - oftype(t, A.t[last_idx])
        oftype(t, A.u[last_idx]) + Δt_right * (oftype(t, α) * Δt_right + oftype(t, β))
    end

    result = ifelse(t < tmin, left_val, result)
    result = ifelse(t > tmax, right_val, result)

    return result
end

end # module
