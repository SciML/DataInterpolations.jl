module DataInterpolationsMooncakeExt

using DataInterpolations, Mooncake, ChainRulesCore, FindFirstFunctions
using DataInterpolations: _interpolate, munge_data, AbstractInterpolation,
    LinearInterpolation, QuadraticInterpolation
import Mooncake: @from_chainrules, @zero_adjoint, MinimalCtx, DefaultCtx

# When the ChainRules pullback for _interpolate returns a Tangent{AbstractInterpolation},
# this tells Mooncake how to accumulate the u-component into the interpolation's fdata.
function Mooncake.increment_and_get_rdata!(
        f::Mooncake.FData{<:NamedTuple},
        r::Mooncake.NoRData,
        t::ChainRulesCore.Tangent{<:AbstractInterpolation}
    )
    u_tang = ChainRulesCore.unthunk(t.u)
    if !(u_tang isa ChainRulesCore.AbstractZero)
        f.data.u .+= u_tang
    end
    return Mooncake.NoRData()
end

# Same, when the cache's bitstype `t_props` scalars (`first_val`/`inv_step`)
# make Mooncake hand back a populated `RData{NamedTuple}`: those fields are
# construction-time constants, so accumulate `u` and return the rdata unchanged.
function Mooncake.increment_and_get_rdata!(
        f::Mooncake.FData{<:NamedTuple},
        r::Mooncake.RData{<:NamedTuple},
        t::ChainRulesCore.Tangent{<:AbstractInterpolation}
    )
    u_tang = ChainRulesCore.unthunk(t.u)
    if !(u_tang isa ChainRulesCore.AbstractZero)
        f.data.u .+= u_tang
    end
    return r
end

# Constructor rules: stop Mooncake recursing into LinearParameterCache and other
# internal structs. The 6-arg and 7-arg forms are the internal constructors that
# have ChainRules rrules defined in DataInterpolationsChainRulesCoreExt.
@from_chainrules MinimalCtx Tuple{Type{LinearInterpolation}, Any, Any, Any, Any, Any, Any} true
@from_chainrules MinimalCtx Tuple{Type{QuadraticInterpolation}, Any, Any, Any, Any, Any, Any, Any} true

# _interpolate: the core computation for all interpolation calls (A(t) dispatches here)
@from_chainrules MinimalCtx Tuple{typeof(_interpolate), LinearInterpolation, Number} true
@from_chainrules MinimalCtx Tuple{typeof(_interpolate), QuadraticInterpolation, Number} true

# munge_data: validates/reshapes u and t - identity in the non-missing case.
# Match all three method dispatches that exist in DataInterpolations.
@from_chainrules MinimalCtx Tuple{typeof(munge_data), AbstractVector, AbstractVector} true
@from_chainrules MinimalCtx Tuple{typeof(munge_data), AbstractMatrix, AbstractVector} true
@from_chainrules MinimalCtx Tuple{typeof(munge_data), AbstractArray, Any} true

# `get_idx` calls these (bare `StrategyKind`, or `Auto` for the uniform path);
# they return integer indices, so zero-adjoint cuts the gradient at the index
# boundary and stops Mooncake recursing into FFF's `llvmcall` SIMD kernels.
@zero_adjoint DefaultCtx Tuple{typeof(FindFirstFunctions.searchsorted_last), FindFirstFunctions.StrategyKind, AbstractVector, Any, Integer}
@zero_adjoint DefaultCtx Tuple{typeof(FindFirstFunctions.searchsorted_first), FindFirstFunctions.StrategyKind, AbstractVector, Any, Integer}
@zero_adjoint DefaultCtx Tuple{typeof(FindFirstFunctions.searchsorted_last), FindFirstFunctions.StrategyKind, AbstractVector, Any}
@zero_adjoint DefaultCtx Tuple{typeof(FindFirstFunctions.searchsorted_first), FindFirstFunctions.StrategyKind, AbstractVector, Any}
@zero_adjoint DefaultCtx Tuple{typeof(FindFirstFunctions.searchsorted_last), FindFirstFunctions.Auto, AbstractVector, Any, Integer}
@zero_adjoint DefaultCtx Tuple{typeof(FindFirstFunctions.searchsorted_first), FindFirstFunctions.Auto, AbstractVector, Any, Integer}
@zero_adjoint DefaultCtx Tuple{typeof(FindFirstFunctions.searchsorted_last), FindFirstFunctions.Auto, AbstractVector, Any}
@zero_adjoint DefaultCtx Tuple{typeof(FindFirstFunctions.searchsorted_first), FindFirstFunctions.Auto, AbstractVector, Any}

end
