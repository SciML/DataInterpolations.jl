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

# `SearchProperties{T}` and `Auto{T}` now carry `first_val::T` / `inv_step::T`
# (precomputed scalar fields for the props-aware UniformStep kernel). These
# fields show up in Mooncake's rdata for the interpolation cache because
# they are bitstype `Float64` (or `Float32`), but they are **not**
# differentiable — they are constants attached to the cache at construction.
# Mooncake doesn't know that from the type alone, so it routes the
# `Tangent{<:AbstractInterpolation}` (which has only `u` populated) through
# `increment_and_get_rdata!` with a populated `RData{NamedTuple{...}}` rather
# than `NoRData`. Tell Mooncake the right thing: accumulate `u`, then return
# the rdata unchanged so its t_props / strategy slots accumulate zero.
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

# Sorted-search dispatched through `Auto{T}` carries the props' `first_val::T`
# and `inv_step::T` Float fields, which Mooncake exposes as rdata. The
# `searchsorted_last` / `searchsorted_first` calls return integer indices — they are
# positional bookkeeping, not differentiable. Declare them as
# zero-adjoint so Mooncake doesn't try to recurse into FFF's strategy
# kernels (which contain `llvmcall` SIMD intrinsics that Mooncake cannot
# differentiate through). DI's interpolation `_interpolate` always feeds
# the search results into integer indexing, so the gradient flow is
# already cut at the index boundary — zero-adjoint here is correct.
@zero_adjoint DefaultCtx Tuple{typeof(FindFirstFunctions.searchsorted_last), FindFirstFunctions.Auto, AbstractVector, Any, Integer}
@zero_adjoint DefaultCtx Tuple{typeof(FindFirstFunctions.searchsorted_first), FindFirstFunctions.Auto, AbstractVector, Any, Integer}
@zero_adjoint DefaultCtx Tuple{typeof(FindFirstFunctions.searchsorted_last), FindFirstFunctions.Auto, AbstractVector, Any}
@zero_adjoint DefaultCtx Tuple{typeof(FindFirstFunctions.searchsorted_first), FindFirstFunctions.Auto, AbstractVector, Any}

end
