abstract type AbstractIntegralInverseInterpolation{T} <: AbstractInterpolation{T} end

_integral(A::AbstractIntegralInverseInterpolation, idx, t) = throw(IntegralNotFoundError())
invert_integral(A::AbstractInterpolation) = throw(IntegralInverseNotFoundError())

struct LinearInterpolationIntInv{uType, tType, itpType, T} <:
       AbstractIntegralInverseInterpolation{T}
    u::uType
    t::tType
    extrapolate::Bool
    idx_prev::Base.RefValue{Int}
    itp::itpType
    function LinearInterpolationIntInv(u, t, A)
        new{typeof(u), typeof(t), typeof(A), eltype(u)}(
            u, t, A.extrapolate, Ref(1), A)
    end
end

function invertible_integral(A::LinearInterpolation{<:AbstractVector{<:Number}})
    return all(A.u .> 0)
end

function invert_integral(A::LinearInterpolation{<:AbstractVector{<:Number}})
    !invertible_integral(A) && throw(IntegralNotInvertibleError())
    return LinearInterpolationIntInv(A.t, A.I, A)
end

function _interpolate(
        A::LinearInterpolationIntInv{<:AbstractVector{<:Number}}, t::Number, iguess)
    idx = get_idx(A.t, t, iguess)
    Δt = t - A.t[idx]
    # TODO: Get from parameter cache: itp.p.slope[idx]
    slope = (A.itp.u[idx + 1] - A.itp.u[idx]) / (A.itp.t[idx + 1] - A.itp.t[idx])
    x = A.itp.u[idx]
    u = A.u[idx] + 2Δt / (x + sqrt(x^2 + 2slope * Δt))
    u, idx
end

struct ConstantInterpolationIntInv{uType, tType, itpType, T} <:
       AbstractIntegralInverseInterpolation{T}
    u::uType
    t::tType
    extrapolate::Bool
    idx_prev::Base.RefValue{Int}
    itp::itpType
    function ConstantInterpolationIntInv(u, t, A)
        new{typeof(u), typeof(t), typeof(A), eltype(u)}(
            u, t, A.extrapolate, Ref(1), A
        )
    end
end

function invertible_integral(A::ConstantInterpolation{<:AbstractVector{<:Number}})
    return all(A.u .> 0)
end

function invert_integral(A::ConstantInterpolation{<:AbstractVector{<:Number}})
    !invertible_integral(A) && throw(IntegralNotInvertibleError())
    return ConstantInterpolationIntInv(A.t, A.I, A)
end

function _interpolate(
        A::ConstantInterpolationIntInv{<:AbstractVector{<:Number}}, t::Number, iguess)
    idx = get_idx(A.t, t, iguess; ub_shift = 0)
    if A.itp.dir === :left
        # :left means that value to the left is used for interpolation
        idx_ = get_idx(A.t, t, idx; lb = 1, ub_shift = 0)
    else
        # :right means that value to the right is used for interpolation
        idx_ = get_idx(A.t, t, idx; side = :first, lb = 1, ub_shift = 0)
    end
    A.u[idx] + (t - A.t[idx]) / A.itp.u[idx_], idx
end
