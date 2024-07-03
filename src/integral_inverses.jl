abstract type AbstractIntegralInverseInterpolation{T} <: AbstractInterpolation{T} end

struct LinearInterpolationIntInv{tType, IType, itpType, T} <:
       AbstractIntegralInverseInterpolation{T}
    u::tType
    t::IType
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

invert_integral(A::AbstractInterpolation) = throw(IntegralInverseNotFoundError())

function invert_integral(A::LinearInterpolation{<:AbstractVector{<:Number}})
    !invertible_integral(A) && throw(IntegralNotInvertibleError())
    return LinearInterpolationIntInv(A.t, A.I, A)
end

function _interpolate(
        A::LinearInterpolationIntInv{<:AbstractVector{<:Number}}, t::Number, iguess)
    (; itp) = A
    idx = get_idx(A.t, t, iguess)
    Δt = t - A.t[idx]
    # TODO: Get from parameter cache: itp.p.slope[idx]
    slope = (itp.u[idx + 1] - itp.u[idx]) / (itp.t[idx + 1] - itp.t[idx])
    x = itp.u[idx]
    u = A.u[idx] + 2Δt / (x + sqrt(x^2 + 2slope * Δt))
    u, idx
end
