module DataInterpolationsCurveFitExt

using DataInterpolations
import DataInterpolations: munge_data,
    Curvefit, CurvefitCache, _interpolate, get_show, derivative,
    ExtrapolationError,
    integral, IntegralNotFoundError, DerivativeNotFoundError

using CurveFit, ForwardDiff

### Curvefit
function Curvefit(
        u,
        t,
        model,
        p0,
        alg = nothing,
        box = false,
        lb = nothing,
        ub = nothing;
        extrapolate = false
    )
    u, t = munge_data(u, t)
    if box && (lb === nothing || ub === nothing)
        error("lower or upper bound should not be nothing")
    end
    prob = CurveFit.NonlinearCurveFitProblem((p, x) -> model(x, p), p0, t, u)
    sol = if isnothing(alg)
        CurveFit.solve(prob)
    else
        CurveFit.solve(prob, alg)
    end
    pmin = sol.u
    CurvefitCache(u, t, model, p0, ub, lb, alg, pmin, extrapolate)
end

# Curvefit
function _interpolate(
        A::CurvefitCache{<:AbstractVector{<:Number}},
        t::Union{AbstractVector{<:Number}, Number}
    )
    ((t < A.t[1] || t > A.t[end]) && !A.extrapolate) &&
        throw(ExtrapolationError())
    return A.m(t, A.pmin)
end

function _interpolate(
        A::CurvefitCache{<:AbstractVector{<:Number}},
        t::Union{AbstractVector{<:Number}, Number},
        i
    )
    return _interpolate(A, t), i
end

function derivative(
        A::CurvefitCache{<:AbstractVector{<:Number}},
        t::Union{AbstractVector{<:Number}, Number}, order = 1
    )
    ((t < A.t[1] || t > A.t[end]) && !A.extrapolate) && throw(ExtrapolationError())
    order > 2 && throw(DerivativeNotFoundError())
    order == 1 && return ForwardDiff.derivative(x -> A.m(x, A.pmin), t)
    return ForwardDiff.derivative(t -> ForwardDiff.derivative(x -> A.m(x, A.pmin), t), t)
end

function get_show(A::CurvefitCache)
    return "Curvefit" *
        " with $(length(A.t)) points.\n"
end

function integral(A::CurvefitCache{<:AbstractVector{<:Number}}, t::Number)
    throw(IntegralNotFoundError())
end

function integral(A::CurvefitCache{<:AbstractVector{<:Number}}, t1::Number, t2::Number)
    throw(IntegralNotFoundError())
end

end # module
