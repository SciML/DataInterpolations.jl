module DataInterpolationsOptimExt

using DataInterpolations
import DataInterpolations: munge_data,
    Curvefit, CurvefitCache, _interpolate, get_show, derivative, ExtrapolationError

isdefined(Base, :get_extension) ? (using Optim, ForwardDiff) :
(using ..Optim, ..ForwardDiff)

### Curvefit
function Curvefit(u,
        t,
        model,
        p0,
        alg,
        box = false,
        lb = nothing,
        ub = nothing;
        extrapolate = false)
    u, t = munge_data(u, t)
    errfun(t, u, p) = sum(abs2.(u .- model(t, p)))
    if box == false
        mfit = optimize(p -> errfun(t, u, p), p0, alg)
    else
        if lb === nothing || ub === nothing
            error("lower or upper bound should not be nothing")
        end
        od = OnceDifferentiable(p -> errfun(t, u, p), p0, autodiff = :finite)
        mfit = optimize(od, lb, ub, p0, Fminbox(alg))
    end
    pmin = Optim.minimizer(mfit)
    CurvefitCache{true}(u, t, model, p0, ub, lb, alg, pmin, extrapolate)
end

# Curvefit
function _interpolate(A::CurvefitCache{<:AbstractVector{<:Number}},
        t::Union{AbstractVector{<:Number}, Number})
    ((t < A.t[1] || t > A.t[end]) && !A.extrapolate) &&
        throw(ExtrapolationError())
    A.m(t, A.pmin)
end

function _interpolate(A::CurvefitCache{<:AbstractVector{<:Number}},
        t::Union{AbstractVector{<:Number}, Number},
        i)
    _interpolate(A, t), i
end

function derivative(A::CurvefitCache{<:AbstractVector{<:Number}},
        t::Union{AbstractVector{<:Number}, Number})
    ((t < A.t[1] || t > A.t[end]) && !A.extrapolate) && throw(ExtrapolationError())
    ForwardDiff.derivative(x -> A.m(x, A.pmin), t)
end

function get_show(interp::CurvefitCache)
    return "Curvefit" *
           " with $(length(interp.t)) points, using $(nameof(typeof(interp.alg)))\n"
end

end # module
