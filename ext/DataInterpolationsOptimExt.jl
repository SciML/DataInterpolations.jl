module DataInterpolationsOptimExt

if isdefined(Base, :get_extension)
    using DataInterpolations: AbstractInterpolation, munge_data
    import DataInterpolations: Curvefit, _interpolate
    using Reexport
else
    using ..DataInterpolations: AbstractInterpolation, munge_data
    import ..DataInterpolations: Curvefit, _interpolate
    using Reexport
end

isdefined(Base, :get_extension) ? (@reexport using Optim) : (@reexport using ..Optim)

### Curvefit
struct CurvefitCache{
    uType,
    tType,
    mType,
    p0Type,
    ubType,
    lbType,
    algType,
    pminType,
    FT,
    T,
} <: AbstractInterpolation{FT, T}
    u::uType
    t::tType
    m::mType        # model type
    p0::p0Type      # intial params
    ub::ubType      # upper bound of params
    lb::lbType      # lower bound of params
    alg::algType    # alg to optimize cost function
    pmin::pminType  # optimized params
    function CurvefitCache{FT}(u, t, m, p0, ub, lb, alg, pmin) where {FT}
        new{typeof(u), typeof(t), typeof(m),
            typeof(p0), typeof(ub), typeof(lb),
            typeof(alg), typeof(pmin), FT, eltype(u)}(u,
            t,
            m,
            p0,
            ub,
            lb,
            alg,
            pmin)
    end
end

function Curvefit(u, t, model, p0, alg, box = false, lb = nothing, ub = nothing)
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
    CurvefitCache{true}(u, t, model, p0, ub, lb, alg, pmin)
end

# Curvefit
function _interpolate(A::CurvefitCache{<:AbstractVector{<:Number}},
    t::Union{AbstractVector{<:Number}, Number})
    A.m(t, A.pmin)
end

function _interpolate(A::CurvefitCache{<:AbstractVector{<:Number}},
    t::Union{AbstractVector{<:Number}, Number},
    i)
    _interpolate(A, t), i
end

end # module
