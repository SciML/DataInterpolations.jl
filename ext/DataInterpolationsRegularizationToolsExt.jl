module DataInterpolationsRegularizationToolsExt

using DataInterpolations
using DataInterpolations: munge_data, _interpolate, RegularizationSmooth
using LinearAlgebra

isdefined(Base, :get_extension) ? (import RegularizationTools as RT) :
(import ..RegularizationTools as RT)

# TODO:
# x midpoint rule
# x scattered/interpolation
# x GCV
# x L-curve
# - bounds on λ
# - initial guess for λ? will require mods to RegularizationTools
# - scaled λ? will need to work out equivalency with λ² formulation, and resolve with
#   derivative rather than difference matrix
# - optimize λ via standard deviation?
# - relative weights?
# - arbitrary weighting -- implemented but not yet tested
# - midpoint rule with scattered?
# x add argument types for `RegularizationSmooth` constructor methods (why isn't this done
#   for the other interpolaters?)
# - make use of `munge_data` features (allow for matrix rather than vector u & t arguments?)
# - validate data and t̂
# x unit tests

const LA = LinearAlgebra

"""
# Arguments
- `u::Vector`:  dependent data
- `t::Vector`:  independent data

# Optional Arguments
- `t̂::Vector`: t-values to use for the smooth curve (useful when data has missing values or
               is "scattered"); if not provided, then `t̂ = t`; must be monotonically
               increasing
- `wls::{Vector,Symbol}`: weights to use with the least-squares fitting term; if set to
                          `:midpoint`, then midpoint-rule integration weights are used for
                          _both_ `wls` and `wr`
- `wr::Vector`: weights to use with the roughness term
- `d::Int = 2`: derivative used to calculate roughness; e.g., when `d = 2`, the 2nd
                derivative (i.e. the curvature) of the data is used to calculate roughness.

# Keyword Arguments
- `λ::{Number,Tuple} = 1.0`: regularization parameter; larger values result in a smoother
                             curve; the provided value is used directly when `alg = :fixed`;
                             otherwise it is used as an initial guess for the optimization
                             method, or as bounds if a 2-tuple is provided (TBD)
- `alg::Symbol = :gcv_svd`: algorithm for determining an optimal value for λ; the provided λ
                            value is used directly if `alg = :fixed`; otherwise `alg =
                            [:gcv_svd, :gcv_tr, :L_curve]` is passed to the
                            RegularizationTools solver

## Example Constructors
Smoothing using all arguments
```julia
A = RegularizationSmooth(u, t, t̂, wls, wr, d; λ=[1.0], alg=[:gcv_svd])
```
"""
function RegularizationSmooth(u::AbstractVector, t::AbstractVector, t̂::AbstractVector,
    wls::AbstractVector, wr::AbstractVector, d::Int = 2;
    λ::Real = 1.0, alg::Symbol = :gcv_svd)
    u, t = munge_data(u, t)
    M = _mapping_matrix(t̂, t)
    Wls½ = LA.diagm(sqrt.(wls))
    Wr½ = LA.diagm(sqrt.(wr))
    û, λ, Aitp = _reg_smooth_solve(u, t̂, d, M, Wls½, Wr½, λ, alg)
    RegularizationSmooth{true}(u, û, t, t̂, wls, wr, d, λ, alg, Aitp)
end
"""
Direct smoothing, no `t̂` or weights
```julia
A = RegularizationSmooth(u, t, d; λ=[1.0], alg=[:gcv_svd])
```
"""
function RegularizationSmooth(u::AbstractVector, t::AbstractVector, d::Int = 2;
    λ::Real = 1.0,
    alg::Symbol = :gcv_svd)
    u, t = munge_data(u, t)
    t̂ = t
    N = length(t)
    M = Array{Float64}(LA.I, N, N)
    Wls½ = Array{Float64}(LA.I, N, N)
    Wr½ = Array{Float64}(LA.I, N - d, N - d)
    û, λ, Aitp = _reg_smooth_solve(u, t̂, d, M, Wls½, Wr½, λ, alg)
    RegularizationSmooth{true}(u, û, t, t̂, LA.diag(Wls½), LA.diag(Wr½), d, λ, alg, Aitp)
end
"""
`t̂` provided, no weights
```julia
A = RegularizationSmooth(u, t, t̂, d; λ=[1.0], alg=[:gcv_svd])
```
"""
function RegularizationSmooth(u::AbstractVector, t::AbstractVector, t̂::AbstractVector,
    d::Int = 2; λ::Real = 1.0, alg::Symbol = :gcv_svd)
    u, t = munge_data(u, t)
    N, N̂ = length(t), length(t̂)
    M = _mapping_matrix(t̂, t)
    Wls½ = Array{Float64}(LA.I, N, N)
    Wr½ = Array{Float64}(LA.I, N̂ - d, N̂ - d)
    û, λ, Aitp = _reg_smooth_solve(u, t̂, d, M, Wls½, Wr½, λ, alg)
    RegularizationSmooth{true}(u, û, t, t̂, LA.diag(Wls½), LA.diag(Wr½), d, λ, alg, Aitp)
end
"""
`t̂` and `wls` provided
```julia
A = RegularizationSmooth(u, t, t̂, wls, d; λ=[1.0], alg=[:gcv_svd])
```
"""
function RegularizationSmooth(u::AbstractVector, t::AbstractVector, t̂::AbstractVector,
    wls::AbstractVector, d::Int = 2; λ::Real = 1.0,
    alg::Symbol = :gcv_svd)
    u, t = munge_data(u, t)
    N, N̂ = length(t), length(t̂)
    M = _mapping_matrix(t̂, t)
    Wls½ = LA.diagm(sqrt.(wls))
    Wr½ = Array{Float64}(LA.I, N̂ - d, N̂ - d)
    û, λ, Aitp = _reg_smooth_solve(u, t̂, d, M, Wls½, Wr½, λ, alg)
    RegularizationSmooth{true}(u, û, t, t̂, wls, LA.diag(Wr½), d, λ, alg, Aitp)
end
"""
`wls` provided, no `t̂`
```julia
A = RegularizationSmooth(u, t, nothing, wls,d; λ=[1.0], alg=[:gcv_svd])
```
"""
function RegularizationSmooth(u::AbstractVector, t::AbstractVector, t̂::Nothing,
    wls::AbstractVector, d::Int = 2; λ::Real = 1.0,
    alg::Symbol = :gcv_svd)
    u, t = munge_data(u, t)
    t̂ = t
    N = length(t)
    M = Array{Float64}(LA.I, N, N)
    Wls½ = LA.diagm(sqrt.(wls))
    Wr½ = Array{Float64}(LA.I, N - d, N - d)
    û, λ, Aitp = _reg_smooth_solve(u, t̂, d, M, Wls½, Wr½, λ, alg)
    RegularizationSmooth{true}(u, û, t, t̂, wls, LA.diag(Wr½), d, λ, alg, Aitp)
end
"""
`wls` and `wr` provided, no `t̂`
```julia
A = RegularizationSmooth(u, t, nothing, wls, wr, d; λ=[1.0], alg=[:gcv_svd])
```
"""
function RegularizationSmooth(u::AbstractVector, t::AbstractVector, t̂::Nothing,
    wls::AbstractVector, wr::AbstractVector, d::Int = 2;
    λ::Real = 1.0, alg::Symbol = :gcv_svd)
    u, t = munge_data(u, t)
    t̂ = t
    N = length(t)
    M = Array{Float64}(LA.I, N, N)
    Wls½ = LA.diagm(sqrt.(wls))
    Wr½ = LA.diagm(sqrt.(wr))
    û, λ, Aitp = _reg_smooth_solve(u, t̂, d, M, Wls½, Wr½, λ, alg)
    RegularizationSmooth{true}(u, û, t, t̂, wls, LA.diag(Wr½), d, λ, alg, Aitp)
end
"""
Keyword provided for `wls`, no `t̂`
```julia
A = RegularizationSmooth(u, t, nothing, :midpoint, d; λ=[1.0], alg=[:gcv_svd])
```
"""
function RegularizationSmooth(u::AbstractVector, t::AbstractVector, t̂::Nothing,
    wls::Symbol, d::Int = 2; λ::Real = 1.0, alg::Symbol = :gcv_svd)
    u, t = munge_data(u, t)
    t̂ = t
    N = length(t)
    M = Array{Float64}(LA.I, N, N)
    wls, wr = _weighting_by_kw(t, d, wls)
    Wls½ = LA.diagm(sqrt.(wls))
    Wr½ = LA.diagm(sqrt.(wr))
    û, λ, Aitp = _reg_smooth_solve(u, t̂, d, M, Wls½, Wr½, λ, alg)
    RegularizationSmooth{true}(u, û, t, t̂, LA.diag(Wls½), LA.diag(Wr½), d, λ, alg, Aitp)
end
# """ t̂ provided and keyword for wls  _TBD_ """
# function RegularizationSmooth(u::AbstractVector, t::AbstractVector, t̂::AbstractVector,
#                               wls::Symbol, d::Int=2; λ::Real=1.0, alg::Symbol=:gcv_svd)

""" Solve for the smoothed dependent variables and create spline interpolator """
function _reg_smooth_solve(u::AbstractVector, t̂::AbstractVector, d::Int, M::AbstractMatrix,
    Wls½::AbstractMatrix, Wr½::AbstractMatrix, λ::Real, alg::Symbol)
    λ = float(λ) # `float` expected by RT
    D = _derivative_matrix(t̂, d)
    Ψ = RT.setupRegularizationProblem(Wls½ * M, Wr½ * D)
    Wls½u = Wls½ * u

    if alg == :fixed
        b̄ = RT.to_standard_form(Ψ, Wls½u) # via b̄
        ū = RT.solve(Ψ, b̄, λ)
        û = RT.to_general_form(Ψ, Wls½u, ū)
    else
        # the provided λ (a scalar) is used as an initial guess; using bounds for Brent()
        # method is TBD, JJS 12/21/21
        result = RT.solve(Ψ, Wls½u; alg = alg, method = RT.NelderMead(), λ₀ = λ)
        û = result.x
        λ = result.λ
    end
    Aitp = CubicSpline(û, t̂)
    # It seems logical to use B-Spline of order d+1, but I am unsure if theory supports the
    # extra computational cost, JJS 12/25/21
    #Aitp = BSplineInterpolation(û,t̂,d+1,:ArcLen,:Average)
    return û, λ, Aitp
end

"""
Order d derivative matrix for the provided t vector
"""
function _derivative_matrix(t::AbstractVector, d::Int)
    N = length(t)
    if d == 0
        return Array{Float64}(LA.I, (N, N))
    end
    dt = t[(d + 1):end] - t[1:(end - d)]
    V = LA.diagm(1 ./ dt)
    Ddm1 = diff(_derivative_matrix(t, d - 1), dims = 1)
    D = d * V * Ddm1
    return D
end

"""Linear interpolation mapping matrix, which maps `û` to `u`."""
function _mapping_matrix(t̂::AbstractVector, t::AbstractVector)
    N = length(t)
    N̂ = length(t̂)
    # map the scattered points to the appropriate index of the smoothed points
    idx = searchsortedlast.(Ref(t̂), t)
    # allow for "extrapolation"; i.e. for t̂ extremum that are interior to t
    idx[idx .== 0] .+= 1
    idx[idx .== N̂] .+= -1
    # create the linear interpolation matrix
    m2 = @. (t - t̂[idx]) / (t̂[idx + 1] - t̂[idx])
    M = zeros(eltype(t), (N, N̂))
    for i in 1:N
        M[i, idx[i]] = 1 - m2[i]
        M[i, idx[i] + 1] = m2[i]
    end
    return M
end

""" Common-use weighting, currently only `:midpoint` for midpoint-rule integration """
function _weighting_by_kw(t::AbstractVector, d::Int, wls::Symbol)
    # `:midpoint` only for now, but plan to add functionality for `:relative` weighting
    N = length(t)
    if wls == :midpoint
        bmp = zeros(N)
        bmp[1] = -t[1] + t[2]
        for i in 2:(N - 1)
            bmp[i] = -t[i - 1] + t[i + 1]
        end
        bmp[N] = -t[N - 1] + t[N]
        # divide by 2 doesn't matter in the minimize step, but keeping for correctness if
        # used as a template elsewhere
        bmp = bmp / 2
        start = floor(Int, d / 2) + 1
        final = iseven(d) ? N - (start - 1) : N - start
        b̃mp = bmp[start:final]
        return bmp, b̃mp
    else
        throw("Unknown `$(wls)` keyword used for weighting, use `:midpoint`.")
    end
end

function DataInterpolations._interpolate(A::RegularizationSmooth{
        <:AbstractVector{<:Number},
    },
    t::Number)
    DataInterpolations._interpolate(A.Aitp, t)
end

end # module
