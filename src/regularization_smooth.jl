# TODO:
# - midpoint rule
# x scattered/interpolation
# x GCV
# x L-curve
# - bounds on λ
# - initial guess for λ? will require mods to RegularizationTools
# - scaled λ? will need to work out equivalency with λ² formulation, and resolve with
#   derivative rather than difference matrix
# - standard deviation
# - relative weights? (not sure if worth implementing)
# - arbitrary weighting


const LA = LinearAlgebra

# might move this to RegularizationTools.jl? JJS 12/6/21
"""
Order d derivative matrix for the provided t vector
"""
function _derivative_matrix(t::AbstractVector, d::Int)
    N = length(t)
    if d == 0
        return Array{Float64}(LA.I, (N, N))
    end
    dt = t[d+1:end] - t[1:end-d]
    V = LA.diagm(1 ./ dt)
    # could use this symbol:  Dᵈ⁻¹ (D\^d\^-\^1) (but a lot of typing...)
    Ddm1 = diff(_derivative_matrix(t, d-1), dims=1)
    D = d * V*Ddm1
    return D
end

### Regularization data smoothing and interpolation
struct RegularizationSmooth{uType,tType,FT,T} <: AbstractInterpolation{FT,T}
    u::uType
    û::uType
    t::tType
    t̂::tType
    d::Int              # derivative degree used to calculate the roughness
    λ::AbstractFloat    # regularization parameter
    alg::Symbol  # how to determine λ: `:fixed`, `:gcv_svd`, `:gcv_tr`, `L_curve`
    Aitp::AbstractInterpolation{FT,T}
    RegularizationSmooth{FT}(u,û,t,t̂,d,λ,alg,Aitp) where FT =
        new{typeof(u),typeof(t),FT,eltype(u)}(u,û,t,t̂,d,λ,alg,Aitp)
end

function RegularizationSmooth(u, t, d=2; λ=1.0, t̂=nothing, alg=:gcv_svd)
    # TBD:  make use of `munge_data`
    u, t = munge_data(u, t)
    λ = float(λ) # in case an integer is given for λ
    N = length(t)
    if isnothing(t̂)
        t̂ = t
        # identity matrix for A (RegularizationTools `A` matrix, not `A` object of
        # DataInterpolations)
        A = Array{Float64}(LA.I, (N, N))
    else
        A = _mapping_matrix(t̂, t) # mapping matrix for A
    end

    # RT.Γ function, used internally there, does not not account for spacing of x 
    #Ψ = RT.setupRegularizationProblem(A, d)
    D = _derivative_matrix(t̂,d)
    Ψ = RT.setupRegularizationProblem(A, D)
    
    if alg == :fixed
        b̄ = RT.to_standard_form(Ψ, u) # via b̄
        ū = RT.solve(Ψ, b̄, λ)
        û = RT.to_general_form(Ψ, u, ū)
    else
        # note that the provided λ is currently ignored
        result = RT.solve(Ψ, u, alg=alg)
        û = result.x
        λ = result.λ
    end
    # I'm not sure about mathematical theory for performing a posteriori interpolation via
    # Regularization smoothing -- applying Cubic-spline interpolation seems reasonable to
    # me, so setting that up here to use after object-creation; could allow users the option
    # of a posteriori interpolation or not, or even a choice for the interpolation,
    # depending on what is needed. B-spline interpolation of order d+1 also seems to be a
    # reasonable choice, but it is not allowing extrapolation and doesn't show much benefit
    # over cubic-spline. JJS 12/3/21
    Aitp = CubicSpline(û,t̂)
    #Aitp = BSplineInterpolation(û,t̂,d+1,:ArcLen,:Average)
    RegularizationSmooth{true}(u,û,t,t̂,d,λ,alg,Aitp)
end

"""Linear interpolation mapping matrix, which maps `û` to `u`."""
function _mapping_matrix(t̂::AbstractVector,t::AbstractVector)
    N = length(t)
    N̂ = length(t̂)
    # map the scattered points to the appropriate index of the smoothed points
    idx = searchsortedlast.(Ref(t̂),t)
    # allow for "extrapolation"; i.e. for t̂ extremum that are interior to t
    idx[idx.==0] .+= 1
    idx[idx.==N̂] .+= -1 
    # create the linear interpolation matrix
    m2 = @. (t - t̂[idx])/(t̂[idx+1] - t̂[idx])
    M = zeros((N,N̂))
    for i in 1:N
        M[i,idx[i]] = 1 - m2[i] 
        M[i,idx[i]+1] = m2[i]
    end
    return M
end

_interpolate(A::RegularizationSmooth{<:AbstractVector{<:Number}}, t::Number) =
    _interpolate(A.Aitp, t)
