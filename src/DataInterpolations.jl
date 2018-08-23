module DataInterpolations

include("interpolation_utils.jl")

### Interface Functionality

abstract type AbstractInterpolation{FT,T} <: AbstractVector{T} end

Base.size(A::AbstractInterpolation) = size(A.u)
Base.size(A::AbstractInterpolation{true}) = length(A.u) .+ size(A.t)
Base.getindex(A::AbstractInterpolation,i) = A.u[i]
Base.getindex(A::AbstractInterpolation{true},i) = i<=length(A.u) ? A.u[i] : A.t[i-length(A.u)]
Base.setindex!(A::AbstractInterpolation,x,i) = A.u[i] = x
Base.setindex!(A::AbstractInterpolation{true},x,i) =
                                   i <= length(A.u) ? (A.u[i] = x) : (A.t[i-length(A.u)] = x)

### Linear Interpolation

struct LinearInterpolation{uType,tType,FT,T} <: AbstractInterpolation{FT,T}
  u::uType
  t::tType
  LinearInterpolation{FT}(u,t) where FT = new{typeof(u),typeof(t),FT,eltype(u)}(u,t)
end
LinearInterpolation(u,t) = LinearInterpolation{true}(u,t)

struct QuadraticInterpolation{uType,tType,FT,T} <: AbstractInterpolation{FT,T}
  u::uType
  t::tType
  QuadraticInterpolation{FT}(u,t) where FT = new{typeof(u),typeof(t),FT,eltype(u)}(u,t)
end
QuadraticInterpolation(u,t) = QuadraticInterpolation{true}(u,t)

function (A::LinearInterpolation{<:AbstractVector{<:Number}})(t::Number)
  idx = findfirst(x->x>=t,A.t)-1
  θ = (t - A.t[idx])/ (A.t[idx+1] - A.t[idx])
  (1-θ)*A.u[idx] + θ*A.u[idx+1]
end

function (A::LinearInterpolation{<:AbstractMatrix{<:Number}})(t::Number)
  idx = findfirst(x->x>=t,A.t)-1
  θ = (t - A.t[idx])/ (A.t[idx+1] - A.t[idx])
  (1-θ)*A.u[:,idx] + θ*A.u[:,idx+1]
end

function (A::QuadraticInterpolation{<:AbstractVector{<:Number}})(t::Number)
  i₀, i₁, i₂ = findRequiredIdxs(A,t)
  l₀ = ((t-A.t[i₁])*(t-A.t[i₂]))/((A.t[i₀]-A.t[i₁])*(A.t[i₀]-A.t[i₂]))
  l₁ = ((t-A.t[i₀])*(t-A.t[i₂]))/((A.t[i₁]-A.t[i₀])*(A.t[i₁]-A.t[i₂]))
  l₂ = ((t-A.t[i₀])*(t-A.t[i₁]))/((A.t[i₂]-A.t[i₀])*(A.t[i₂]-A.t[i₁]))
  A.u[i₀]*l₀ + A.u[i₁]*l₁ + A.u[i₂]*l₂
end

export LinearInterpolation, QuadraticInterpolation

end # module
