# Linear Interpolation
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

# Quadratic Interpolation
function (A::QuadraticInterpolation{<:AbstractVector{<:Number}})(t::Number)
  i₀, i₁, i₂ = findRequiredIdxs(A,t)
  l₀ = ((t-A.t[i₁])*(t-A.t[i₂]))/((A.t[i₀]-A.t[i₁])*(A.t[i₀]-A.t[i₂]))
  l₁ = ((t-A.t[i₀])*(t-A.t[i₂]))/((A.t[i₁]-A.t[i₀])*(A.t[i₁]-A.t[i₂]))
  l₂ = ((t-A.t[i₀])*(t-A.t[i₁]))/((A.t[i₂]-A.t[i₀])*(A.t[i₂]-A.t[i₁]))
  A.u[i₀]*l₀ + A.u[i₁]*l₁ + A.u[i₂]*l₂
end

function (A::QuadraticInterpolation{<:AbstractMatrix{<:Number}})(t::Number)
  i₀, i₁, i₂ = findRequiredIdxs(A,t)
  l₀ = ((t-A.t[i₁])*(t-A.t[i₂]))/((A.t[i₀]-A.t[i₁])*(A.t[i₀]-A.t[i₂]))
  l₁ = ((t-A.t[i₀])*(t-A.t[i₂]))/((A.t[i₁]-A.t[i₀])*(A.t[i₁]-A.t[i₂]))
  l₂ = ((t-A.t[i₀])*(t-A.t[i₁]))/((A.t[i₂]-A.t[i₀])*(A.t[i₂]-A.t[i₁]))
  A.u[:,i₀]*l₀ + A.u[:,i₁]*l₁ + A.u[:,i₂]*l₂
end
