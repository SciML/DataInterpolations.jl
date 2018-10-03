# Linear Interpolation
function (A::LinearInterpolation{<:AbstractVector{<:Number}})(t::Number)
  idx = findfirst(x->x>=t,A.t)-1
  idx == 0 ? idx += 1 : nothing
  θ = (t - A.t[idx])/ (A.t[idx+1] - A.t[idx])
  (1-θ)*A.u[idx] + θ*A.u[idx+1]
end

function (A::LinearInterpolation{<:AbstractMatrix{<:Number}})(t::Number)
  idx = findfirst(x->x>=t,A.t)-1
  idx == 0 ? idx += 1 : nothing
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

# Lagrange Interpolation
function (A::LagrangeInterpolation{<:AbstractVector{<:Number}})(t::Number)
  idxs = findRequiredIdxs(A,t)
  if A.t[idxs[1]] == t
    return A.u[idxs[1]]
  end
  N = zero(A.u[1]); D = zero(A.t[1]); tmp = N
  for i = 1:length(idxs)
    mult = one(A.t[1])
    for j = 1:(i-1)
      mult *= (A.t[idxs[i]] - A.t[idxs[j]])
    end
    for j = (i+1):length(idxs)
      mult *= (A.t[idxs[i]] - A.t[idxs[j]])
    end
    tmp = inv((t - A.t[idxs[i]]) * mult)
    D += tmp
    N += (tmp * A.u[idxs[i]])
  end
  N/D
end

function (A::LagrangeInterpolation{<:AbstractMatrix{<:Number}})(t::Number)
  idxs = findRequiredIdxs(A,t)
  if A.t[idxs[1]] == t
    return A.u[:,idxs[1]]
  end
  N = zero(A.u[:,1]); D = zero(A.t[1]); tmp = D
  for i = 1:length(idxs)
    mult = one(A.t[1])
    for j = 1:(i-1)
      mult *= (A.t[idxs[i]] - A.t[idxs[j]])
    end
    for j = (i+1):length(idxs)
      mult *= (A.t[idxs[i]] - A.t[idxs[j]])
    end
    tmp = inv((t - A.t[idxs[i]]) * mult)
    D += tmp
    @. N += (tmp * A.u[:,idxs[i]])
  end
  N/D
end

# QuadraticSpline Interpolation
function (A::QuadraticSpline{<:AbstractVector{<:Number}})(t::Number)
  i = findfirst(x->x>=t,A.t)
  i == 1 ? i += 1 : nothing
  Cᵢ = A.u[i-1]
  σ = 1//2 * (A.z[i] - A.z[i-1])/(A.t[i] - A.t[i-1])
  A.z[i-1] * (t - A.t[i-1]) + σ * (t - A.t[i-1])^2 + Cᵢ
end

# BSpline Interpolation
function (A::BSpline{<:AbstractVector{<:Number}})(t::Number)
  # change t into param [0 1]
  t = (t-A.t[1])/(A.t[end]-A.t[1])
  B = compute_splines(A, t)
  ucum = zero(eltype(A.u))
  for i = 1:length(A.t)
    ucum += B[i] * A.u[i]
  end
  ucum
end
