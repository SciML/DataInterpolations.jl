# Linear Interpolation
function _interpolate(A::LinearInterpolation{<:AbstractVector}, t::Number)
  idx = max(1, min(searchsortedlast(A.t, t), length(A.t) - 1))
  θ = (t - A.t[idx])/(A.t[idx + 1] - A.t[idx])
  return (1 - θ)*A.u[idx] + θ*A.u[idx+1]
end

function _interpolate(A::LinearInterpolation{<:AbstractMatrix}, t::Number)
  idx = max(1, min(searchsortedlast(A.t, t), length(A.t) - 1))
  θ = (t - A.t[idx])/(A.t[idx + 1] - A.t[idx])
  return (1 - θ)*A.u[:,idx] + θ*A.u[:,idx + 1]
end

# Quadratic Interpolation
function _interpolate(A::QuadraticInterpolation{<:AbstractVector}, t::Number)
  idx = max(1, min(searchsortedlast(A.t, t), length(A.t) - 2))
  i₀, i₁, i₂ = idx, idx + 1, idx + 2
  l₀ = ((t - A.t[i₁])*(t - A.t[i₂]))/((A.t[i₀] - A.t[i₁])*(A.t[i₀] - A.t[i₂]))
  l₁ = ((t - A.t[i₀])*(t - A.t[i₂]))/((A.t[i₁] - A.t[i₀])*(A.t[i₁] - A.t[i₂]))
  l₂ = ((t - A.t[i₀])*(t - A.t[i₁]))/((A.t[i₂] - A.t[i₀])*(A.t[i₂] - A.t[i₁]))
  return A.u[i₀]*l₀ + A.u[i₁]*l₁ + A.u[i₂]*l₂
end

function _interpolate(A::QuadraticInterpolation{<:AbstractMatrix}, t::Number)
  idx = max(1, min(searchsortedlast(A.t, t), length(A.t) - 2))
  i₀, i₁, i₂ = idx, idx + 1, idx + 2
  l₀ = ((t - A.t[i₁])*(t - A.t[i₂]))/((A.t[i₀] - A.t[i₁])*(A.t[i₀] - A.t[i₂]))
  l₁ = ((t - A.t[i₀])*(t - A.t[i₂]))/((A.t[i₁] - A.t[i₀])*(A.t[i₁] - A.t[i₂]))
  l₂ = ((t - A.t[i₀])*(t - A.t[i₁]))/((A.t[i₂] - A.t[i₀])*(A.t[i₂] - A.t[i₁]))
  return A.u[:,i₀]*l₀ + A.u[:,i₁]*l₁ + A.u[:,i₂]*l₂
end

# Lagrange Interpolation
function _interpolate(A::LagrangeInterpolation{<:AbstractVector}, t::Number)
  idxs = findRequiredIdxs(A,t)
  if A.t[idxs[1]] == t
    return A.u[idxs[1]]
  end
  N = zero(A.u[1]); D = zero(A.t[1]); tmp = N
  for i = 1:length(idxs)
    if isnan(A.bcache[idxs[i]])
      mult = one(A.t[1])
      for j = 1:(i-1)
        mult *= (A.t[idxs[i]] - A.t[idxs[j]])
      end
      for j = (i+1):length(idxs)
        mult *= (A.t[idxs[i]] - A.t[idxs[j]])
      end
      A.bcache[idxs[i]] = mult
    else
      mult = A.bcache[idxs[i]]
    end
    tmp = inv((t - A.t[idxs[i]]) * mult)
    D += tmp
    N += (tmp * A.u[idxs[i]])
  end
  N/D
end

function _interpolate(A::LagrangeInterpolation{<:AbstractMatrix}, t::Number)
  idxs = findRequiredIdxs(A,t)
  if A.t[idxs[1]] == t
    return A.u[:,idxs[1]]
  end
  N = zero(A.u[:,1]); D = zero(A.t[1]); tmp = D
  for i = 1:length(idxs)
    if isnan(A.bcache[idxs[i]])
      mult = one(A.t[1])
      for j = 1:(i-1)
        mult *= (A.t[idxs[i]] - A.t[idxs[j]])
      end
      for j = (i+1):length(idxs)
        mult *= (A.t[idxs[i]] - A.t[idxs[j]])
      end
      A.bcache[idxs[i]] = mult
    else
      mult = A.bcache[idxs[i]]
    end
    tmp = inv((t - A.t[idxs[i]]) * mult)
    D += tmp
    @. N += (tmp * A.u[:,idxs[i]])
  end
  N/D
end

function _interpolate(A::AkimaInterpolation{<:AbstractVector}, t::Number)
  i = searchsortedlast(A.t, t)
  i == 0 && return A.u[1]
  i == length(A.t) && return A.u[end]
  wj = t - A.t[i]
  @evalpoly wj A.u[i] A.b[i] A.c[i] A.d[i]
end

# ConstantInterpolation Interpolation
function _interpolate(A::ConstantInterpolation{<:AbstractVector}, t::Number)
  if A.dir === :left
    # :left means that value to the left is used for interpolation
    i = searchsortedlast(A.t, t)
    return A.u[max(1, i)]
  else
    # :right means that value to the right is used for interpolation
    i = searchsortedfirst(A.t, t)
    return A.u[min(length(A.t), i)]
  end
end

function _interpolate(A::ConstantInterpolation{<:AbstractMatrix}, t::Number)
  if A.dir === :left
    # :left means that value to the left is used for interpolation
    i = searchsortedlast(A.t, t)
    return A.u[:, max(1, i)]
  else
    # :right means that value to the right is used for interpolation
    i = searchsortedfirst(A.t, t)
    return A.u[:, min(length(A.t), i)]
  end
end

# QuadraticSpline Interpolation
function _interpolate(A::QuadraticSpline{<:AbstractVector{<:Number}}, t::Number)
  i = findfirst(x->x>=t,A.t)
  i == 1 ? i += 1 : nothing
  Cᵢ = A.u[i-1]
  σ = 1//2 * (A.z[i] - A.z[i-1])/(A.t[i] - A.t[i-1])
  A.z[i-1] * (t - A.t[i-1]) + σ * (t - A.t[i-1])^2 + Cᵢ
end

# CubicSpline Interpolation
function _interpolate(A::CubicSpline{<:AbstractVector{<:Number}}, t::Number)
  i = findfirst(x->x>=t,A.t)
  i == nothing ? i = length(A.t) - 1 : i -= 1
  i == 0 ? i += 1 : nothing
  I = A.z[i] * (A.t[i+1] - t)^3 / (6A.h[i+1]) + A.z[i+1] * (t - A.t[i])^3 / (6A.h[i+1])
  C = (A.u[i+1]/A.h[i+1] - A.z[i+1]*A.h[i+1]/6)*(t - A.t[i])
  D = (A.u[i]/A.h[i+1] - A.z[i]*A.h[i+1]/6)*(A.t[i+1] - t)
  I + C + D
end

# BSpline Curve Interpolation
function _interpolate(A::BSplineInterpolation{<:AbstractVector{<:Number}}, t::Number)
  # change t into param [0 1]
  idx = searchsortedlast(A.t,t)
  idx == length(A.t) ? idx -= 1 : nothing
  t = A.p[idx] + (t - A.t[idx])/(A.t[idx+1] - A.t[idx]) * (A.p[idx+1] - A.p[idx])
  n = length(A.t)
  N = spline_coefficients(n,A.d,A.k,t)
  ucum = zero(eltype(A.u))
  for i = 1:n
    ucum += N[i] * A.c[i]
  end
  ucum
end

# BSpline Curve Approx
function _interpolate(A::BSplineApprox{<:AbstractVector{<:Number}}, t::Number)
  # change t into param [0 1]
  idx = searchsortedlast(A.t,t)
  idx == length(A.t) ? idx -= 1 : nothing
  t = A.p[idx] + (t - A.t[idx])/(A.t[idx+1] - A.t[idx]) * (A.p[idx+1] - A.p[idx])
  n = length(A.t)
  N = spline_coefficients(A.h,A.d,A.k,t)
  ucum = zero(eltype(A.u))
  for i = 1:A.h
    ucum += N[i] * A.c[i]
  end
  ucum
end

# Curvefit
function _interpolate(A::CurvefitCache{<:AbstractVector{<:Number}}, t::Union{AbstractVector{<:Number},Number})
  A.m(t,A.pmin)
end
