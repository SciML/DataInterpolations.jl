derivative(A, t) = derivative(A, t, firstindex(A.t)-1)[1]


function derivative(A::LinearInterpolation{<:AbstractVector}, t::Number, iguess)
  idx = searchsortedfirst(A.t, t, iguess)
  if A.t[idx] >= t
    idx -= 1
  end
  idx == 0 ? idx += 1 : nothing
  θ = 1 / (A.t[idx+1] - A.t[idx])
  (A.u[idx+1] - A.u[idx]) / (A.t[idx+1] - A.t[idx]), idx
end

function derivative(A::LinearInterpolation{<:AbstractMatrix}, t::Number, iguess)
  idx = searchsortedfirst(A.t, t, iguess)
  if A.t[idx] >= t
    idx -= 1
  end
  idx == 0 ? idx += 1 : nothing
  θ = 1 / (A.t[idx+1] - A.t[idx])
  (@views @. (A.u[:, idx+1] - A.u[:, idx]) / (A.t[idx+1] - A.t[idx])), idx
end

function derivative(A::QuadraticInterpolation{<:AbstractVector}, t::Number, iguess)
  i₀, i₁, i₂ = _quad_interp_indices(A, t, iguess)
  dl₀ = (2t - A.t[i₁] - A.t[i₂]) / ((A.t[i₀] - A.t[i₁]) * (A.t[i₀] - A.t[i₂]))
  dl₁ = (2t - A.t[i₀] - A.t[i₂]) / ((A.t[i₁] - A.t[i₀]) * (A.t[i₁] - A.t[i₂]))
  dl₂ = (2t - A.t[i₀] - A.t[i₁]) / ((A.t[i₂] - A.t[i₀]) * (A.t[i₂] - A.t[i₁]))
  A.u[i₀] * dl₀ + A.u[i₁] * dl₁ + A.u[i₂] * dl₂, i₀
end

function derivative(A::QuadraticInterpolation{<:AbstractMatrix}, t::Number, iguess)
  idx = searchsortedfirst(A.t, t, iguess)
  if A.t[idx] >= t
      idx -= 1
  end
  idx == 0 ? idx += 1 : nothing
  if idx == length(A.t) - 1
    i₀ = idx - 1; i₁ = idx; i₂ = i₁ + 1;
  else
    i₀ = idx; i₁ = i₀ + 1; i₂ = i₁ + 1;
  end
  dl₀ = (2t - A.t[i₁] - A.t[i₂]) / ((A.t[i₀] - A.t[i₁]) * (A.t[i₀] - A.t[i₂]))
  dl₁ = (2t - A.t[i₀] - A.t[i₂]) / ((A.t[i₁] - A.t[i₀]) * (A.t[i₁] - A.t[i₂]))
  dl₂ = (2t - A.t[i₀] - A.t[i₁]) / ((A.t[i₂] - A.t[i₀]) * (A.t[i₂] - A.t[i₁]))
  (@views @. A.u[:, i₀] * dl₀ + A.u[:, i₁] * dl₁ + A.u[:, i₂] * dl₂), idx
end

function derivative(A::LagrangeInterpolation{<:AbstractVector}, t::Number)
  idxs = findRequiredIdxs(A, t)
  if A.t[idxs[1]] == t
    return zero(A.u[idxs[1]])
  end
  G = zero(A.u[1]); F = zero(A.t[1])
  DG = zero(A.u[1]); DF = zero(A.t[1])
  tmp = G
  for i = 1:length(idxs)
    if isnan(A.bcache[idxs[i]])
      mult = one(A.t[1])
      for j = 1:(i - 1)
        mult *= (A.t[idxs[i]] - A.t[idxs[j]])
      end
      for j = (i+1):length(idxs)
        mult *= (A.t[idxs[i]] - A.t[idxs[j]])
      end
      A.bcache[idxs[i]] = mult
    else
      mult = A.bcache[idxs[i]]
    end
    wi = inv(mult)
    tti = t - A.t[idxs[i]]
    tmp = wi / (t - A.t[idxs[i]])
    g = tmp * A.u[idxs[i]]
    G += g
    DG -= g / (t - A.t[idxs[i]])
    F += tmp
    DF -= tmp / (t - A.t[idxs[i]])
  end
  (DG * F - G * DF) / (F ^ 2)
end

function derivative(A::LagrangeInterpolation{<:AbstractMatrix}, t::Number)
  idxs = findRequiredIdxs(A, t)
  if A.t[idxs[1]] == t
    return zero(A.u[:, idxs[1]])
  end
  G = zero(A.u[:, 1]); F = zero(A.t[1])
  DG = zero(A.u[:, 1]); DF = zero(A.t[1])
  tmp = G
  for i = 1:length(idxs)
    if isnan(A.bcache[idxs[i]])
      mult = one(A.t[1])
      for j = 1:(i - 1)
        mult *= (A.t[idxs[i]] - A.t[idxs[j]])
      end
      for j = (i+1):length(idxs)
        mult *= (A.t[idxs[i]] - A.t[idxs[j]])
      end
      A.bcache[idxs[i]] = mult
    else
      mult = A.bcache[idxs[i]]
    end
    wi = inv(mult)
    tti = t - A.t[idxs[i]]
    tmp = wi / (t - A.t[idxs[i]])
    g = tmp * A.u[:, idxs[i]]
    @. G += g
    @. DG -= g / (t - A.t[idxs[i]])
    F += tmp
    DF -= tmp / (t - A.t[idxs[i]])
  end
  @. (DG * F - G * DF) / (F ^ 2)
end

function derivative(A::AkimaInterpolation{<:AbstractVector}, t::Number, iguess)
  t < A.t[1] && return zero(A.u[1]), 1
  t > A.t[end] && return zero(A.u[end]), lastindex(t)
  i = searchsortedlast(A.t, t, iguess)
  j = min(i, length(A.c))  # for smooth derivative at A.t[end]
  wj = t - A.t[i]
  (@evalpoly wj A.b[i] 2A.c[j] 3A.d[j]), i
end

function derivative(A::ConstantInterpolation{<:AbstractVector}, t::Number)
  return isempty(searchsorted(A.t, t)) ? zero(A.u[1]) : eltype(A.u)(NaN)
end

function derivative(A::ConstantInterpolation{<:AbstractMatrix}, t::Number)
  return isempty(searchsorted(A.t, t)) ? zero(A.u[:, 1]) : eltype(A.u)(NaN) .* A.u[:, 1]
end

# QuadraticSpline Interpolation
function derivative(A::QuadraticSpline{<:AbstractVector}, t::Number, iguess)
  i = searchsortedfirst(A.t, t, iguess)
  i == 1 ? i += 1 : nothing
  σ = 1//2 * (A.z[i] - A.z[i - 1]) / (A.t[i] - A.t[i - 1])
  A.z[i-1] + 2σ * (t - A.t[i-1]), i
end

# CubicSpline Interpolation
function derivative(A::CubicSpline{<:AbstractVector}, t::Number, iguess)
  i = searchsortedfirst(A.t, t, iguess)
  isnothing(i) ? i = length(A.t) - 1 : i -= 1
  i == 0 ? i += 1 : nothing
  dI = -3A.z[i] * (A.t[i + 1] - t)^2 / (6A.h[i + 1]) + 3A.z[i + 1] * (t - A.t[i])^2 / (6A.h[i + 1])
  dC = A.u[i + 1] / A.h[i + 1] - A.z[i + 1] * A.h[i + 1] / 6
  dD = -(A.u[i] / A.h[i + 1] - A.z[i] * A.h[i + 1] / 6)
  dI + dC + dD, i
end

function derivative(A::BSplineInterpolation{<:AbstractVector{<:Number}}, t::Number, iguess)
  # change t into param [0 1]
  idx = searchsortedlast(A.t,t, iguess)
  idx == length(A.t) ? idx -= 1 : nothing
  n = length(A.t)
  scale = (A.p[idx+1] - A.p[idx]) / (A.t[idx+1] - A.t[idx])
  t_ = A.p[idx] + (t - A.t[idx]) * scale
  N = DataInterpolations.spline_coefficients(n, A.d-1, A.k, t_)
  ducum = zero(eltype(A.u))
  for i = 1:(n - 1)
    ducum += N[i + 1] * (A.c[i + 1] - A.c[i]) / (A.k[i + A.d + 1] - A.k[i + 1])
  end
  ducum * A.d * scale, idx
end

# BSpline Curve Approx
function derivative(A::BSplineApprox{<:AbstractVector{<:Number}}, t::Number, iguess)
  # change t into param [0 1]
  idx = searchsortedlast(A.t,t, iguess)
  idx == 0 ? idx += 1 : nothing
  scale = (A.p[idx+1] - A.p[idx]) / (A.t[idx+1] - A.t[idx])
  t_ = A.p[idx] + (t - A.t[idx]) * scale
  N = spline_coefficients(A.h, A.d-1, A.k, t_)
  ducum = zero(eltype(A.u))
  for i = 1:(A.h - 1)
    ducum += N[i + 1] * (A.c[i + 1] - A.c[i]) / (A.k[i + A.d + 1] - A.k[i + 1])
  end
  ducum * A.d * scale, idx
end
