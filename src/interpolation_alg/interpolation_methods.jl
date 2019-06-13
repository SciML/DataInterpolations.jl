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
  idx = findfirst(x->x>=t,A.t)-1
  idx == 0 ? idx += 1 : nothing
  if idx == length(A.t) - 1
    i₀ = idx - 1; i₁ = idx; i₂ = i₁ + 1;
  else
    i₀ = idx; i₁ = i₀ + 1; i₂ = i₁ + 1;
  end
  l₀ = ((t-A.t[i₁])*(t-A.t[i₂]))/((A.t[i₀]-A.t[i₁])*(A.t[i₀]-A.t[i₂]))
  l₁ = ((t-A.t[i₀])*(t-A.t[i₂]))/((A.t[i₁]-A.t[i₀])*(A.t[i₁]-A.t[i₂]))
  l₂ = ((t-A.t[i₀])*(t-A.t[i₁]))/((A.t[i₂]-A.t[i₀])*(A.t[i₂]-A.t[i₁]))
  A.u[i₀]*l₀ + A.u[i₁]*l₁ + A.u[i₂]*l₂
end

function (A::QuadraticInterpolation{<:AbstractMatrix{<:Number}})(t::Number)
  idx = findfirst(x->x>=t,A.t)-1
  idx == 0 ? idx += 1 : nothing
  if idx == length(A.t) - 1
    i₀ = idx - 1; i₁ = idx; i₂ = i₁ + 1;
  else
    i₀ = idx; i₁ = i₀ + 1; i₂ = i₁ + 1;
  end
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

# ZeroSpline Interpolation
function (A::ZeroSpline{<:AbstractVector{<:Number}})(t::Number)
  i = findfirst(x->x>=t,A.t)
  i == 1 ? (i += 1) : nothing
  if A.dir == :left
    _i = i == nothing ? length(A.t)+1 : (i == 1 ? 2 : i)
    return A.u[_i-1]
  else
    _i = i == nothing ? length(A.t) : i
    return A.u[_i]
  end
end
 function (A::ZeroSpline{<:AbstractMatrix{<:Number}})(t::Number)
  i = findfirst(x->x>=t,A.t)
  i == 1 ? (i += 1) : nothing
  if A.dir == :left
    _i = i == nothing ? length(x)+1 : (i == 1 ? 2 : i)
    return A.u[:,i-1]
  else
    _i = i == nothing ? length(x) : i
    return A.u[:,i]
  end
end

# QuadraticSpline Interpolation
function (A::QuadraticSpline{<:AbstractVector{<:Number}})(t::Number)
  i = findfirst(x->x>=t,A.t)
  i == 1 ? i += 1 : nothing
  Cᵢ = A.u[i-1]
  σ = 1//2 * (A.z[i] - A.z[i-1])/(A.t[i] - A.t[i-1])
  A.z[i-1] * (t - A.t[i-1]) + σ * (t - A.t[i-1])^2 + Cᵢ
end

# CubicSpline Interpolation
function (A::CubicSpline{<:AbstractVector{<:Number}})(t::Number)
  i = findfirst(x->x>=t,A.t)
  i == nothing ? i = length(A.t) - 1 : i -= 1
  i == 0 ? i += 1 : nothing
  I = A.z[i] * (A.t[i+1] - t)^3 / (6A.h[i+1]) + A.z[i+1] * (t - A.t[i])^3 / (6A.h[i+1])
  C = (A.u[i+1]/A.h[i+1] - A.z[i+1]*A.h[i+1]/6)*(t - A.t[i])
  D = (A.u[i]/A.h[i+1] - A.z[i]*A.h[i+1]/6)*(A.t[i+1] - t)
  I + C + D
end

# BSpline Curve Interpolation
function (A::BSplineInterpolation{<:AbstractVector{<:Number}})(t::Number)
  # change t into param [0 1]
  idx = findfirst(x->x>=t,A.t) - 1
  idx == 0 ? idx += 1 : nothing
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
function (A::BSplineApprox{<:AbstractVector{<:Number}})(t::Number)
  # change t into param [0 1]
  idx = findfirst(x->x>=t,A.t) - 1
  idx == 0 ? idx += 1 : nothing
  t = A.p[idx] + (t - A.t[idx])/(A.t[idx+1] - A.t[idx]) * (A.p[idx+1] - A.p[idx])
  n = length(A.t)
  N = spline_coefficients(A.h,A.d,A.k,t)
  ucum = zero(eltype(A.u))
  for i = 1:A.h
    ucum += N[i] * A.c[i]
  end
  ucum
end

# Loess
function (A::Loess{<:AbstractVector{<:Number}})(t::Number)
  tmp = sort(abs.(A.t .- t))
  w = abs.(A.t .- t) ./ tmp[A.q]
  for i = 1:length(A.t)
    if w[i] <= one(A.t[1])
      w[i] = (1 - (w[i] ^ 3)) ^ 3
    else
      w[i] = zero(A.t[1])
    end
  end
  w = Diagonal(w)
  b = inv(transpose(A.x) * w * A.x) * transpose(A.x) * w * A.u
  u = zero(t[1])
  for (idx,v) in enumerate(b)
    u += v*(t^(idx-1))
  end
  u
end

# GaussianProcess
function (A::GPInterpolation{<:AbstractVector{<:Number}})(t::AbstractVector{<:Number})
  rand(A.gp,t,1)
end

function (A::GPInterpolation{<:AbstractVector{<:Number}})(t::Number)
  rand(A.gp,[t],1)[1]
end

# Curvefit
function (A::CurvefitCache{<:AbstractVector{<:Number}})(t::Union{AbstractVector{<:Number},Number})
  A.m(t,A.pmin)
end
