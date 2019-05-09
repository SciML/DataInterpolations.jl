### Linear Interpolation
struct LinearInterpolation{uType,tType,FT,T} <: AbstractInterpolation{FT,T}
  u::uType
  t::tType
  LinearInterpolation{FT}(u,t) where FT = new{typeof(u),typeof(t),FT,eltype(u)}(u,t)
end

function LinearInterpolation(u,t)
  u, t = munge_data(u, t)
  LinearInterpolation{true}(u,t)
end

### Quadratic Interpolation
struct QuadraticInterpolation{uType,tType,FT,T} <: AbstractInterpolation{FT,T}
  u::uType
  t::tType
  QuadraticInterpolation{FT}(u,t) where FT = new{typeof(u),typeof(t),FT,eltype(u)}(u,t)
end

function QuadraticInterpolation(u,t)
  u, t = munge_data(u, t)
  QuadraticInterpolation{true}(u,t)
end

### Lagrange Interpolation
struct LagrangeInterpolation{uType,tType,FT,T} <: AbstractInterpolation{FT,T}
  u::uType
  t::tType
  n::Int
  LagrangeInterpolation{FT}(u,t,n) where FT = new{typeof(u),typeof(t),FT,eltype(u)}(u,t,n)
end

function LagrangeInterpolation(u,t,n=nothing)
  u, t = munge_data(u, t)
  if n == nothing
    n = length(t) - 1 # degree
  end
  if n != length(t) - 1
    error("Currently only n=length(t) - 1 is supported")
  end
  LagrangeInterpolation{true}(u,t,n)
end

### ZeroSpline Interpolation
struct ZeroSpline{uType,tType,dirType,FT,T} <: AbstractInterpolation{FT,T}
  u::uType
  t::tType
  dir::Symbol
  ZeroSpline{FT}(u,t,dir) where FT = new{typeof(u),typeof(t),typeof(dir),FT,eltype(u)}(u,t,dir)
end

function ZeroSpline(u,t;dir=:left)
  u, t = munge_data(u, t)
  ZeroSpline{true}(u,t,dir)
end

### QuadraticSpline Interpolation
struct QuadraticSpline{uType,tType,tAType,dType,zType,FT,T} <: AbstractInterpolation{FT,T}
  u::uType
  t::tType
  tA::tAType
  d::dType
  z::zType
  QuadraticSpline{FT}(u,t,tA,d,z) where FT = new{typeof(u),typeof(t),typeof(tA),
                                                  typeof(d),typeof(z),FT,eltype(u)}(u,t,tA,d,z)
end

function QuadraticSpline(u,t)
  u, t = munge_data(u, t)
  s = length(t)
  dl = ones(eltype(t),s-1)
  d = ones(eltype(t),s)
  du = zeros(eltype(t),s-1)
  tA = Tridiagonal(dl,d,du)
  d = zero(t)
  for i = 2:length(d)
    d[i] = 2//1 * (u[i] - u[i-1])/(t[i] - t[i-1])
  end
  z = tA\d
  QuadraticSpline{true}(u,t,tA,d,z)
end

# Cubic Spline Interpolation
struct CubicSpline{uType,tType,hType,zType,FT,T} <: AbstractInterpolation{FT,T}
  u::uType
  t::tType
  h::hType
  z::zType
  CubicSpline{FT}(u,t,h,z) where FT = new{typeof(u),typeof(t),typeof(h),typeof(z),FT,eltype(u)}(u,t,h,z)
end

function CubicSpline(u,t)
  u, t = munge_data(u, t)
  n = length(t) - 1
  h = vcat(0, diff(t), 0)
  dl = h[2:n+1]
  d = 2 .* (h[1:n+1] .+ h[2:n+2])
  du = h[2:n+1]
  tA = LinearAlgebra.Tridiagonal(dl,d,du)
  d = zero(t)
  for i = 2:n
    d[i] = 6(u[i+1] - u[i]) / h[i+1] - 6(u[i] - u[i-1]) / h[i]
  end
  z = tA\d
  CubicSpline{true}(u,t,h[1:n+1],z)
end

### BSpline Curve Interpolation
struct BSplineInterpolation{uType,tType,pType,kType,cType,FT,T} <: AbstractInterpolation{FT,T}
  u::uType
  t::tType
  d::Int    # degree
  p::pType  # params vector
  k::kType  # knot vector
  c::cType  # control points
  pVecType::Symbol
  knotVecType::Symbol
  BSplineInterpolation{FT}(u,t,d,p,k,c,pVecType,knotVecType) where FT =  new{typeof(u),typeof(t),typeof(p),typeof(k),typeof(c),FT,eltype(u)}(u,t,d,p,k,c,pVecType,knotVecType)
end

function BSplineInterpolation(u,t,d,pVecType,knotVecType)
  u, t = munge_data(u, t)
  n = length(t)
  s = zero(eltype(u))
  p = zero(t)
  k = zeros(eltype(t),n+d+1)
  l = zeros(eltype(u),n-1)
  p[1] = zero(eltype(t))
  p[end] = one(eltype(t))

  for i = 2:n
    s += √((t[i] - t[i-1])^2 + (u[i] - u[i-1])^2)
    l[i-1] = s
  end
  if pVecType == :Uniform
    for i = 2:(n-1)
      p[i] = p[1] + (i-1)*(p[end]-p[1])/(n-1)
    end
  elseif pVecType == :ArcLen
    for i = 2:(n-1)
      p[i] = p[1] + l[i-1]/s * (p[end]-p[1])
    end
  end

  lidx = 1
  ridx = length(k)
  while lidx <= (d+1) && ridx >= (length(k)-d)
    k[lidx] = p[1]
    k[ridx] = p[end]
    lidx += 1
    ridx -= 1
  end

  ps = zeros(eltype(t),n-2)
  s = zero(eltype(t))
  for i = 2:n-1
    s += p[i]
    ps[i-1] = s
  end

  if knotVecType == :Uniform
    # uniformly spaced knot vector
    # this method is not recommended because, if it is used with the chord length method for global interpolation,
    # the system of linear equations would be singular.
    for i = (d+2):n
      k[i] = k[1] + (i-d-1)//(n-d) * (k[end]-k[1])
    end
  elseif knotVecType == :Average
    # average spaced knot vector
    idx = 1
    if d+2 <= n
      k[d+2] = 1//d * ps[d]
    end
    for i = (d+3):n
      k[i] = 1//d * (ps[idx+d] - ps[idx])
      idx += 1
    end
  end
  # control points
  N = spline_coefficients(n,d,k,p)
  c = vec(N\u[:,:])
  BSplineInterpolation{true}(u,t,d,p,k,c,pVecType,knotVecType)
end

### BSpline Curve Approx
struct BSplineApprox{uType,tType,pType,kType,cType,FT,T} <: AbstractInterpolation{FT,T}
  u::uType
  t::tType
  d::Int    # degree
  h::Int    # number of control points (n => h >= d >= 1)
  p::pType  # params vector
  k::kType  # knot vector
  c::cType  # control points
  pVecType::Symbol
  knotVecType::Symbol
  BSplineApprox{FT}(u,t,d,h,p,k,c,pVecType,knotVecType) where FT =  new{typeof(u),typeof(t),typeof(p),typeof(k),typeof(c),FT,eltype(u)}(u,t,d,h,p,k,c,pVecType,knotVecType)
end

function BSplineApprox(u,t,d,h,pVecType,knotVecType)
  u, t = munge_data(u, t)
  n = length(t)
  s = zero(eltype(u))
  p = zero(t)
  k = zeros(eltype(t),h+d+1)
  l = zeros(eltype(u),n-1)
  p[1] = zero(eltype(t))
  p[end] = one(eltype(t))

  for i = 2:n
    s += √((t[i] - t[i-1])^2 + (u[i] - u[i-1])^2)
    l[i-1] = s
  end
  if pVecType == :Uniform
    for i = 2:(n-1)
      p[i] = p[1] + (i-1)*(p[end]-p[1])/(n-1)
    end
  elseif pVecType == :ArcLen
    for i = 2:(n-1)
      p[i] = p[1] + l[i-1]/s * (p[end]-p[1])
    end
  end

  lidx = 1
  ridx = length(k)
  while lidx <= (d+1) && ridx >= (length(k)-d)
    k[lidx] = p[1]
    k[ridx] = p[end]
    lidx += 1
    ridx -= 1
  end

  ps = zeros(eltype(t),n-2)
  s = zero(eltype(t))
  for i = 2:n-1
    s += p[i]
    ps[i-1] = s
  end

  if knotVecType == :Uniform
    # uniformly spaced knot vector
    # this method is not recommended because, if it is used with the chord length method for global interpolation,
    # the system of linear equations would be singular.
    for i = (d+2):h
      k[i] = k[1] + (i-d-1)//(h-d) * (k[end]-k[1])
    end
  elseif knotVecType == :Average
    # NOTE: verify that average method can be applied when size of k is less than size of p
    # average spaced knot vector
    idx = 1
    if d+2 <= h
      k[d+2] = 1//d * ps[d]
    end
    for i = (d+3):h
      k[i] = 1//d * (ps[idx+d] - ps[idx])
      idx += 1
    end
  end
  # control points
  c = zeros(eltype(u), h)
  c[1] = u[1]
  c[end] = u[end]
  q = zeros(eltype(u), n)
  N = zeros(eltype(t), n, h)
  for i = 1:n
    N[i, :] .= spline_coefficients(h,d,k,p[i])
  end
  for k = 2:n-1
    q[k] = u[k] - N[k, 1]*u[1] - N[k, h]*u[end]
  end
  Q = Matrix{eltype(u)}(undef, h-2, 1)
  for i = 2:h-1
    s = 0.0
    for k = 2:n-1
      s += N[k, i]*q[k]
    end
    Q[i-1] = s
  end
  N = N[2:end-1,2:h-1]
  M = transpose(N) * N
  P = M\Q
  c[2:end-1] .= vec(P)
  BSplineApprox{true}(u,t,d,h,p,k,c,pVecType,knotVecType)
end

### Loess
struct Loess{uType,tType,αType,xType,FT,T} <: AbstractInterpolation{FT,T}
  u::uType
  t::tType
  d::Int
  α::αType
  q::Int
  x::xType
  Loess{FT}(u,t,d,α,q,x) where FT = new{typeof(u),typeof(t),typeof(α),typeof(x),FT,eltype(u)}(u,t,d,α,q,x)
end

function Loess(u,t,d,α)
  u, t = munge_data(u, t)
  n = length(t)
  q = floor(Int,n*α)
  x = Matrix{eltype(t)}(undef,n,d+1)
  x[:,1] .= one(t[1])
  for i = 2:(d+1)
    x[:,i] = t .^ (i-1)
  end
  Loess{true}(u,t,d,α,q,x)
end

### GaussianProcesses
struct GPInterpolation{uType,tType,gpType,FT,T} <: AbstractInterpolation{FT,T}
  u::uType
  t::tType
  gp::gpType
  GPInterpolation{FT}(u,t,gp) where FT = new{typeof(u),typeof(t),typeof(gp),FT,eltype(u)}(u,t,gp)
end

function GPInterpolation(u,t,m,k,n=-2.0)
  u, t = munge_data(u, t)
  gp = GP(t,u,m,k,n)
  GPInterpolation{true}(u,t,gp)
end

### Curvefit
struct CurvefitCache{uType,tType,mType,p0Type,ubType,lbType,algType,pminType,FT,T} <: AbstractInterpolation{FT,T}
  u::uType
  t::tType
  m::mType        # model type
  p0::p0Type      # intial params
  ub::ubType      # upper bound of params
  lb::lbType      # lower bound of params
  alg::algType    # alg to optimize cost function
  pmin::pminType  # optimized params
  CurvefitCache{FT}(u,t,m,p0,ub,lb,alg,pmin) where FT = new{typeof(u),typeof(t),typeof(m),
                                      typeof(p0),typeof(ub),typeof(lb),
                                      typeof(alg),typeof(pmin),FT,eltype(u)}(u,t,m,p0,ub,lb,alg,pmin)
end

function Curvefit(u, t, model, p0, alg, box=false, lb=nothing, ub=nothing)
  u, t = munge_data(u, t)
  errfun(t,u,p) = sum(abs2.(u .- model(t,p)))
  if box == false
      mfit = optimize(p->errfun(t,u,p), p0, alg)
  else
      if lb == nothing || ub == nothing
          error("lower or upper bound should not be nothing")
      end
      od = OnceDifferentiable(p->errfun(t,u,p), p0, autodiff=:finite)
      mfit = optimize(od, lb, ub, p0, Fminbox(alg))
  end
  pmin = Optim.minimizer(mfit)
  CurvefitCache{true}(u,t,model,p0,ub,lb,alg,pmin)
end
