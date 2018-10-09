### Linear Interpolation
struct LinearInterpolation{uType,tType,FT,T} <: AbstractInterpolation{FT,T}
  u::uType
  t::tType
  LinearInterpolation{FT}(u,t) where FT = new{typeof(u),typeof(t),FT,eltype(u)}(u,t)
end
LinearInterpolation(u,t) = LinearInterpolation{true}(u,t)

### Quadratic Interpolation
struct QuadraticInterpolation{uType,tType,FT,T} <: AbstractInterpolation{FT,T}
  u::uType
  t::tType
  QuadraticInterpolation{FT}(u,t) where FT = new{typeof(u),typeof(t),FT,eltype(u)}(u,t)
end
QuadraticInterpolation(u,t) = QuadraticInterpolation{true}(u,t)

### Lagrange Interpolation
struct LagrangeInterpolation{uType,tType,FT,T} <: AbstractInterpolation{FT,T}
  u::uType
  t::tType
  n::Int
  LagrangeInterpolation{FT}(u,t,n) where FT = new{typeof(u),typeof(t),FT,eltype(u)}(u,t,n)
end
LagrangeInterpolation(u,t,n) = LagrangeInterpolation{true}(u,t,n)

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
  n = length(t)
  q = floor(Int,n*α)
  w = Matrix{eltype(t)}(I,n,n)
  x = Matrix{eltype(t)}(undef,n,d+1)
  x[:,1] .= one(t[1])
  for i = 2:(d+1)
    x[:,i] = u .^ (i-1)
  end
  Loess{true}(u,t,d,α,q,x)
end
