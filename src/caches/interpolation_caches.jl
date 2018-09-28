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

### BSpline Interpolation
struct BSpline{uType,tType,pType,kType,FT,T} <: AbstractInterpolation{FT,T}
  u::uType
  t::tType
  d::Int    # degree
  p::pType  # params vector
  k::kType  # knot vector
  BSpline{FT}(u,t) where FT =  new{typeof(u),typeof(t),typeof(p),typeof(k),FT,eltype(u)}(u,t,d)
end

function BSpline(u,t,d)
  n = length(t)
  s = zero(eltype(u))
  p = zero(t)
  l = zeros(eltype(u),n-1)
  for i = 2:n
    s += sqrt((t[i] - t[i-1])^2 + (u[i] - u[i-1])^2)
    l[i-1] = s
  end
  a = p[1] = pdomain[1]; b = p[end] = pdomain[2]
  for i = 2:(n-1)
    p[i] = a + l[i-1]/s * (b-a)
  end
  lk = n + d + 1
  k = zeros(eltype(t),lk)
  for i = (n+1):lk
    k[i] = one(t[1])
  end
  if knotVec == :Uniform
    # uniformly spaced knot vector
    for i = (d+1):n
      k[i] = (i-d)//(n-d)
    end
  elseif knotVec == :Average
    # average spaced knot vector
    for i = (d+1):n
      k[i] = 1//d * (ts[d] - ts[i-1])
    end
  end
  BSpline{true}(u,t,d,p,k)
end
