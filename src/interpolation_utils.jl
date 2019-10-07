function findRequiredIdxs(A::LagrangeInterpolation, t)
  idxs = sortperm(A.t,by=x->abs(t-x))
  idxs[1:(A.n+1)]
end

function spline_coefficients(n, d, k, u::Number)
  N = zeros(n)
  if u == k[1]
    N[1] = one(u)
  elseif u == k[end]
    N[end] = one(u)
  else
    i = findfirst(x->x>u,k) - 1
    N[i] = one(u)
    for deg = 1:d
      N[i-deg] = (k[i+1]-u)/(k[i+1]-k[i-deg+1]) * N[i-deg+1]
      for j = (i-deg+1):(i-1)
        N[j] = (u-k[j])/(k[j+deg]-k[j]) * N[j] + (k[j+deg+1]-u)/(k[j+deg+1]-k[j+1]) * N[j+1]
      end
      N[i] = (u-k[i])/(k[i+deg]-k[i]) * N[i]
    end
  end
  N
end

function spline_coefficients(n, d, k, u::AbstractVector)
  N = zeros(eltype(u),n,n)
  for i = 1:n
    N[i,:] .= spline_coefficients(n,d,k,u[i])
  end
  N
end

function weibull_fun(x,p)
  ans = copy(x)
  for i = 1:length(x)
    if x[i] >= 0
      ans[i] = p[1] * (1 - exp(-1 * (x[i]/p[2])^p[3]))
    else
      ans[i] = 0.0
    end
  end
  ans
end

# helper functions to get t_max and c_max (OneCompartmentPK)
function OneCompartmentPK_tmax(A)
  p = A.c_f.param
  2.303 * log(10.0,p[1]/p[4]) / (p[1] - p[4])
end

function OneCompartmentPK_cmax(A)
  p = A.c_f.param
  t = OneCompartmentPK_tmax(A)
  p[1] * p[2] * p[3] * (exp(-p[4]*t) - exp(-p[1]*t)) / (p[5] * (p[1] - p[4]))
end

# helper function for data manipulation
function munge_data(u::AbstractVector, t::AbstractVector)
  Tu = Base.nonmissingtype(eltype(u))
  Tt = Base.nonmissingtype(eltype(t))
  newu = Tu[]
  newt = Tt[]
  @assert length(t) == length(u)
  @inbounds for i in eachindex(t)
    ui = u[i]
    ti = t[i]
    if !ismissing(ui) && !ismissing(ti)
      push!(newu, ui)
      push!(newt, ti)
    end
  end
  return newu, newt
end

function munge_data(U::StridedMatrix, t::AbstractVector)
  TU = Base.nonmissingtype(eltype(U))
  Tt = Base.nonmissingtype(eltype(t))
  newUs = [TU[] for i in 1:size(U, 1)]
  newt  = Tt[]
  @assert length(t) == size(U,2)
  @inbounds for j in eachindex(t)
    tj = t[j]
    if ismissing(tj) || any(ismissing, view(U, :, j))
      continue
    end

    push!(newt, tj)
    for i in 1:length(newUs)
      push!(newUs[i], U[i,j])
    end
  end

  return vcat(adjoint.(newUs)...), newt
end
