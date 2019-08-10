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
function munge_data(u::AbstractVector, t)
  df = DataFrame(u = u, t = t)
  dropmissing!(df, disallowmissing=true)
  df.u, df.t
end

function munge_data(u::AbstractMatrix, t)
  df = convert(DataFrame, u')
  df.t = t
  dropmissing!(df, disallowmissing=true)
  t = df.t
  select!(df, Not(:t))
  convert(Matrix, df)', t
end
