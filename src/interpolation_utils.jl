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

# helper function for data manipulation
function munge_data(u::AbstractVector, t::AbstractVector)
  newu = DataFrames.dropmissing(u)
  newt = DataFrames.dropmissing(t)
  @assert length(t) == length(u)
  return newu, newt
end

function munge_data(U::StridedMatrix, t::AbstractVector)
  tmpU = [U[:,j] for j un length(t)]
  newUs = DataFrames.dropmissing(tmpU)
  newt  = DataFrames.dropmissing(t)
  @assert length(t) == size(U,2)
  return hcat(newUs...), newt
end
