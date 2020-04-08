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
  newUs = []
  newt  = Tt[]
  @assert length(t) == size(U,2)
  @inbounds for (j, tj) in enumerate(t)

    vUj = view(U, :, j)
    if ismissing(tj) || any(ismissing, vUj)
      continue
    end

    push!(newt, tj)
    push!(newUs, vUj)
  end

  return hcat(newUs...), newt
end
