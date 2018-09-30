function findRequiredIdxs(A::QuadraticInterpolation, t)
  idxs = sortperm(A.t,by=x->abs(t-x))
  return idxs[1], idxs[2], idxs[3]
end

function findRequiredIdxs(A::LagrangeInterpolation, t)
  idxs = sortperm(A.t,by=x->abs(t-x))
  idxs[1:(A.n+1)]
end

function compute_splines(A::BSpline, t)
  n = length(A.t)
  d,k = A.d,A.k
  B = zero(A.t)
  if t == k[1]
    B[1] = one(A.t[1])
  elseif t == k[end]
    B[end] = one(A.t[1])
  else
    idx = findfirst(x->x>t,A.k) - 1
    B[idx] = one(A.t[1])
    for deg = 1:d
      B[idx-deg+1] = (k[idx+2]-t)/(k[idx+2]-k[idx-deg+2]) * B[idx-deg+2]
      for i = idx-deg+2:idx
        B[i] = (t-k[i])/(k[i+deg]-k[i]) * B[i] + (k[i+deg+1]-t)/(k[i+deg+1]-k[i+1]) * B[i+1]
      end
      B[idx+1] = (t-k[idx+1])/(k[idx+1+deg]-k[idx+1]) * B[idx+1]
    end
  end
  B    
end
