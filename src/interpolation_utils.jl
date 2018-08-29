function findRequiredIdxs(A::QuadraticInterpolation, t)
  idxs = sortperm(A.t,by=x->abs(t-x))
  return idxs[1], idxs[2], idxs[3]
end

function findRequiredIdxs(A::LagrangeInterpolation, t)
  idxs = sortperm(A.t,by=x->abs(t-x))
  idxs[1:(A.n+1)]
end
