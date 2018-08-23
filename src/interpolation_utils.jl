function findRequiredIdxs(A, t)
  idx = findfirst(x->x>=t,A.t)  # assuming A.t is already sorted
  if idx == 1 || idx == 2
    x₀ = A.t[1]; x₁ = A.t[2]; x₂ = A.t[3]
  elseif idx == length(A.t)
    x₀ = A.t[end-2]; x₁ = A.t[end-1]; x₂ = A.t[end]
  else
    x = A.t[idx-1]; y = A.t[idx]
    if(abs(t-A.t[idx-2])<=abs(t-A.t[idx+1]))
      x₀ = A.t[idx-2]; x₁ = x; x₂ = y
    else
      x₀ = x; x₁ = y; x₂ = A.t[idx+1]
    end
  end
  i₀, i₁, i₂ = 0, 0, 0
  for i = 1:length(A.t)
    if A.t[i] == x₀
        i₀ = i
    elseif A.t[i] == x₁
        i₁ = i
    elseif A.t[i] == x₂
        i₂ = i
    end
  end
  return i₀, i₂, i₁
end
