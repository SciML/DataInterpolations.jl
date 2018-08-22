function findRequiredIdxs(A, t)
  tmp = sort(A.t)
  idx = findfirst(x->x>=t,tmp)
  if idx == 1
    x₀ = tmp[1]; x₁ = tmp[2]; x₂ = tmp[3]
  elseif idx == length(tmp)
    x₀ = tmp[end-2]; x₁ = tmp[end-1]; x₂ = tmp[end]
  else
    x = tmp[idx-1]; y = tmp[idx]
    if(abs(t-tmp[idx-2])<=abs(t-tmp[idx+1]))
      x₀ = tmp[idx-2]; x₁ = x; x₂ = y
    else
      x₀ = x; x₁ = y; x₂ = tmp[idx+1]
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
