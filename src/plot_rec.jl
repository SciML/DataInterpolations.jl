@recipe function f(A::AbstractInterpolation; plotdensity = 10_000, denseplot = true)
  t = sort(A.t)
  start = t[1]; stop = t[end]
  if denseplot
    plott = collect(range(start,stop=stop,length=plotdensity))
  else
    plott = t
  end
  output = A.(plott)
  if !(eltype(output) <: AbstractVector)
    output = [[x] for x in output]
  end
  DiffEqArray(output, plott)
end

@recipe function f(A::GPInterpolation; plotdensity = 10_000, denseplot = true)
  t = sort(A.t)
  start = t[1]; stop = t[end]
  if denseplot
    plott = collect(range(start,stop=stop,length=plotdensity))
  else
    plott = t
  end
  plott, A(plott)
end
