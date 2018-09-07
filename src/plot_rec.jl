@recipe function f(A::AbstractInterpolation; plotdensity = 10_000, denseplot = true)
  start = A.t[1]; stop = A.t[end]
  if denseplot
    plott = collect(range(start,stop=stop,length=plotdensity))
  else
    plott = A.t
  end
  output = A.(plott)
  if !(eltype(output) <: AbstractVector)
    output = [[x] for x in output]
  end
  DiffEqArray(output, plott)
end
