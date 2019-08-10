################################################################################
#                                 Type recipes                                 #
################################################################################

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

################################################################################
#                                Series recipes                                #
################################################################################


############################################################
#                      Interpolations                      #
############################################################

########################################
#         Linear Interpolation         #
########################################

@recipe function f(::Type{Val{:linear_interp}}, x, y, z)

    label --> "Linear fit"

    LinearInterpolation(y, x)
end

########################################
#       Quadratic Interpolation        #
########################################

@recipe function f(::Type{Val{:quadratic_interp}}, x, y, z) # TODO are these good values?

    label --> "Quadratic fit"

    QuadraticInterpolation(y, x)
end

########################################
#           Quadratic Spline           #
########################################

@recipe function f(::Type{Val{:quadratic_spline}}, x, y, z)

    label --> "Quadratic Spline"

    QuadraticSpline(y, x)
end

########################################
#        Lagrange Interpolation        #
########################################

@recipe function f(::Type{Val{:lagrange_interp}}, x, y, z; n = length(x) - 1)

    label --> "Lagrange Fit"

    LagrangeInterpolation(y, x, n)
end

########################################
#             Cubic Spline             #
########################################

@recipe function f(::Type{Val{:cubic_spline}}, x, y, z)

    label --> "Cubic Spline"

    CubicSpline(y, x)
end

########################################
#                Loess                 #
########################################

@recipe function f(::Type{Val{:loess}}, x, y, z; d = 5, α = 0.75) # TODO are these good values?

    label --> "LOESS fit"

    Loess(y, x, d, α)
end

# @recipe function f(::Type{Val{:bspline_interp}}, x, y; pVec=:)
