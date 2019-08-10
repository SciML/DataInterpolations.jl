################################################################################
#                                 Type recipes                                 #
################################################################################

function to_plottable(A::AbstractInterpolation; plotdensity = 10_000, denseplot = true)
  t = sort(A.t)
  start = t[1]; stop = t[end]
  if denseplot
    plott = collect(range(start,stop=stop,length=plotdensity))
  else
    plott = t
  end
  output = A.(plott)
  plott, output
end

function to_plottable(A::GPInterpolation; plotdensity = 10_000, denseplot = true)
  t = sort(A.t)
  start = t[1]; stop = t[end]
  if denseplot
    plott = collect(range(start,stop=stop,length=plotdensity))
  else
    plott = t
  end
  plott, A(plott)
end

@recipe function f(A::AbstractInterpolation; plotdensity = 10_000, denseplot = true)
    to_plottable(A; plotdensity = plotdensity, denseplot=denseplot)
end

@recipe function f(A::GPInterpolation; plotdensity = 10_000, denseplot = true)
    to_plottable(A; plotdensity = plotdensity, denseplot=denseplot)
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

@recipe function f(::Type{Val{:linear_interp}}, x, y, z; plotdensity = 10_000, denseplot = true)

    seriestype := :path

    label --> "Linear fit"

    nx, ny = to_plottable(LinearInterpolation(y, x); plotdensity = plotdensity, denseplot = denseplot)

    x := nx
    y := ny
end

########################################
#       Quadratic Interpolation        #
########################################

@recipe function f(::Type{Val{:quadratic_interp}}, x, y, z; plotdensity = 10_000, denseplot = true)

    seriestype := :path

    label --> "Quadratic fit"

    nx, ny = to_plottable(
        QuadraticInterpolation(
            T.(y),
            T.(x)
        );
        plotdensity = plotdensity,
        denseplot = denseplot
    )
    x := nx
    y := ny
end

########################################
#           Quadratic Spline           #
########################################

@recipe function f(::Type{Val{:quadratic_spline}}, x, y, z; plotdensity = 10_000, denseplot = true)

    seriestype := :path

    label --> "Quadratic Spline"

    T = promote_type(eltype(y), eltype(x))

    nx, ny = to_plottable(
        QuadraticSpline(
            T.(y),
            T.(x)
        );
        plotdensity = plotdensity,
        denseplot = denseplot
    )

    x := nx
    y := ny
end

########################################
#        Lagrange Interpolation        #
########################################

@recipe function f(::Type{Val{:lagrange_interp}}, x, y, z; n = length(x) - 1, plotdensity = 10_000, denseplot = true)

    seriestype := :path

    label --> "Lagrange Fit"

    T = promote_type(eltype(y), eltype(x))

    nx, ny = to_plottable(
        LagrangeInterpolation(
            T.(y),
            T.(x),
            n
        );
        plotdensity = plotdensity,
        denseplot = denseplot
    )

    x := nx
    y := ny
end

########################################
#             Cubic Spline             #
########################################

@recipe function f(::Type{Val{:cubic_spline}}, x, y, z; plotdensity = 10_000, denseplot = true)

    seriestype := :path

    label --> "Cubic Spline"

    T = promote_type(eltype(y), eltype(x))

    nx, ny = to_plottable(
        CubicSpline(
            T.(y),
            T.(x)
        );
        plotdensity = plotdensity,
        denseplot = denseplot
    )
    x := nx
    y := ny
end

########################################
#                Loess                 #
########################################

@recipe function f(::Type{Val{:loess}}, x, y, z; d = 5, α = 0.75, plotdensity = length(x) * 6, denseplot = true)

    seriestype := :path

    label --> "LOESS fit"

    nx, ny = to_plottable(Loess(y, x, d, α); plotdensity = plotdensity, denseplot = denseplot)

    x := nx
    y := ny
end

@recipe function f(::Type{Val{:bspline_interp}}, x, y; d = 5, pVec=:ArcLen, knotVec = :Average)

        seriestype := :path

        label --> "B-Spline"

        T = promote_type(eltype(y), eltype(x))

        nx, ny = to_plottable(
            BSplineInterpolation(
                T.(y),
                T.(x),
                d,
                pVec,
                knotVec
            );
            plotdensity = plotdensity,
            denseplot = denseplot
        )
        x := nx
        y := ny
end

@recipe function f(::Type{Val{:bspline_approx}}, x, y; d = 5, h = length(x)-1, pVec=:ArcLen, knotVec = :Average)

        seriestype := :path

        label --> "B-Spline"

        T = promote_type(eltype(y), eltype(x))

        nx, ny = to_plottable(
            BSplineApprox(
                T.(y),
                T.(x),
                d,
                h,
                pVec,
                knotVec
            );
            plotdensity = plotdensity,
            denseplot = denseplot
        )
        x := nx
        y := ny
end
