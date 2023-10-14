################################################################################
#                                 Type recipes                                 #
################################################################################

function to_plottable(A::AbstractInterpolation; plotdensity = 10_000, denseplot = true)
    t = sort(A.t)
    start = t[1]
    stop = t[end]
    if denseplot
        plott = collect(range(start, stop = stop, length = plotdensity))
    else
        plott = t
    end
    output = A.(plott)
    plott, output
end

@recipe function f(A::AbstractInterpolation; plotdensity = 10_000, denseplot = true)
    to_plottable(A; plotdensity = plotdensity, denseplot = denseplot)
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

@recipe function f(::Type{Val{:linear_interp}},
    x,
    y,
    z;
    plotdensity = 10_000,
    denseplot = true)
    seriestype := :path

    label --> "Linear fit"

    nx, ny = to_plottable(LinearInterpolation(y, x);
        plotdensity = plotdensity,
        denseplot = denseplot)

    x := nx
    y := ny
end

########################################
#       Quadratic Interpolation        #
########################################

@recipe function f(::Type{Val{:quadratic_interp}},
    x,
    y,
    z;
    plotdensity = 10_000,
    denseplot = true)
    seriestype := :path

    label --> "Quadratic fit"

    nx, ny = to_plottable(QuadraticInterpolation(T.(y),
            T.(x));
        plotdensity = plotdensity,
        denseplot = denseplot)
    x := nx
    y := ny
end

########################################
#           Quadratic Spline           #
########################################

@recipe function f(::Type{Val{:quadratic_spline}},
    x,
    y,
    z;
    plotdensity = 10_000,
    denseplot = true)
    seriestype := :path

    label --> "Quadratic Spline"

    T = promote_type(eltype(y), eltype(x))

    nx, ny = to_plottable(QuadraticSpline(T.(y),
            T.(x));
        plotdensity = plotdensity,
        denseplot = denseplot)

    x := nx
    y := ny
end

########################################
#        Lagrange Interpolation        #
########################################

@recipe function f(::Type{Val{:lagrange_interp}},
    x, y, z;
    n = length(x) - 1,
    plotdensity = 10_000,
    denseplot = true)
    seriestype := :path

    label --> "Lagrange Fit"

    T = promote_type(eltype(y), eltype(x))

    nx, ny = to_plottable(LagrangeInterpolation(T.(y),
            T.(x),
            n);
        plotdensity = plotdensity,
        denseplot = denseplot)

    x := nx
    y := ny
end

########################################
#             Cubic Spline             #
########################################

@recipe function f(::Type{Val{:cubic_spline}},
    x,
    y,
    z;
    plotdensity = 10_000,
    denseplot = true)
    seriestype := :path

    label --> "Cubic Spline"

    T = promote_type(eltype(y), eltype(x))

    nx, ny = to_plottable(CubicSpline(T.(y),
            T.(x));
        plotdensity = plotdensity,
        denseplot = denseplot)
    x := nx
    y := ny
end

@recipe function f(::Type{Val{:bspline_interp}},
    x, y, z;
    d = 5,
    pVec = :ArcLen,
    knotVec = :Average,
    plotdensity = length(x) * 6,
    denseplot = true)
    seriestype := :path

    label --> "B-Spline"

    @show x y eltype(x)

    # T = promote_type(eltype(y), eltype(x))

    nx, ny = to_plottable(BSplineInterpolation(T.(y),
            T.(x),
            d,
            pVec,
            knotVec);
        plotdensity = plotdensity,
        denseplot = denseplot)
    x := nx
    y := ny
end

########################################
#       B-spline (approximation)       #
########################################

@recipe function f(::Type{Val{:bspline_approx}},
    x, y, z;
    d = 5,
    h = length(x) - 1,
    pVec = :ArcLen,
    knotVec = :Average,
    plotdensity = length(x) * 6,
    denseplot = true)
    seriestype := :path

    label --> "B-Spline"

    T = promote_type(eltype(y), eltype(x))

    nx, ny = to_plottable(BSplineApprox(T.(y),
            T.(x),
            d,
            h,
            pVec,
            knotVec);
        plotdensity = plotdensity,
        denseplot = denseplot)
    x := nx
    y := ny
end

########################################
#          Akima interpolation          #
########################################

@recipe function f(::Type{Val{:akima}},
    x,
    y,
    z;
    plotdensity = length(x) * 6,
    denseplot = true)
    seriestype := :path

    label --> "Akima"

    T = promote_type(eltype(y), eltype(x))

    nx, ny = to_plottable(AkimaInterpolation(T.(y),
            T.(x));
        plotdensity = plotdensity,
        denseplot = denseplot)
    x := nx
    y := ny
end

################################################################################
#                                  Shorthands                                  #
################################################################################

"""
    akima(u, t)
    akima!(u, t)

Plot the Akima interpolation of the given data.
"""
@shorthands akima
