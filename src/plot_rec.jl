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

@recipe function f(A::AbstractInterpolation; plotdensity = 10_000, denseplot = true,
        label_interp = string(nameof(typeof(A))), label_data = "Data points")
    seriescolor = get(plotattributes, :seriescolor, :blue)
    @series begin
        seriestype := :path
        seriescolor --> seriescolor
        label --> label_interp
        to_plottable(A; plotdensity = plotdensity, denseplot = denseplot)
    end
    @series begin
        seriestype := :scatter
        seriescolor --> seriescolor
        label --> label_data
        A.t, A.u
    end
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
        extrapolate = false,
        safetycopy = false,
        plotdensity = 10_000,
        denseplot = true,
        label_interp = "LinearInterpolation",
        label_data = "Data points")
    T = promote_type(eltype(y), eltype(x))
    nx, ny = to_plottable(LinearInterpolation(T.(y), T.(x); extrapolate, safetycopy);
        plotdensity = plotdensity,
        denseplot = denseplot)
    seriescolor = get(plotattributes, :seriescolor, :blue)
    @series begin
        seriestype := :path
        seriescolor --> seriescolor
        label --> label_interp
        x := nx
        y := ny
    end
    @series begin
        seriestype := :scatter
        seriescolor --> seriescolor
        label --> label_data
        x := x
        y := y
    end
end

########################################
#       Quadratic Interpolation        #
########################################

@recipe function f(::Type{Val{:quadratic_interp}},
        x,
        y,
        z;
        mode = :Forward,
        extrapolate = false,
        safetycopy = false,
        plotdensity = 10_000,
        denseplot = true,
        label_interp = "QuadraticInterpolation",
        label_data = "Data points")
    T = promote_type(eltype(y), eltype(x))
    nx, ny = to_plottable(
        QuadraticInterpolation(T.(y),
            T.(x), mode; extrapolate, safetycopy);
        plotdensity = plotdensity,
        denseplot = denseplot)
    seriescolor = get(plotattributes, :seriescolor, :blue)
    @series begin
        seriestype := :path
        seriescolor --> seriescolor
        label --> label_interp
        x := nx
        y := ny
    end
    @series begin
        seriestype := :scatter
        seriescolor --> seriescolor
        label --> label_data
        x := x
        y := y
    end
end

########################################
#        Lagrange Interpolation        #
########################################

@recipe function f(::Type{Val{:lagrange_interp}},
        x, y, z;
        n = length(x) - 1,
        extrapolate = false,
        safetycopy = false,
        plotdensity = 10_000,
        denseplot = true,
        label_interp = "LagrangeInterpolation",
        label_data = "Data points")
    T = promote_type(eltype(y), eltype(x))
    nx, ny = to_plottable(LagrangeInterpolation(T.(y),
            T.(x),
            n; extrapolate, safetycopy);
        plotdensity = plotdensity,
        denseplot = denseplot)
    seriescolor = get(plotattributes, :seriescolor, :blue)
    @series begin
        seriestype := :path
        seriescolor --> seriescolor
        label --> label_interp
        x := nx
        y := ny
    end
    @series begin
        seriestype := :scatter
        seriescolor --> seriescolor
        label --> label_data
        x := x
        y := y
    end
end

########################################
#           Quadratic Spline           #
########################################

@recipe function f(::Type{Val{:quadratic_spline}},
        x,
        y,
        z;
        extrapolate = false,
        safetycopy = false,
        plotdensity = 10_000,
        denseplot = true,
        label_interp = "QuadraticSpline",
        label_data = "Data points")
    T = promote_type(eltype(y), eltype(x))
    nx, ny = to_plottable(QuadraticSpline(T.(y),
            T.(x); extrapolate, safetycopy);
        plotdensity = plotdensity,
        denseplot = denseplot)
    seriescolor = get(plotattributes, :seriescolor, :blue)
    @series begin
        seriestype := :path
        seriescolor --> seriescolor
        label --> label_interp
        x := nx
        y := ny
    end
    @series begin
        seriestype := :scatter
        seriescolor --> seriescolor
        label --> label_data
        x := x
        y := y
    end
end

########################################
#             Cubic Spline             #
########################################

@recipe function f(::Type{Val{:cubic_spline}},
        x,
        y,
        z;
        extrapolate = false,
        safetycopy = false,
        plotdensity = 10_000,
        denseplot = true,
        label_interp = "CubicSpline",
        label_data = "Data points")
    T = promote_type(eltype(y), eltype(x))
    nx, ny = to_plottable(CubicSpline(T.(y),
            T.(x); extrapolate, safetycopy);
        plotdensity = plotdensity,
        denseplot = denseplot)
    seriescolor = get(plotattributes, :seriescolor, :blue)
    @series begin
        seriestype := :path
        seriescolor --> seriescolor
        label --> label_interp
        x := nx
        y := ny
    end
    @series begin
        seriestype := :scatter
        seriescolor --> seriescolor
        label --> label_data
        x := x
        y := y
    end
end

########################################
#          Akima interpolation          #
########################################

@recipe function f(::Type{Val{:akima_interp}},
        x,
        y,
        z;
        extrapolate = false,
        safetycopy = false,
        plotdensity = 10_000,
        denseplot = true,
        label_interp = "AkimaInterpolation",
        label_data = "Data points")
    T = promote_type(eltype(y), eltype(x))
    nx, ny = to_plottable(AkimaInterpolation(T.(y),
            T.(x); extrapolate, safetycopy);
        plotdensity = plotdensity,
        denseplot = denseplot)
    seriescolor = get(plotattributes, :seriescolor, :blue)
    @series begin
        seriestype := :path
        seriescolor --> seriescolor
        label --> label_interp
        x := nx
        y := ny
    end
    @series begin
        seriestype := :scatter
        seriescolor --> seriescolor
        label --> label_data
        x := x
        y := y
    end
end

########################################
#       B-spline Interpolation       #
########################################

@recipe function f(::Type{Val{:bspline_interp}},
        x, y, z;
        d = 5,
        pVecType = :ArcLen,
        knotVecType = :Average,
        extrapolate = false,
        safetycopy = false,
        plotdensity = 10_000,
        denseplot = true,
        label_interp = "BSplineInterpolation",
        label_data = "Data points")
    T = promote_type(eltype(y), eltype(x))
    nx, ny = to_plottable(
        BSplineInterpolation(T.(y),
            T.(x),
            d,
            pVecType,
            knotVecType; extrapolate, safetycopy);
        plotdensity = plotdensity,
        denseplot = denseplot)
    seriescolor = get(plotattributes, :seriescolor, :blue)
    @series begin
        seriestype := :path
        seriescolor --> seriescolor
        label --> label_interp
        x := nx
        y := ny
    end
    @series begin
        seriestype := :scatter
        seriescolor --> seriescolor
        label --> label_data
        x := x
        y := y
    end
end

########################################
#       B-spline (approximation)       #
########################################

@recipe function f(::Type{Val{:bspline_approx}},
        x, y, z;
        d = 5,
        h = length(x) - 1,
        pVecType = :ArcLen,
        knotVecType = :Average,
        extrapolate = false,
        safetycopy = false,
        plotdensity = 10_000,
        denseplot = true,
        label_interp = "BSplineApprox",
        label_data = "Data points")
    T = promote_type(eltype(y), eltype(x))
    nx, ny = to_plottable(
        BSplineApprox(T.(y),
            T.(x),
            d,
            h,
            pVecType,
            knotVecType; extrapolate, safetycopy);
        plotdensity = plotdensity,
        denseplot = denseplot)
    seriescolor = get(plotattributes, :seriescolor, :blue)
    @series begin
        seriestype := :path
        seriescolor --> seriescolor
        label --> label_interp
        x := nx
        y := ny
    end
    @series begin
        seriestype := :scatter
        seriescolor --> seriescolor
        label --> label_data
        x := x
        y := y
    end
end

########################################
#          Cubic Hermite Spline          #
########################################

@recipe function f(::Type{Val{:cubic_hermite_spline}},
        x,
        y,
        z;
        du = nothing,
        extrapolate = false,
        safetycopy = false,
        plotdensity = 10_000,
        denseplot = true,
        label_interp = "CubicHermiteSpline",
        label_data = "Data points")
    isnothing(du) && error("Provide `du` as a keyword argument.")
    T = promote_type(eltype(y), eltype(x))
    nx, ny = to_plottable(
        CubicHermiteSpline(T.(du), T.(y),
            T.(x); extrapolate, safetycopy);
        plotdensity = plotdensity,
        denseplot = denseplot)
    seriescolor = get(plotattributes, :seriescolor, :blue)
    @series begin
        seriestype := :path
        seriescolor --> seriescolor
        label --> label_interp
        x := nx
        y := ny
    end
    @series begin
        seriestype := :scatter
        seriescolor --> seriescolor
        label --> label_data
        x := x
        y := y
    end
end

########################################
#        PCHIP Interpolation           #
########################################

@recipe function f(::Type{Val{:pchip_interp}},
        x,
        y,
        z;
        extrapolate = false,
        safetycopy = false,
        plotdensity = 10_000,
        denseplot = true,
        label_interp = "PCHIPInterpolation",
        label_data = "Data points")
    T = promote_type(eltype(y), eltype(x))
    nx, ny = to_plottable(PCHIPInterpolation(T.(y),
            T.(x); extrapolate, safetycopy);
        plotdensity = plotdensity,
        denseplot = denseplot)
    seriescolor = get(plotattributes, :seriescolor, :blue)
    @series begin
        seriestype := :path
        seriescolor --> seriescolor
        label --> label_interp
        x := nx
        y := ny
    end
    @series begin
        seriestype := :scatter
        seriescolor --> seriescolor
        label --> label_data
        x := x
        y := y
    end
end

########################################
#          Quintic Hermite Spline          #
########################################

@recipe function f(::Type{Val{:quintic_hermite_spline}},
        x,
        y,
        z;
        du = nothing,
        ddu = nothing,
        extrapolate = false,
        safetycopy = false,
        plotdensity = 10_000,
        denseplot = true,
        label_interp = "QuinticHermiteSpline",
        label_data = "Data points")
    (isnothing(du) || isnothing(ddu)) &&
        error("Provide `du` and `ddu` as keyword arguments.")
    T = promote_type(eltype(y), eltype(x))
    nx, ny = to_plottable(
        QuinticHermiteSpline(T.(ddu), T.(du), T.(y),
            T.(x); extrapolate, safetycopy);
        plotdensity = plotdensity,
        denseplot = denseplot)
    seriescolor = get(plotattributes, :seriescolor, :blue)
    @series begin
        seriestype := :path
        seriescolor --> seriescolor
        label --> label_interp
        x := nx
        y := ny
    end
    @series begin
        seriestype := :scatter
        seriescolor --> seriescolor
        label --> label_data
        x := x
        y := y
    end
end
