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
