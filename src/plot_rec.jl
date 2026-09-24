################################################################################
#                                 Type recipes                                 #
################################################################################

# `A.(t)` returns a plain `Vector` for scalar-valued interpolations (e.g. `Number` or
# `String` data) — nothing to reshape. For vector/array-valued interpolations it returns
# one array per time point (`Vector{<:AbstractArray}`), regardless of whether `A.u` itself
# is stored as a `Matrix`, a `Vector` of `Vector`s, or a higher-rank `Array` — but Plots
# needs a single `(nt, ndim)` matrix: an x vector of length `nt` paired with a y matrix
# with `nt` rows, one column per output series. Flatten each output point (so array-valued
# points still reduce to a flat set of series) and stack into that orientation.
_plot_reshape(output::AbstractVector) = output
function _plot_reshape(output::AbstractVector{<:AbstractArray})
    return permutedims(reduce(hcat, vec.(output)))
end

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
    return plott, _plot_reshape(output)
end

@recipe function f(
        A::AbstractInterpolation; plotdensity = 10_000, denseplot = true,
        label_interp = string(nameof(typeof(A))), label_data = "Data points"
    )
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
        # `get_data` (from show.jl) reshapes `A.u` — whatever its native container — into
        # the same `(npoints, ndim)` orientation Plots needs. Using it here (rather than
        # re-evaluating `A.(A.t)`) matters for approximating types like `BSplineApprox`,
        # whose curve doesn't pass through the raw data: re-evaluating would silently plot
        # the smoothed curve's values as if they were the data points.
        A.t, get_data(A.u)
    end
end
