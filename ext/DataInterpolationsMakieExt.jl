module DataInterpolationsMakieExt

using DataInterpolations
using DataInterpolations: AbstractInterpolation
using Makie

# Define the default type of plot that you want
Makie.plottype(::AbstractInterpolation) = Makie.ScatterLines

# Define the attributes that you want to use
function Makie.used_attributes(::Makie.PointBased, ::AbstractInterpolation)
    (:plotdensity, :denseplot)
end
function Makie.used_attributes(::Type{<:Makie.ScatterLines}, ::AbstractInterpolation)
    (:plotdensity, :denseplot)
end

# Define the conversion of the data to the plot
function Makie.convert_arguments(
        ::Makie.PointBased,
        A::AbstractInterpolation;
        plotdensity = 10_000,
        denseplot = true
)
    DataInterpolations.to_plottable(A; plotdensity = plotdensity, denseplot = denseplot)
end

# Define the conversion of the data to the plot for the ScatterLines type
# Note that this is a bit of a hack
# and should actually be handled by a plot! method,
# except that doesn't work anymore or does it?
function Makie.convert_arguments(
        ::Type{<:Makie.ScatterLines},
        A::AbstractInterpolation;
        plotdensity = 10_000,
        denseplot = true
)
    densex,
    densey = convert_arguments(
        Makie.PointBased(), A; plotdensity = plotdensity, denseplot = denseplot)
    return [
        Makie.SpecApi.Lines(densex, densey),
        Makie.SpecApi.Scatter(A.t, A.u)
    ]
end

end # module
