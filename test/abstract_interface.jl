using DataInterpolations, Test
using FindFirstFunctions

struct GenericTestInterpolation <: DataInterpolations.AbstractInterpolation{Float64}
    u::Vector{Float64}
    t::Vector{Float64}
    I::Vector{Float64}
    iguesser
    extrapolation_left::DataInterpolations.ExtrapolationType.T
    extrapolation_right::DataInterpolations.ExtrapolationType.T
    kind::FindFirstFunctions.StrategyKind
    t_props::FindFirstFunctions.SearchProperties
    cache_parameters::Bool
end

# The contract test supplies only the developer hooks. All assertions below exercise the
# public generic evaluation, derivative, integral, and shape interfaces.
function DataInterpolations._interpolate(
        ::GenericTestInterpolation, t::Number, ::Any
    )
    return 2t
end

function DataInterpolations._derivative(
        ::GenericTestInterpolation, ::Number, ::Any
    )
    return 2.0
end

function DataInterpolations._integral(
        ::GenericTestInterpolation, ::Integer, t1::Number, t2::Number
    )
    return 2 * (t2 - t1)
end

@testset "AbstractInterpolation contract" begin
    t = [0.0, 1.0, 2.0]
    t_props = FindFirstFunctions.SearchProperties(t)
    A = GenericTestInterpolation(
        [0.0, 2.0, 4.0], t, Float64[], FindFirstFunctions.Guesser(t),
        DataInterpolations.ExtrapolationType.None,
        DataInterpolations.ExtrapolationType.None,
        FindFirstFunctions.Auto(t, t_props).kind, t_props, false,
    )

    @test A isa DataInterpolations.AbstractInterpolation
    @test A(0.5) == 1.0
    @test A([0.5, 1.5]) == [1.0, 3.0]

    out = zeros(2)
    @test A(out, [0.5, 1.5]) === out
    @test out == [1.0, 3.0]

    @test DataInterpolations.derivative(A, 0.5) == 2.0
    @test DataInterpolations.derivative(A, 0.5, 2) == 0.0
    @test DataInterpolations.integral(A, 1.5) == 3.0
    @test DataInterpolations.integral(A, 0.5, 1.5) == 2.0

    @test DataInterpolations.output_dim(A) == 0
    @test DataInterpolations.output_size(A) == ()
end
