using DataInterpolations, BenchmarkTools
using StableRNGs

const SUITE = BenchmarkGroup()
const rng = StableRNG(123)

# Smooth 1D data
t = collect(range(0.0, 10.0, length = 200))
u = sin.(t) .+ 0.5 .* cos.(3 .* t)
t_eval = collect(range(0.05, 9.95, length = 1000))

# =============================================================================
# Construction
# =============================================================================

SUITE["construct"] = BenchmarkGroup()

SUITE["construct"]["LinearInterpolation"] = @benchmarkable LinearInterpolation($u, $t)
SUITE["construct"]["QuadraticInterpolation"] = @benchmarkable QuadraticInterpolation(
    $u, $t
)
SUITE["construct"]["CubicSpline"] = @benchmarkable CubicSpline($u, $t)
SUITE["construct"]["AkimaInterpolation"] = @benchmarkable AkimaInterpolation($u, $t)
SUITE["construct"]["PCHIPInterpolation"] = @benchmarkable PCHIPInterpolation($u, $t)
SUITE["construct"]["ConstantInterpolation"] = @benchmarkable ConstantInterpolation($u, $t)
SUITE["construct"]["LagrangeInterpolation"] = @benchmarkable LagrangeInterpolation(
    $(u[1:20]), $(t[1:20])
)

# =============================================================================
# Evaluation
# =============================================================================

SUITE["eval"] = BenchmarkGroup()

lin = LinearInterpolation(u, t)
quad = QuadraticInterpolation(u, t)
cubic = CubicSpline(u, t)
akima = AkimaInterpolation(u, t)
const_interp = ConstantInterpolation(u, t)

SUITE["eval"]["linear_scalar"] = @benchmarkable $lin(5.4321)
SUITE["eval"]["linear_vector"] = @benchmarkable $lin($t_eval)
SUITE["eval"]["quadratic_vector"] = @benchmarkable $quad($t_eval)
SUITE["eval"]["cubic_vector"] = @benchmarkable $cubic($t_eval)
SUITE["eval"]["akima_vector"] = @benchmarkable $akima($t_eval)
SUITE["eval"]["constant_vector"] = @benchmarkable $const_interp($t_eval)

function sum_scalar_queries(interpolation, queries)
    result = zero(eltype(queries))
    for query in queries
        result += interpolation(query)
    end
    return result
end

for n in (33, 513)
    knots = collect(range(0.0, 8.0; length = n))
    sites = knots .+ (1 .- cos.(knots)) ./ 1000
    spline = BSplineInterpolation(sin.(knots) .+ knots ./ 16, sites, 3, :Average)
    queries = collect(range(0.1, 7.9; length = 1024))
    out = similar(queries)
    SUITE["eval"]["bspline_scalar_sorted_$n"] = @benchmarkable sum_scalar_queries($spline, $queries)
    SUITE["eval"]["bspline_batch_sorted_$n"] = @benchmarkable $spline($out, $queries)
end

# =============================================================================
# Derivatives
# =============================================================================

SUITE["derivative"] = BenchmarkGroup()
SUITE["derivative"]["cubic"] = @benchmarkable DataInterpolations.derivative(
    $cubic, 5.4321
)
SUITE["derivative"]["linear"] = @benchmarkable DataInterpolations.derivative(
    $lin, 5.4321
)

# =============================================================================
# Vector-valued output
# =============================================================================

SUITE["vector_valued"] = BenchmarkGroup()

umat = rand(rng, 5, 200)
lin_vv = LinearInterpolation(umat, t)
SUITE["vector_valued"]["construct"] = @benchmarkable LinearInterpolation($umat, $t)
SUITE["vector_valued"]["eval"] = @benchmarkable $lin_vv(5.4321)
