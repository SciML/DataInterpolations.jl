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

# Large-n BSpline derivatives: only the O(degree) nonzero window should contribute.
const n_bs = 10_000
const t_bs = collect(range(0.0, 100.0; length = n_bs))
const u_bs = sin.(t_bs)
const bs_interp = BSplineInterpolation(u_bs, t_bs, 3, :Average)
const bs_approx = BSplineApprox(u_bs, t_bs, 3, n_bs ÷ 3, :Average)
const t_bs_mid = t_bs[n_bs ÷ 2] + 0.01
SUITE["derivative"]["bspline_interp_d1"] = @benchmarkable DataInterpolations.derivative(
    $bs_interp, $t_bs_mid, 1
)
SUITE["derivative"]["bspline_interp_d2"] = @benchmarkable DataInterpolations.derivative(
    $bs_interp, $t_bs_mid, 2
)
SUITE["derivative"]["bspline_approx_d1"] = @benchmarkable DataInterpolations.derivative(
    $bs_approx, $t_bs_mid, 1
)

# =============================================================================
# Vector-valued output
# =============================================================================

SUITE["vector_valued"] = BenchmarkGroup()

umat = rand(rng, 5, 200)
lin_vv = LinearInterpolation(umat, t)
SUITE["vector_valued"]["construct"] = @benchmarkable LinearInterpolation($umat, $t)
SUITE["vector_valued"]["eval"] = @benchmarkable $lin_vv(5.4321)
