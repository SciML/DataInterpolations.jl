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

# =============================================================================
# Vector-valued output
# =============================================================================

SUITE["vector_valued"] = BenchmarkGroup()

umat = rand(rng, 5, 200)
lin_vv = LinearInterpolation(umat, t)
SUITE["vector_valued"]["construct"] = @benchmarkable LinearInterpolation($umat, $t)
SUITE["vector_valued"]["eval"] = @benchmarkable $lin_vv(5.4321)

# =============================================================================
# Knots that do not look linear: by default each lookup starts from the previous
# interval; `cache_search = false` uses an unhinted search (#618). The default cases
# compare against earlier revisions; the `cache_search = false` ones need the keyword.
# =============================================================================

SUITE["irregular"] = BenchmarkGroup()

for n in (100, 10_000)
    t_irr = collect(100 .* range(0.0, 1.0; length = n) .^ 3)
    u_irr = sin.(t_irr ./ 7)
    x_random = 100 .* rand(rng, 1000)
    x_sorted = sort(x_random)
    for (name, method) in (("linear", LinearInterpolation), ("cubic", CubicSpline))
        modes = Any[("default", (;))]
        hasfield(method, :cache_search) &&
            push!(modes, ("cache_search_false", (; cache_search = false)))
        for (mode, kw) in modes
            A = method(u_irr, t_irr; kw...)
            SUITE["irregular"]["$(name)_$(n)_$(mode)_sorted"] = @benchmarkable sum($A, $x_sorted)
            SUITE["irregular"]["$(name)_$(n)_$(mode)_random"] = @benchmarkable sum($A, $x_random)
        end
    end
end
