#=
Cross-library 1D interpolation benchmark for DataInterpolations.jl.

Compares DataInterpolations.jl (PR #529 branch, with cached Auto(t_props)) against
Interpolations.jl, Dierckx.jl, BasicInterpolators.jl, and PCHIPInterpolation.jl.

Usage:
    julia +1.11 --project=bench bench/cross_library_comparison.jl

The script writes a fresh `bench/cross_library_comparison.md` with the report.
=#

import Pkg
const BENCH_DIR = @__DIR__
const REPO_ROOT = dirname(BENCH_DIR)
Pkg.activate(BENCH_DIR)

using Printf
using Random
using Statistics
using BenchmarkTools
using LinearAlgebra
using InteractiveUtils: versioninfo

using DataInterpolations
using Interpolations
using Dierckx
using BasicInterpolators
using PCHIPInterpolation

const DI = DataInterpolations
const ITP = Interpolations
const DRX = Dierckx
const BI = BasicInterpolators
const PCHIP = PCHIPInterpolation

const RNG = MersenneTwister(0x00C0FFEE)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

"""
Build the test data `(u, t)` of length `n`. `pattern` is `:uniform` or `:nonuniform`.
The samples are deterministic for a given (n, pattern).
"""
function make_knots(n::Int, pattern::Symbol)
    rng = MersenneTwister(0x0BADBEEF + n + (pattern === :uniform ? 0 : 1))
    if pattern === :uniform
        # Keep as a `range` so Interpolations.cubic_spline_interpolation accepts it.
        # Concrete Vector form for libraries that need it is taken via `collect(t)` later.
        t = range(0.0, 1.0; length = n)
    elseif pattern === :nonuniform
        t = sort(rand(rng, n))
        t .= (t .- first(t)) ./ (last(t) - first(t))
    else
        error("unknown pattern $(pattern)")
    end
    u = @. sin(2π * t) + 0.3 * cos(7 * t)
    return u, t
end

"""
Build a query batch of length `m` in the same domain. `pattern` is `:sorted`,
`:random`, or `:chained` (monotone, ODE-style).
"""
function make_queries(m::Int, pattern::Symbol)
    rng = MersenneTwister(0x00C0FFEE + m + (pattern === :sorted ? 0 : pattern === :random ? 1 : 2))
    if pattern === :sorted
        tt = sort(rand(rng, m))
    elseif pattern === :random
        tt = rand(rng, m)
    elseif pattern === :chained
        steps = rand(rng, m)
        tt = cumsum(steps)
        tt .= (tt .- first(tt)) ./ (last(tt) - first(tt)) .* 0.999 .+ 0.0005
    else
        error("unknown query pattern $(pattern)")
    end
    return tt
end

# Pretty-print BenchmarkTools.Trial median and IQR (q3-q1) in human units
function fmt_trial(t::BenchmarkTools.Trial)
    med = median(t).time # ns
    q1 = quantile(t.times, 0.25)
    q3 = quantile(t.times, 0.75)
    iqr = q3 - q1
    return string(BenchmarkTools.prettytime(med), " (IQR ", BenchmarkTools.prettytime(iqr), ")")
end

fmt_trial(::Nothing) = "—"

# Result store: Dict{(algorithm, case, library) => Dict{params => Trial}}
# We'll just use nested dictionaries to keep it simple
const RESULTS = Dict{String, Any}()

function record!(category::String, algo::String, lib::String, params::String, trial)
    key = string(category, " | ", algo)
    d = get!(RESULTS, key, Dict{Tuple{String, String}, Any}())
    d[(lib, params)] = trial
    return nothing
end

# Limit each individual benchmark; full sweep is large.
const BMK_SECONDS = 0.5
const BMK_SAMPLES = 100

function bench(expr_setup::Function, args...; seconds = BMK_SECONDS, samples = BMK_SAMPLES)
    b = @benchmarkable $(expr_setup)($(args)...) evals = 1 samples = samples seconds = seconds
    return run(b)
end

# -----------------------------------------------------------------------------
# Library adapter functions
#
# Each adapter returns the interpolator for the given (u, t).
# We separate construction from evaluation so we can benchmark each.
# -----------------------------------------------------------------------------

# --- Linear --------------------------------------------------------------
# Helpers: many third-party libraries want `Vector{Float64}` knots, not a `range`.
_vec(t) = collect(Float64, t)

build_di_linear(u, t) = DI.LinearInterpolation(_vec(u), _vec(t))
build_itp_linear_uniform(u, t) = ITP.linear_interpolation(t, u) # t may be a range
build_itp_linear_nonuniform(u, t) = ITP.interpolate((_vec(t),), _vec(u), ITP.Gridded(ITP.Linear()))
build_drx_linear(u, t) = DRX.Spline1D(_vec(t), _vec(u); k = 1, bc = "extrapolate")
build_bi_linear(u, t) = BI.LinearInterpolator(_vec(t), _vec(u), BI.WeakBoundaries())

# --- Cubic spline (natural BC for DI, Interpolations) --------------------
build_di_cubic(u, t) = DI.CubicSpline(_vec(u), _vec(t))
build_itp_cubic_uniform(u, t) = ITP.cubic_spline_interpolation(t, u)
build_drx_cubic(u, t) = DRX.Spline1D(_vec(t), _vec(u); k = 3, bc = "extrapolate")
build_bi_cubic(u, t) = BI.CubicSplineInterpolator(_vec(t), _vec(u), BI.WeakBoundaries())

# --- Quadratic spline ----------------------------------------------------
build_di_quadratic(u, t) = DI.QuadraticSpline(_vec(u), _vec(t))
build_drx_quadratic(u, t) = DRX.Spline1D(_vec(t), _vec(u); k = 2, bc = "extrapolate")

# --- Akima ---------------------------------------------------------------
build_di_akima(u, t) = DI.AkimaInterpolation(_vec(u), _vec(t))

# --- PCHIP / monotone cubic Hermite -------------------------------------
function build_di_cubic_hermite(u, t)
    uv = _vec(u)
    tv = _vec(t)
    n = length(tv)
    du = similar(uv)
    @inbounds for i in 2:(n - 1)
        du[i] = (uv[i + 1] - uv[i - 1]) / (tv[i + 1] - tv[i - 1])
    end
    du[1] = (uv[2] - uv[1]) / (tv[2] - tv[1])
    du[n] = (uv[n] - uv[n - 1]) / (tv[n] - tv[n - 1])
    return DI.CubicHermiteSpline(du, uv, tv)
end
build_pchip(u, t) = PCHIP.Interpolator(_vec(t), _vec(u))

# --- Single-eval dispatch -----------------------------------------------
# DI, Dierckx, BasicInterpolators all use `A(x)`
# Interpolations also uses `A(x)`
# PCHIP uses `A(x)`
single_eval(A, x) = A(x)

# --- Batched eval --------------------------------------------------------
# DI: `A(out, tt)` in-place; uses sorted-batch fast path when sorted.
batched_eval_di!(out, A, tt) = A(out, tt)
# Interpolations: broadcast (no native batched call)
batched_eval_itp!(out, A, tt) = (out .= A.(tt); out)
# Dierckx: has Spline1D batched evaluation (returns a new array). We'll call the
# in-place form if we can find it; otherwise wrap broadcast.
batched_eval_drx!(out, A::DRX.Spline1D, tt) = (out .= DRX.evaluate.(Ref(A), tt); out)
# BasicInterpolators: broadcast
batched_eval_bi!(out, A, tt) = (out .= A.(tt); out)
# PCHIP: broadcast
batched_eval_pchip!(out, A, tt) = (out .= A.(tt); out)

# -----------------------------------------------------------------------------
# Verification: every library should agree on the same query to within tol
# -----------------------------------------------------------------------------

function verify_agreement(algo, builders, u, t, tt; tol = 1.0e-6, ref_name = first(keys(builders)))
    ref_build = builders[ref_name]
    ref = ref_build(u, t)
    ref_vals = [single_eval(ref, x) for x in tt]
    for (name, build) in builders
        name == ref_name && continue
        A = build(u, t)
        vals = [single_eval(A, x) for x in tt]
        diff = maximum(abs.(ref_vals .- vals))
        if diff > tol
            @warn "Library $name disagrees with $ref_name on $algo by max diff $diff (tol=$tol)"
        end
    end
    return nothing
end

# -----------------------------------------------------------------------------
# Run sweep
# -----------------------------------------------------------------------------

const SMOKE = get(ENV, "BENCH_SMOKE", "0") != "0"

const N_VALUES = SMOKE ? [100, 1_000] : [100, 1_000, 10_000, 100_000]
const M_VALUES = SMOKE ? [10, 1_000] : [1, 10, 1_000, 100_000]
const KNOT_PATTERNS = [:uniform, :nonuniform]

# Algorithm spec:
#   builders :: Dict{library => build_fun(u, t)}
#   batched! :: Dict{library => (out, A, tt) -> ...}
#   supports :: Dict{library => set of knot_patterns it supports}
const ALGORITHMS = let
    d = Dict{String, NamedTuple}()
    d["Linear"] = (
        builders = Dict(
            "DataInterpolations" => build_di_linear,
            "Interpolations (uniform)" => build_itp_linear_uniform,
            "Interpolations (gridded)" => build_itp_linear_nonuniform,
            "Dierckx (k=1)" => build_drx_linear,
            "BasicInterpolators" => build_bi_linear,
        ),
        batched! = Dict(
            "DataInterpolations" => batched_eval_di!,
            "Interpolations (uniform)" => batched_eval_itp!,
            "Interpolations (gridded)" => batched_eval_itp!,
            "Dierckx (k=1)" => batched_eval_drx!,
            "BasicInterpolators" => batched_eval_bi!,
        ),
        supports = Dict(
            "DataInterpolations" => [:uniform, :nonuniform],
            "Interpolations (uniform)" => [:uniform],
            "Interpolations (gridded)" => [:uniform, :nonuniform],
            "Dierckx (k=1)" => [:uniform, :nonuniform],
            "BasicInterpolators" => [:uniform, :nonuniform],
        ),
    )
    d["CubicSpline"] = (
        builders = Dict(
            "DataInterpolations" => build_di_cubic,
            "Interpolations (uniform)" => build_itp_cubic_uniform,
            "Dierckx (k=3)" => build_drx_cubic,
            "BasicInterpolators" => build_bi_cubic,
        ),
        batched! = Dict(
            "DataInterpolations" => batched_eval_di!,
            "Interpolations (uniform)" => batched_eval_itp!,
            "Dierckx (k=3)" => batched_eval_drx!,
            "BasicInterpolators" => batched_eval_bi!,
        ),
        supports = Dict(
            "DataInterpolations" => [:uniform, :nonuniform],
            "Interpolations (uniform)" => [:uniform],
            "Dierckx (k=3)" => [:uniform, :nonuniform],
            "BasicInterpolators" => [:uniform, :nonuniform],
        ),
    )
    d["QuadraticSpline"] = (
        builders = Dict(
            "DataInterpolations" => build_di_quadratic,
            "Dierckx (k=2)" => build_drx_quadratic,
        ),
        batched! = Dict(
            "DataInterpolations" => batched_eval_di!,
            "Dierckx (k=2)" => batched_eval_drx!,
        ),
        supports = Dict(
            "DataInterpolations" => [:uniform, :nonuniform],
            "Dierckx (k=2)" => [:uniform, :nonuniform],
        ),
    )
    d["Akima"] = (
        builders = Dict("DataInterpolations" => build_di_akima),
        batched! = Dict("DataInterpolations" => batched_eval_di!),
        supports = Dict("DataInterpolations" => [:uniform, :nonuniform]),
    )
    d["MonotoneCubic"] = (
        builders = Dict(
            "DataInterpolations (CubicHermite)" => build_di_cubic_hermite,
            "PCHIPInterpolation" => build_pchip,
        ),
        batched! = Dict(
            "DataInterpolations (CubicHermite)" => batched_eval_di!,
            "PCHIPInterpolation" => batched_eval_pchip!,
        ),
        supports = Dict(
            "DataInterpolations (CubicHermite)" => [:uniform, :nonuniform],
            "PCHIPInterpolation" => [:uniform, :nonuniform],
        ),
    )
    d
end

# Quick agreement check (use small n)
function run_agreement_checks()
    println("\n== Cross-library agreement check ==")
    tt_check = collect(range(0.05, 0.95; length = 11))
    for knot in (:uniform, :nonuniform)
        u, t = make_knots(50, knot)
        for (algo, spec) in ALGORITHMS
            # Filter builders that support this knot pattern
            supported = Dict(k => v for (k, v) in spec.builders if knot in spec.supports[k])
            length(supported) < 2 && continue
            # Different algorithms can disagree slightly (different BCs). We use
            # a moderately loose tolerance only meant to catch outright bugs.
            tol = algo == "Linear" ? 1.0e-10 : (algo in ("CubicSpline", "QuadraticSpline") ? 5.0e-2 : 1.0e-1)
            try
                verify_agreement(algo, supported, u, t, tt_check; tol = tol)
            catch e
                @warn "Agreement check failed: $algo $knot" exception = (e, catch_backtrace())
            end
        end
    end
    return nothing
end

# -----------------------------------------------------------------------------
# Bench cases
# -----------------------------------------------------------------------------

function warmup()
    # Force compilation of every builder + dispatch path so the actual benchmarks
    # measure work, not first-time compile.
    println("Warming up…")
    for knot in KNOT_PATTERNS
        u, t = make_knots(100, knot)
        tt = make_queries(8, :sorted)
        out = similar(tt)
        for (_algo, spec) in ALGORITHMS
            for (lib, build) in spec.builders
                knot in spec.supports[lib] || continue
                A = try
                    build(u, t)
                catch
                    continue
                end
                single_eval(A, 0.5)
                try
                    spec.batched![lib](out, A, tt)
                catch
                end
            end
        end
    end
    return nothing
end

function run_construction()
    println("\n== Construction time ==")
    for (algo, spec) in ALGORITHMS
        for n in N_VALUES
            for knot in KNOT_PATTERNS
                u, t = make_knots(n, knot)
                for (lib, build) in spec.builders
                    knot in spec.supports[lib] || continue
                    trial = try
                        bench(build, u, t)
                    catch e
                        @warn "construction failed: $algo $lib n=$n knot=$knot" exception = e
                        nothing
                    end
                    record!(
                        "construction", algo, lib,
                        string("n=", n, ",", knot),
                        trial,
                    )
                    println(
                        @sprintf(
                            "  %-25s | %-30s | n=%-7d | %s | %s",
                            algo, lib, n, String(knot), fmt_trial(trial)
                        )
                    )
                end
            end
        end
    end
    return nothing
end

function run_single_query()
    println("\n== Single-query latency ==")
    x_query = 0.42718
    for (algo, spec) in ALGORITHMS
        for n in N_VALUES
            for knot in KNOT_PATTERNS
                u, t = make_knots(n, knot)
                for (lib, build) in spec.builders
                    knot in spec.supports[lib] || continue
                    A = try
                        build(u, t)
                    catch e
                        @warn "skipping single-query (build failed): $algo $lib n=$n" exception = e
                        continue
                    end
                    trial = try
                        bench(single_eval, A, x_query)
                    catch e
                        @warn "single-query failed: $algo $lib" exception = e
                        nothing
                    end
                    record!(
                        "single-query", algo, lib,
                        string("n=", n, ",", knot),
                        trial,
                    )
                    println(
                        @sprintf(
                            "  %-25s | %-30s | n=%-7d | %s | %s",
                            algo, lib, n, String(knot), fmt_trial(trial)
                        )
                    )
                end
            end
        end
    end
    return nothing
end

# For batched/random/chained, we hold n at a single representative knot pattern
# (uniform) — exploring the cross-product of n × knot-pattern × m × query-pattern
# would explode the runtime. The batched-call story is the same on non-uniform
# (DI's sorted fast-path still kicks in).

function run_batched(query_pattern::Symbol, m_values = M_VALUES)
    label = query_pattern === :sorted ? "Sorted batch" :
        query_pattern === :random ? "Random batch" :
        "Chained ODE-style"
    println("\n== $label ==")
    knot = :uniform
    for (algo, spec) in ALGORITHMS
        for n in N_VALUES
            u, t = make_knots(n, knot)
            for m in m_values
                tt = make_queries(m, query_pattern)
                for (lib, build) in spec.builders
                    knot in spec.supports[lib] || continue
                    A = try
                        build(u, t)
                    catch e
                        continue
                    end
                    out = similar(tt)
                    batched! = spec.batched![lib]
                    if query_pattern === :chained
                        # Sequential single-eval loop is what ODE solvers do.
                        # Don't use the batched interface for this case.
                        f = (A, tt) -> begin
                            s = 0.0
                            for x in tt
                                s += single_eval(A, x)
                            end
                            s
                        end
                        trial = try
                            bench(f, A, tt)
                        catch e
                            @warn "chained eval failed: $algo $lib" exception = e
                            nothing
                        end
                    else
                        trial = try
                            bench(batched!, out, A, tt)
                        catch e
                            @warn "batched eval failed: $algo $lib" exception = e
                            nothing
                        end
                    end
                    record!(
                        string(label), algo, lib,
                        string("n=", n, ",m=", m),
                        trial,
                    )
                    println(
                        @sprintf(
                            "  %-25s | %-30s | n=%-7d | m=%-7d | %s",
                            algo, lib, n, m, fmt_trial(trial)
                        )
                    )
                end
            end
        end
    end
    return nothing
end

# -----------------------------------------------------------------------------
# Markdown output
# -----------------------------------------------------------------------------

function library_order(libs)
    # Stable order: DI first, then others alphabetically
    di = filter(l -> startswith(l, "DataInterpolations"), libs)
    others = sort(setdiff(libs, di))
    return vcat(sort(di), others)
end

function unique_params(d)
    return sort!(unique([p for ((_, p), _) in d]))
end

function unique_libs(d)
    return sort!(unique([l for ((l, _), _) in d]))
end

function write_table(io::IO, key::AbstractString, d::Dict)
    libs = library_order(unique_libs(d))
    params = unique_params(d)
    print(io, "\n### ", key, "\n\n")
    # Header
    print(io, "| Library | ")
    for p in params
        print(io, p, " | ")
    end
    print(io, "\n")
    print(io, "|---|")
    for _ in params
        print(io, "---|")
    end
    print(io, "\n")
    for lib in libs
        print(io, "| ", lib, " | ")
        for p in params
            trial = get(d, (lib, p), nothing)
            print(io, fmt_trial(trial), " | ")
        end
        print(io, "\n")
    end
    return nothing
end

function write_report(path::String; total_seconds = 0.0)
    open(path, "w") do io
        println(io, "# Cross-library 1D interpolation benchmark")
        println(io)
        println(io, "## Setup")
        println(io)
        # Host/Julia info
        println(io, "```")
        io_buf = IOBuffer()
        versioninfo(io_buf)
        print(io, String(take!(io_buf)))
        println(io, "```")
        println(io)
        println(io, "Bench harness: `BenchmarkTools.@benchmark` with `evals=1`, max samples=$(BMK_SAMPLES), max seconds=$(BMK_SECONDS).")
        println(io)
        commit = read(`git -C $REPO_ROOT rev-parse HEAD`, String) |> strip
        println(io, "Commit: `$commit`")
        println(io)
        println(io, "Library versions:")
        println(io, "```")
        deps = Pkg.project().dependencies
        all_info = Pkg.dependencies()
        for pkg in ("DataInterpolations", "Interpolations", "Dierckx", "BasicInterpolators", "PCHIPInterpolation", "BenchmarkTools")
            if haskey(deps, pkg)
                v = all_info[deps[pkg]].version
                println(io, "  ", pkg, " ", v)
            end
        end
        println(io, "```")
        println(io)
        println(io, "Total bench time: ", round(total_seconds; digits = 1), " s")
        println(io)

        # Section 1: construction
        println(io, "## Construction time")
        println(io)
        println(io, "Rows = library, columns = (n, knot pattern). Values = median wall time (IQR).")
        for key in sort!(collect(keys(RESULTS)))
            startswith(key, "construction") || continue
            algo = split(key, " | ")[2]
            write_table(io, algo, RESULTS[key])
        end

        # Section 2: single-query latency
        println(io, "\n## Single-query latency")
        println(io)
        println(io, "Cold single evaluation `A(x_query)`. Rows = library, columns = (n, knot pattern).")
        for key in sort!(collect(keys(RESULTS)))
            startswith(key, "single-query") || continue
            algo = split(key, " | ")[2]
            write_table(io, algo, RESULTS[key])
        end

        # Section 3: sorted batch
        println(io, "\n## Sorted batch")
        println(io)
        println(io, "`A(out, tt)` where `tt` is sorted random points in domain. (knot pattern = uniform)")
        for key in sort!(collect(keys(RESULTS)))
            startswith(key, "Sorted batch") || continue
            algo = split(key, " | ")[2]
            write_table(io, algo, RESULTS[key])
        end

        # Section 4: random batch
        println(io, "\n## Random batch")
        println(io)
        println(io, "`A(out, tt)` where `tt` is unsorted. (knot pattern = uniform)")
        for key in sort!(collect(keys(RESULTS)))
            startswith(key, "Random batch") || continue
            algo = split(key, " | ")[2]
            write_table(io, algo, RESULTS[key])
        end

        # Section 5: chained
        println(io, "\n## Chained ODE-style")
        println(io)
        println(io, "Sequential `for x in tt; A(x); end` over a monotone sequence. (knot pattern = uniform)")
        for key in sort!(collect(keys(RESULTS)))
            startswith(key, "Chained ODE-style") || continue
            algo = split(key, " | ")[2]
            write_table(io, algo, RESULTS[key])
        end

        println(io, "\n## Reproducer")
        println(io)
        println(io, "Bench script: `bench/cross_library_comparison.jl`")
        println(io)
        println(io, "Bench Project.toml: `bench/Project.toml` (devs DI from `..`).")
        println(io)
        println(
            io,
            """
            To rerun:
            ```bash
            cd /home/crackauc/sandbox/tmp_20260515_091703_4914/DataInterpolations.jl
            git checkout fff-strategy-batched-evals
            julia +1.11 --project=bench bench/cross_library_comparison.jl
            ```
            """,
        )
    end
    return nothing
end

# -----------------------------------------------------------------------------
# Top-level
# -----------------------------------------------------------------------------

function main()
    println("Starting cross-library benchmark sweep…")
    println("Bench dir: ", BENCH_DIR)

    t_start = time()

    run_agreement_checks()
    warmup()

    run_construction()
    run_single_query()
    run_batched(:sorted)
    run_batched(:random)
    # Chained is most interesting at moderate m; cap to avoid 100k single-evals on Dierckx etc.
    run_batched(:chained, [1_000])

    t_end = time()
    total = t_end - t_start
    println("\nTotal bench time: ", round(total; digits = 1), " s")

    report_path = joinpath(BENCH_DIR, "cross_library_comparison.md")
    write_report(report_path; total_seconds = total)
    println("Wrote report to ", report_path)
    return nothing
end

main()
