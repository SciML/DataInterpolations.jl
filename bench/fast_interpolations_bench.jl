#=
FastInterpolations.jl benchmark, ported from their `benchmark/interpolation_benchmark.jl`
(`ProjectTorreyPines/FastInterpolations.jl`, commit 616b106b at the time of import).

This is the comparison they advertise on their README:
  - Compares Interpolations.jl, DataInterpolations.jl, FastInterpolations.jl (+ their
    Series interpolant), and Dierckx.jl.
  - Workload: `mpert × mpert` independent 1D interpolants over the same uniform
    `range(0.0, 1.0; length = npsi)` grid, evaluated at `n_eval` query points clustered
    near psi = 0 (cubic spacing). Mimics fusion-physics matrix-of-interpolants workloads.
  - Default config (matching their `--default`): `npsi = 64`, `mpert = 100` → 10_000
    interpolants per package, `n_eval = 1000` query points → 10⁷ total scalar evaluations.

Usage:
    julia +1.11 --project=bench bench/fast_interpolations_bench.jl
    julia +1.11 --project=bench bench/fast_interpolations_bench.jl --linear --small
    julia +1.11 --project=bench bench/fast_interpolations_bench.jl --cubic --default

This emits stdout-only output (no markdown report). The numbers feed the comparison
table in `bench/cross_library_comparison.md`.
=#

import Pkg
const BENCH_DIR = @__DIR__
Pkg.activate(BENCH_DIR)

using BenchmarkTools
using Interpolations
using DataInterpolations
using FastInterpolations
using Dierckx
using Random
using Printf
using Statistics

const SIZE_PRESETS = Dict(
    :tiny => (16, 2, 5),
    :small => (64, 5, 100),
    :default => (64, 100, 1000),
    :large => (64, 200, 4000),
)

const METHOD_OPTIONS = [:constant, :linear, :quadratic, :cubic]

function parse_args(args)
    size_key = :default
    method_key = :cubic
    for arg in args
        if startswith(arg, "--")
            key = Symbol(arg[3:end])
            if haskey(SIZE_PRESETS, key)
                size_key = key
            elseif key in METHOD_OPTIONS
                method_key = key
            end
        end
    end
    return size_key, method_key
end

const (SIZE_KEY, METHOD_KEY) = parse_args(ARGS)
const (NPSI, MPERT, N_EVAL_POINTS) = SIZE_PRESETS[SIZE_KEY]

function generate_test_data(npsi::Int, mpert::Int; seed::Int = 42)
    Random.seed!(seed)
    psi_grid = range(0.0, 1.0, length = npsi)
    data = rand(npsi, mpert, mpert)
    return psi_grid, data
end

function generate_evaluation_points(n_points::Int)
    return collect(range(0.0, 1.0, length = n_points)) .^ 3
end

# ---- Interpolations.jl -----------------------------------------------------
function init_interpolations(::Val{:linear}, psi_grid::AbstractRange, data::Array{Float64, 3})
    _, mpert, _ = size(data)
    first_itp = Interpolations.linear_interpolation(psi_grid, data[:, 1, 1])
    interps = Matrix{typeof(first_itp)}(undef, mpert, mpert)
    interps[1, 1] = first_itp
    for m1 in 1:mpert, m2 in 1:mpert
        (m1 == 1 && m2 == 1) && continue
        interps[m1, m2] = Interpolations.linear_interpolation(psi_grid, data[:, m1, m2])
    end
    return interps
end

function init_interpolations(::Val{:cubic}, psi_grid::AbstractRange, data::Array{Float64, 3})
    _, mpert, _ = size(data)
    first_itp = Interpolations.cubic_spline_interpolation(psi_grid, data[:, 1, 1])
    interps = Matrix{typeof(first_itp)}(undef, mpert, mpert)
    interps[1, 1] = first_itp
    for m1 in 1:mpert, m2 in 1:mpert
        (m1 == 1 && m2 == 1) && continue
        interps[m1, m2] = Interpolations.cubic_spline_interpolation(psi_grid, data[:, m1, m2])
    end
    return interps
end

function init_interpolations(::Val{:constant}, psi_grid::AbstractRange, data::Array{Float64, 3})
    _, mpert, _ = size(data)
    first_itp = Interpolations.constant_interpolation(psi_grid, data[:, 1, 1])
    interps = Matrix{typeof(first_itp)}(undef, mpert, mpert)
    interps[1, 1] = first_itp
    for m1 in 1:mpert, m2 in 1:mpert
        (m1 == 1 && m2 == 1) && continue
        interps[m1, m2] = Interpolations.constant_interpolation(psi_grid, data[:, m1, m2])
    end
    return interps
end

function init_interpolations(::Val{:quadratic}, psi_grid::AbstractRange, data::Array{Float64, 3})
    _, mpert, _ = size(data)
    knots = (psi_grid,)
    first_itp = Interpolations.extrapolate(
        Interpolations.scale(
            Interpolations.interpolate(
                data[:, 1, 1],
                Interpolations.BSpline(Interpolations.Quadratic(Interpolations.Reflect(Interpolations.OnCell()))),
            ),
            knots,
        ),
        Interpolations.Throw(),
    )
    interps = Matrix{typeof(first_itp)}(undef, mpert, mpert)
    interps[1, 1] = first_itp
    for m1 in 1:mpert, m2 in 1:mpert
        (m1 == 1 && m2 == 1) && continue
        interps[m1, m2] = Interpolations.extrapolate(
            Interpolations.scale(
                Interpolations.interpolate(
                    data[:, m1, m2],
                    Interpolations.BSpline(Interpolations.Quadratic(Interpolations.Reflect(Interpolations.OnCell()))),
                ),
                knots,
            ),
            Interpolations.Throw(),
        )
    end
    return interps
end

# Scalar-loop and broadcast evaluators
function eval_interpolations!(A, interps, psi)
    @inbounds for m2 in axes(interps, 2), m1 in axes(interps, 1)
        A[m1, m2] = interps[m1, m2](psi)
    end
    return A
end

run_interpolations_loop!(A, interps, psis) = (
    for psi in psis
        eval_interpolations!(A, interps, psi)
    end; A
)
function run_interpolations_broadcast!(A_all, interps, psis)
    @inbounds for m2 in axes(interps, 2), m1 in axes(interps, 1)
        @. A_all[:, m1, m2] = interps[m1, m2](psis)
    end
    return A_all
end

# ---- FastInterpolations.jl -----------------------------------------------
const FI_INIT_TABLE = Dict(
    :linear => FastInterpolations.linear_interp,
    :cubic => FastInterpolations.cubic_interp,
    :quadratic => FastInterpolations.quadratic_interp,
    :constant => FastInterpolations.constant_interp,
)

function init_fi(method::Symbol, psi_grid, data::Array{Float64, 3})
    f = FI_INIT_TABLE[method]
    _, mpert, _ = size(data)
    first_itp = f(psi_grid, data[:, 1, 1])
    interps = Matrix{typeof(first_itp)}(undef, mpert, mpert)
    interps[1, 1] = first_itp
    for m1 in 1:mpert, m2 in 1:mpert
        (m1 == 1 && m2 == 1) && continue
        interps[m1, m2] = f(psi_grid, data[:, m1, m2])
    end
    return interps
end

function eval_fi!(A, interps, psi)
    @inbounds for m2 in axes(interps, 2), m1 in axes(interps, 1)
        A[m1, m2] = interps[m1, m2](psi)
    end
    return A
end

run_fi_loop!(A, interps, psis) = (
    for psi in psis
        eval_fi!(A, interps, psi)
    end; A
)
function run_fi_vector!(A_all, interps, psis)
    @inbounds for m2 in axes(interps, 2), m1 in axes(interps, 1)
        @views interps[m1, m2](A_all[:, m1, m2], psis)
    end
    return A_all
end

# Series API (one anchor per query, shared across series)
function init_fi_series(method::Symbol, psi_grid, data::Array{Float64, 3})
    f = FI_INIT_TABLE[method]
    _, mpert, _ = size(data)
    ys = Series([data[:, m1, m2] for m2 in 1:mpert for m1 in 1:mpert])
    return f(psi_grid, ys)
end

run_fi_series_loop!(A, sitp, psis) = (
    for psi in psis
        sitp(A, psi)
    end; A
)
run_fi_series_vector!(A_all, sitp, psis) = (sitp(A_all, psis); A_all)

# ---- DataInterpolations.jl -----------------------------------------------
const DI_INIT_TABLE = Dict(
    :linear => DataInterpolations.LinearInterpolation,
    :cubic => DataInterpolations.CubicSpline,
    :quadratic => DataInterpolations.QuadraticInterpolation,
    :constant => DataInterpolations.ConstantInterpolation,
)

function init_di(method::Symbol, psi_grid, data::Array{Float64, 3})
    f = DI_INIT_TABLE[method]
    _, mpert, _ = size(data)
    t = collect(psi_grid)
    first_itp = f(data[:, 1, 1], t)
    interps = Matrix{typeof(first_itp)}(undef, mpert, mpert)
    interps[1, 1] = first_itp
    for m1 in 1:mpert, m2 in 1:mpert
        (m1 == 1 && m2 == 1) && continue
        interps[m1, m2] = f(data[:, m1, m2], t)
    end
    return interps
end

eval_di!(A, interps, psi) = (
    @inbounds for m2 in axes(interps, 2), m1 in axes(interps, 1)
        A[m1, m2] = interps[m1, m2](psi)
    end; A
)
run_di_loop!(A, interps, psis) = (
    for psi in psis
        eval_di!(A, interps, psi)
    end; A
)
function run_di_vector!(A_all, interps, psis)
    @inbounds for m2 in axes(interps, 2), m1 in axes(interps, 1)
        @views interps[m1, m2](A_all[:, m1, m2], psis)
    end
    return A_all
end

# ---- Dierckx.jl ---------------------------------------------------------
function init_dierckx(method::Symbol, psi_grid, data::Array{Float64, 3})
    method == :constant && return nothing  # Dierckx has no k=0
    k = Dict(:linear => 1, :quadratic => 2, :cubic => 3)[method]
    _, mpert, _ = size(data)
    t = collect(psi_grid)
    first_itp = Dierckx.Spline1D(t, data[:, 1, 1]; k = k, s = 0.0)
    interps = Matrix{typeof(first_itp)}(undef, mpert, mpert)
    interps[1, 1] = first_itp
    for m1 in 1:mpert, m2 in 1:mpert
        (m1 == 1 && m2 == 1) && continue
        interps[m1, m2] = Dierckx.Spline1D(t, data[:, m1, m2]; k = k, s = 0.0)
    end
    return interps
end

eval_dierckx!(A, interps, psi) = (
    @inbounds for m2 in axes(interps, 2), m1 in axes(interps, 1)
        A[m1, m2] = interps[m1, m2](psi)
    end; A
)
run_dierckx_loop!(A, interps, psis) = (
    for psi in psis
        eval_dierckx!(A, interps, psi)
    end; A
)
function run_dierckx_vector!(A_all, interps, psis)
    @inbounds for m2 in axes(interps, 2), m1 in axes(interps, 1)
        @views A_all[:, m1, m2] .= interps[m1, m2](psis)
    end
    return A_all
end

# ---- Driver ----------------------------------------------------------------
function bench_one(label, f; samples = 5, evals = 2, seconds = 120)
    f()  # warm-up
    GC.gc()
    b = @benchmark $f() samples = samples evals = evals seconds = seconds
    return median(b).time / 1.0e6   # ms
end

function run_fi_bench()
    method = METHOD_KEY
    psi_grid, data = generate_test_data(NPSI, MPERT)
    psis = generate_evaluation_points(N_EVAL_POINTS)
    n_interps = MPERT * MPERT
    total_evals = N_EVAL_POINTS * n_interps

    println("="^80)
    println("FastInterpolations.jl-style benchmark · method=$(method) · npsi=$(NPSI), mpert=$(MPERT), n_eval=$(N_EVAL_POINTS)")
    println("Reproducing `ProjectTorreyPines/FastInterpolations.jl/benchmark/interpolation_benchmark.jl`")
    println("="^80)
    println()

    A = Matrix{Float64}(undef, MPERT, MPERT)
    A_all = Array{Float64, 3}(undef, N_EVAL_POINTS, MPERT, MPERT)
    A_series = Vector{Float64}(undef, n_interps)
    A_series_all = [Vector{Float64}(undef, N_EVAL_POINTS) for _ in 1:n_interps]

    results = Dict{String, Tuple{Float64, Float64}}()  # name => (init_ms, eval_ms)

    function report(name, init_ms, eval_ms)
        results[name] = (init_ms, eval_ms)
        total = init_ms + eval_ms
        evs = total_evals / (eval_ms / 1.0e3)
        return @printf("  %-44s  init %8.3f ms  eval %8.3f ms  total %8.3f ms  evals/s %.2e\n", name, init_ms, eval_ms, total, evs)
    end

    # Interpolations.jl
    init_ms = bench_one("ITP init", () -> init_interpolations(Val(method), psi_grid, data))
    itp_interps = init_interpolations(Val(method), psi_grid, data)
    eval_ms = bench_one("ITP scalar", () -> run_interpolations_loop!(A, itp_interps, psis))
    report("Interpolations.jl (scalar)", init_ms, eval_ms)
    eval_ms = bench_one("ITP broadcast", () -> run_interpolations_broadcast!(A_all, itp_interps, psis))
    report("Interpolations.jl (broadcast)", init_ms, eval_ms)

    # FastInterpolations.jl
    init_ms = bench_one("FI init", () -> init_fi(method, psi_grid, data))
    fi_interps = init_fi(method, psi_grid, data)
    eval_ms = bench_one("FI scalar", () -> run_fi_loop!(A, fi_interps, psis))
    report("FastInterpolations.jl (scalar)", init_ms, eval_ms)
    eval_ms = bench_one("FI vector", () -> run_fi_vector!(A_all, fi_interps, psis))
    report("FastInterpolations.jl (vector)", init_ms, eval_ms)

    # FastInterpolations Series
    init_ms = bench_one("FI Series init", () -> init_fi_series(method, psi_grid, data))
    sitp = init_fi_series(method, psi_grid, data)
    eval_ms = bench_one("FI Series scalar", () -> run_fi_series_loop!(A_series, sitp, psis))
    report("FastInterpolations.jl (Series+scalar)", init_ms, eval_ms)
    eval_ms = bench_one("FI Series vector", () -> run_fi_series_vector!(A_series_all, sitp, psis))
    report("FastInterpolations.jl (Series+vector)", init_ms, eval_ms)

    # DataInterpolations.jl
    init_ms = bench_one("DI init", () -> init_di(method, psi_grid, data))
    di_interps = init_di(method, psi_grid, data)
    eval_ms = bench_one("DI scalar", () -> run_di_loop!(A, di_interps, psis))
    report("DataInterpolations.jl (scalar)", init_ms, eval_ms)
    eval_ms = bench_one("DI vector", () -> run_di_vector!(A_all, di_interps, psis))
    report("DataInterpolations.jl (vector)", init_ms, eval_ms)

    # Dierckx (if applicable)
    if method != :constant
        init_ms = bench_one("Dierckx init", () -> init_dierckx(method, psi_grid, data))
        drx_interps = init_dierckx(method, psi_grid, data)
        eval_ms = bench_one("Dierckx scalar", () -> run_dierckx_loop!(A, drx_interps, psis))
        report("Dierckx.jl (scalar)", init_ms, eval_ms)
        eval_ms = bench_one("Dierckx vector", () -> run_dierckx_vector!(A_all, drx_interps, psis))
        report("Dierckx.jl (vector)", init_ms, eval_ms)
    end

    # Summary
    println()
    println("="^80)
    println("Summary: speedup vs DataInterpolations.jl (scalar)")
    println("="^80)
    baseline_total = sum(results["DataInterpolations.jl (scalar)"])
    @printf("  %-44s  %10s  %10s  %10s  %s\n", "Package", "Init (ms)", "Eval (ms)", "Total (ms)", "Speedup")
    for (name, (init_ms, eval_ms)) in sort(collect(results); by = x -> sum(x[2]))
        total = init_ms + eval_ms
        @printf("  %-44s  %10.3f  %10.3f  %10.3f  %.2fx\n", name, init_ms, eval_ms, total, baseline_total / total)
    end
    println()

    return results
end

run_fi_bench()
