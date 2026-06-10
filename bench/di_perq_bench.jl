#=
DI per-query micro-bench over the knot/query regimes where the search
strategy dominates:

   Workload                                | DI before | DI after | FastInterp
   --------------------------------------- | --------- | -------- | ----------
   Range knots, sorted queries             |    77     |    ?     |   3.5
   Uniform Vector knots, sorted queries    |    38     |    ?     |   27
   Non-uniform Vector knots, sorted        |    78     |    ?     |   n/a
   Range knots, shuffled random queries    |    87     |    ?     |   3.5
   Range knots, monotone ODE-chain         |    ?      |    ?     |   ?

n = 10_000, m = 1000, Float64.

Usage:
    julia +1.11 --project=bench bench/di_perq_bench.jl
=#

import Pkg
const BENCH_DIR = @__DIR__
Pkg.activate(BENCH_DIR)

using BenchmarkTools
using Random
using Statistics

using DataInterpolations
using FastInterpolations

const RNG = MersenneTwister(0x00C0FFEE)

const N = 10_000
const M = 1_000

# Knots
const range_knots = range(0.0, 1.0; length = N)
const uniform_vec_knots = collect(range_knots)
const nonuniform_vec_knots = sort!(rand(MersenneTwister(0x0BADBEEF), N))
const u_range = sin.(2π .* range_knots) .+ 0.3 .* cos.(7 .* range_knots)
const u_uniform = collect(u_range)
const u_nonuniform = sin.(2π .* nonuniform_vec_knots) .+ 0.3 .* cos.(7 .* nonuniform_vec_knots)

# Queries — keep in [knot_min, knot_max] so all libraries' interpolations
# stay in-domain. Clamp into the non-uniform vector range below.
function clamp_to(x, t)
    a, b = first(t), last(t)
    return @. clamp(x, a + (b - a) * 1.0e-6, b - (b - a) * 1.0e-6)
end

const queries_sorted_raw = sort!(rand(MersenneTwister(0x00C0FFEE + M), M))
const queries_random_raw = rand(MersenneTwister(0x00C0FFEE + M + 1), M)
const queries_chained_raw = let steps = rand(MersenneTwister(0x00C0FFEE + M + 2), M)
    tt = cumsum(steps)
    tt = (tt .- first(tt)) ./ (last(tt) - first(tt)) .* 0.999 .+ 0.0005
    tt
end

# Library builders
build_di_range(u, t) = DataInterpolations.LinearInterpolation(u, t)
build_di_uniform_vec(u, t) = DataInterpolations.LinearInterpolation(u, t)
build_di_nonuniform(u, t) = DataInterpolations.LinearInterpolation(u, t)
build_fi(u, t) = FastInterpolations.linear_interp(t, u)

# Per-query bench: tight loop with no batched call. Measures latency of A(x).
function per_query_loop(A, queries)
    s = 0.0
    @inbounds for k in eachindex(queries)
        s += A(queries[k])
    end
    return s
end

function bench_pq(A, queries; samples = 500, seconds = 5.0)
    # Warm-up call before benchmark — avoids first-call effects from carrying
    # over Guesser state across benchmarks for different `A`s.
    per_query_loop(A, queries)
    b = @benchmarkable per_query_loop($A, $queries) evals = 1 samples = samples seconds = seconds
    return run(b)
end

fmt_pq(t) = string(round(median(t).time / M; digits = 2), " ns/q")

println("== n=$N, m=$M, single-query latency (ns per query) ==\n")

println("Range knots:")
A_di_r = build_di_range(u_range, range_knots)
A_fi_r = build_fi(u_range, range_knots)
queries_sorted_r = clamp_to(queries_sorted_raw, range_knots)
queries_random_r = clamp_to(queries_random_raw, range_knots)
queries_chained_r = clamp_to(queries_chained_raw, range_knots)
println("  DataInterp (Auto)         | sorted   : ", fmt_pq(bench_pq(A_di_r, queries_sorted_r)))
println("  DataInterp (Auto)         | random   : ", fmt_pq(bench_pq(A_di_r, queries_random_r)))
println("  DataInterp (Auto)         | chained  : ", fmt_pq(bench_pq(A_di_r, queries_chained_r)))
println("  FastInterp                | sorted   : ", fmt_pq(bench_pq(A_fi_r, queries_sorted_r)))
println("  FastInterp                | random   : ", fmt_pq(bench_pq(A_fi_r, queries_random_r)))
println("  FastInterp                | chained  : ", fmt_pq(bench_pq(A_fi_r, queries_chained_r)))

println("\nUniform Vector knots (collect(range)):")
A_di_uv = build_di_uniform_vec(u_uniform, uniform_vec_knots)
A_fi_uv = build_fi(u_uniform, uniform_vec_knots)
queries_sorted_u = clamp_to(queries_sorted_raw, uniform_vec_knots)
queries_random_u = clamp_to(queries_random_raw, uniform_vec_knots)
queries_chained_u = clamp_to(queries_chained_raw, uniform_vec_knots)
println("  DataInterp (Auto)         | sorted   : ", fmt_pq(bench_pq(A_di_uv, queries_sorted_u)))
println("  DataInterp (Auto)         | random   : ", fmt_pq(bench_pq(A_di_uv, queries_random_u)))
println("  DataInterp (Auto)         | chained  : ", fmt_pq(bench_pq(A_di_uv, queries_chained_u)))
println("  FastInterp                | sorted   : ", fmt_pq(bench_pq(A_fi_uv, queries_sorted_u)))
println("  FastInterp                | random   : ", fmt_pq(bench_pq(A_fi_uv, queries_random_u)))

println("\nNon-uniform Vector knots:")
A_di_nv = build_di_nonuniform(u_nonuniform, nonuniform_vec_knots)
A_fi_nv = build_fi(u_nonuniform, nonuniform_vec_knots)
queries_sorted_n = clamp_to(queries_sorted_raw, nonuniform_vec_knots)
queries_random_n = clamp_to(queries_random_raw, nonuniform_vec_knots)
queries_chained_n = clamp_to(queries_chained_raw, nonuniform_vec_knots)
println("  DataInterp (Auto)         | sorted   : ", fmt_pq(bench_pq(A_di_nv, queries_sorted_n)))
println("  DataInterp (Auto)         | random   : ", fmt_pq(bench_pq(A_di_nv, queries_random_n)))
println("  DataInterp (Auto)         | chained  : ", fmt_pq(bench_pq(A_di_nv, queries_chained_n)))
println("  FastInterp                | sorted   : ", fmt_pq(bench_pq(A_fi_nv, queries_sorted_n)))
println("  FastInterp                | random   : ", fmt_pq(bench_pq(A_fi_nv, queries_random_n)))

println("\nReporting Auto kind selection:")
println("  Range knots                : Auto kind = ", A_di_r.strategy.kind)
println("  Uniform Vector knots       : Auto kind = ", A_di_uv.strategy.kind)
println("  Non-uniform Vector knots   : Auto kind = ", A_di_nv.strategy.kind)
