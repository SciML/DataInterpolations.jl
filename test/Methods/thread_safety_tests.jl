using DataInterpolations
using DataInterpolations: derivative, integral, invert_integral
using Test
using Base.Threads

# Regression tests for issue #532: a `BSplineInterpolation`, `BSplineApprox` or
# (default, uncached) `QuadraticSpline` kept a scratch buffer in its `sc` field
# and overwrote it in place on every evaluation. Because a "read" (evaluation)
# mutated shared state, evaluating one shared interpolant from several threads
# silently returned wrong values. Evaluation (and differentiation) must be
# reentrant: it must not write to any state stored in the interpolant, except
# the search hint `iguesser.idx_prev` that the default `cache_search = true`
# keeps on knots that do not look linear. With `cache_search = false` it writes
# nothing at all (issue #618).

function thread_safety_interpolants()
    u = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
    t = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    umat = [
        0.0 1.0 0.0 1.0 0.0 1.0
        1.0 0.0 1.0 0.0 1.0 0.0
    ]            # 2 × length(t): last axis indexes t
    return [
        ("BSplineInterpolation", BSplineInterpolation(u, t, 3, :Average)),
        # degree 1 == "linear" B-spline: same `_interpolate` method, so equally affected
        ("BSplineInterpolation (degree 1)", BSplineInterpolation(u, t, 1, :Average)),
        ("BSplineApprox", BSplineApprox(u, t, 3, 4, :Average)),
        ("BSplineInterpolation (array)", BSplineInterpolation(umat, t, 3, :Average)),
        ("BSplineApprox (array)", BSplineApprox(umat, t, 3, 4, :Average)),
        ("QuadraticSpline", QuadraticSpline(u, t)),                       # cache_parameters = false
        ("QuadraticSpline (cache_parameters)", QuadraticSpline(u, t; cache_parameters = true)),
    ]
end

@testset "Thread safety / reentrancy (issue #532)" begin
    pts = collect(range(1.0, 6.0; length = 41))   # includes the knots themselves

    # Deterministic guard (catches the bug even on a single thread): repeated and
    # interleaved evaluation must leave the interpolant's scratch buffer — and its
    # results — unchanged.
    @testset "evaluation does not mutate the interpolant: $name" for (name, A) in
        thread_safety_interpolants()
        sc_before = copy(A.sc)
        ref = [A(x) for x in pts]
        dref = [derivative(A, x) for x in pts]
        for _ in 1:200, x in pts
            A(x)
            derivative(A, x)
        end
        @test A.sc == sc_before
        @test [A(x) for x in pts] == ref
        @test [derivative(A, x) for x in pts] == dref
    end

    # End-to-end check: many threads hammering one shared interpolant must all
    # agree with the single-threaded reference. Only meaningful with >1 thread.
    if Threads.nthreads() > 1
        @testset "concurrent evaluation is race-free: $name" for (name, A) in
            thread_safety_interpolants()
            ref = [A(x) for x in pts]
            wrong = Threads.Atomic{Int}(0)
            Threads.@threads for _ in 1:50_000
                for (j, x) in enumerate(pts)
                    A(x) == ref[j] || Threads.atomic_add!(wrong, 1)
                end
            end
            @test wrong[] == 0
        end
    else
        @info "issue #532: concurrent stress test skipped; start Julia with more " *
            "than one thread (e.g. `julia -t auto`) to exercise it."
    end
end

# Issue #618: the cached search hint `iguesser.idx_prev` is only read when the
# knots do not look linear, so evaluation on linear-looking knots must not write it.
@testset "evaluation on linear-looking knots leaves idx_prev untouched (issue #618)" begin
    u = sin.(1:12)
    for t in (collect(1.0:12.0), collect(1.0:12.0) .+ 0.001 .* sin.(1:12)),
            A in (
                LinearInterpolation(u, t), LinearInterpolation(u, t; cache_parameters = true),
                QuadraticInterpolation(u, t), CubicSpline(u, t),
                ConstantInterpolation(u, t), AkimaInterpolation(u, t),
            )
        @test A.iguesser.linear_lookup
        A.iguesser.idx_prev[] = 1
        for x in range(first(t), last(t); length = 61)
            A(x)
            derivative(A, x)
        end
        @test A.iguesser.idx_prev[] == 1
    end
end

# Issue #618, knots that do not look linear: the default starts each lookup from the
# interval found by the previous one, stored in `iguesser.idx_prev`. With
# `cache_search = false` the lookup is unhinted and evaluation writes nothing.

# Whether every value reachable from `a` and `b` (fields, array contents, `Ref`
# contents) is equal, so comparing with a `deepcopy` snapshot detects any write.
# Seed the hint with `0` first: lookups store indices `>= 1`, so a write cannot
# leave it unchanged by coincidence.
function same_state(a, b)
    typeof(a) === typeof(b) || return false
    a isa AbstractArray && return isequal(a, b)
    a isa Base.RefValue && return same_state(a[], b[])
    (isstructtype(typeof(a)) && nfields(a) > 0) || return isequal(a, b)
    return all(i -> same_state(getfield(a, i), getfield(b, i)), 1:nfields(a))
end

# Knots clustered toward the left do not look linear, so lookups use the stored
# interval: n = 12 searches with `KIND_LINEAR_SCAN`, n = 40 with `KIND_BRACKET_GALLOP`.
clustered_knots(n) = collect(100 .* range(0, 1; length = n) .^ 3)

# Deterministic, unsorted points spread over `[lo, hi]`.
spread(lo, hi, m) = lo .+ (hi - lo) .* mod.(0.6180339887 .* (1:m), 1)

function irregular_knot_interpolants(n)
    t = clustered_knots(n)
    u = sin.(t ./ 7) .+ 2   # positive, so the integral is invertible
    du = cos.(t ./ 7) ./ 7
    ddu = -sin.(t ./ 7) ./ 49
    umat = vcat(u', 2 .* u')
    uvv = [[x, 2x] for x in u]
    return [
        "LinearInterpolation" => kw -> LinearInterpolation(u, t; kw...),
        "LinearInterpolation (cache_parameters)" =>
            kw -> LinearInterpolation(u, t; cache_parameters = true, kw...),
        "LinearInterpolation (matrix)" => kw -> LinearInterpolation(umat, t; kw...),
        "LinearInterpolation (vector of vectors)" => kw -> LinearInterpolation(uvv, t; kw...),
        "QuadraticInterpolation" => kw -> QuadraticInterpolation(u, t; kw...),
        "QuadraticInterpolation (:Backward)" =>
            kw -> QuadraticInterpolation(u, t, :Backward; kw...),
        "QuadraticInterpolation (matrix)" => kw -> QuadraticInterpolation(umat, t; kw...),
        "LagrangeInterpolation" => kw -> LagrangeInterpolation(u, t; kw...),
        "AkimaInterpolation" => kw -> AkimaInterpolation(u, t; kw...),
        "AkimaInterpolation (matrix)" => kw -> AkimaInterpolation(umat, t; kw...),
        "ConstantInterpolation" => kw -> ConstantInterpolation(u, t; kw...),
        "ConstantInterpolation (dir = :right)" =>
            kw -> ConstantInterpolation(u, t; dir = :right, kw...),
        "ConstantInterpolation (matrix)" => kw -> ConstantInterpolation(umat, t; kw...),
        "SmoothedConstantInterpolation" => kw -> SmoothedConstantInterpolation(u, t; kw...),
        "QuadraticSpline" => kw -> QuadraticSpline(u, t; kw...),
        "QuadraticSpline (matrix)" => kw -> QuadraticSpline(umat, t; kw...),
        "CubicSpline" => kw -> CubicSpline(u, t; kw...),
        "CubicSpline (matrix)" => kw -> CubicSpline(umat, t; kw...),
        "CubicSpline (vector of vectors)" => kw -> CubicSpline(uvv, t; kw...),
        "BSplineInterpolation" => kw -> BSplineInterpolation(u, t, 3, :Average; kw...),
        "BSplineApprox" => kw -> BSplineApprox(u, t, 3, n ÷ 2, :Average; kw...),
        "CubicHermiteSpline" => kw -> CubicHermiteSpline(du, u, t; kw...),
        "CubicHermiteSpline (matrix)" =>
            kw -> CubicHermiteSpline(vcat(du', du'), umat, t; kw...),
        "PCHIPInterpolation" => kw -> PCHIPInterpolation(u, t; kw...),
        "QuinticHermiteSpline" => kw -> QuinticHermiteSpline(ddu, du, u, t; kw...),
        "SmoothArcLengthInterpolation" => kw -> SmoothArcLengthInterpolation(
            [0.3 -1.5 3.1; -0.2 0.2 -1.5; 0.0 0.0 0.0]'; m = 25, kw...
        ),
        "invert_integral(LinearInterpolation)" =>
            kw -> invert_integral(LinearInterpolation(u, t; kw...)),
        "invert_integral(ConstantInterpolation)" =>
            kw -> invert_integral(ConstantInterpolation(u, t; kw...)),
    ]
end

# Every evaluation entry point: values, first and second derivatives, vectorized and
# in-place evaluation and, where defined, integrals (over in-range points).
function evaluate_all(A, xs, xs_in)
    out = Any[A.(xs), derivative.(Ref(A), xs), derivative.(Ref(A), xs, 2), A(xs)]
    buf = similar(out[end])
    push!(out, A(buf, xs))
    if !(A isa Union{LagrangeInterpolation, DataInterpolations.AbstractIntegralInverseInterpolation})
        push!(out, [integral(A, first(A.t), x) for x in xs_in])
    end
    return out
end

@testset "cache_search = false: evaluation writes nothing (issue #618)" begin
    modes = (
        ExtrapolationType.Extension, ExtrapolationType.Periodic,
        ExtrapolationType.Reflective, ExtrapolationType.Linear, ExtrapolationType.Constant,
    )
    @testset "$name, n = $n, $mode" for n in (12, 40),
            (name, build) in irregular_knot_interpolants(n), mode in modes
        name == "SmoothedConstantInterpolation" && mode == ExtrapolationType.Linear && continue
        A = build((; cache_search = false, extrapolation = mode))
        @test !A.cache_search
        @test !A.iguesser.linear_lookup   # the stored-interval branch is the one under test
        lo, hi = first(A.t), last(A.t)
        xs_in = [A.t; (A.t[1:(end - 1)] .+ A.t[2:end]) ./ 2; spread(lo, hi, 200)]
        xs = [xs_in; lo - 5; hi + 5]
        # Far outside its range, a linear integral inverse has no real preimage (`sqrt` of
        # a negative number), with or without `cache_search`.
        A isa DataInterpolations.AbstractIntegralInverseInterpolation ||
            append!(xs, (lo - (hi - lo) / 3, hi + (hi - lo) / 3))
        A.iguesser.idx_prev[] = 0
        snapshot = deepcopy(A)
        results = evaluate_all(A, xs, xs_in)
        @test same_state(A, snapshot)
        # The interval found does not depend on the hint, so the default agrees exactly.
        @test isequal(results, evaluate_all(build((; extrapolation = mode)), xs, xs_in))
    end
end

@testset "unhinted lookup finds the same interval as the hinted one (issue #618)" begin
    sides = (
        (:last, (;)), (:last, (; lb = 1, ub_shift = 0)),
        (:first, (; idx_shift = -1)), (:first, (; lb = 1, ub_shift = 0)),
    )
    # Repeated knots encode jumps in `LinearInterpolation` (#610).
    repeated = [0.0, 1.0, 1.0, 2.0, 5.0, 5.0, 9.0, 20.0, 40.0, 41.0, 80.0, 100.0]
    @testset "$(length(t)) knots" for t in (clustered_knots(12), clustered_knots(40), repeated)
        u = sin.(t)
        A = LinearInterpolation(u, t; cache_search = false)
        B = LinearInterpolation(u, t)
        @test !A.iguesser.linear_lookup
        xs = [t; nextfloat.(t); prevfloat.(t); spread(first(t), last(t), 200); -1.0; 101.0]
        bad = []
        for x in xs, (side, kw) in sides
            lb = get(kw, :lb, 1)
            ub_shift = get(kw, :ub_shift, -1)
            idx_shift = get(kw, :idx_shift, 0)
            raw = side === :last ? searchsortedlast(t, x) : searchsortedfirst(t, x)
            ref = clamp(raw + idx_shift, lb, length(t) + ub_shift)
            got = (
                DataInterpolations.get_idx(A, x, A.iguesser; side, kw...),
                DataInterpolations.get_idx(B, x, B.iguesser; side, kw...),
            )
            got == (ref, ref) || push!(bad, (; x, side, kw, ref, got))
        end
        # Non-finite queries: both modes agree (whatever FindFirstFunctions defines).
        for x in (NaN, Inf, -Inf), (side, kw) in sides
            got = (
                DataInterpolations.get_idx(A, x, A.iguesser; side, kw...),
                DataInterpolations.get_idx(B, x, B.iguesser; side, kw...),
            )
            got[1] == got[2] || push!(bad, (; x, side, kw, got))
        end
        @test isempty(bad)
    end
end

@testset "default cache_search keeps the previous interval (issue #618)" begin
    t = clustered_knots(40)
    A = CubicSpline(sin.(t), t)
    @test A.cache_search
    @test !A.iguesser.linear_lookup
    A.iguesser.idx_prev[] = 1
    A((t[30] + t[31]) / 2)
    @test A.iguesser.idx_prev[] == 30
end

@testset "cache_search = false survives push! and append! (issue #618)" begin
    t = clustered_knots(40)
    u = sin.(t ./ 7) .+ 2
    @testset "$(nameof(typeof(A)))" for A in (
            LinearInterpolation(copy(u), copy(t); cache_search = false),
            QuadraticInterpolation(copy(u), copy(t); cache_search = false),
            ConstantInterpolation(copy(u), copy(t); cache_search = false),
        )
        push!(A, 2.5, 101.0)
        append!(A, [2.0, 3.0], [102.0, 104.0])
        @test !A.cache_search
        A.iguesser.idx_prev[] = 0
        snapshot = deepcopy(A)
        for x in spread(first(A.t), last(A.t), 200)
            A(x)
            derivative(A, x)
        end
        @test same_state(A, snapshot)
    end
end

if Threads.nthreads() > 1
    @testset "shared cache_search = false interpolant across threads (issue #618)" begin
        @testset "$name" for (name, build) in irregular_knot_interpolants(40)
            A = build((; cache_search = false))
            xs = spread(first(A.t), last(A.t), 61)
            ref = A.(xs)
            dref = derivative.(Ref(A), xs)
            A.iguesser.idx_prev[] = 0
            snapshot = deepcopy(A)
            wrong = Threads.Atomic{Int}(0)
            Threads.@threads for _ in 1:2_000
                for (j, x) in enumerate(xs)
                    isequal(A(x), ref[j]) && isequal(derivative(A, x), dref[j]) ||
                        Threads.atomic_add!(wrong, 1)
                end
            end
            @test wrong[] == 0
            @test same_state(A, snapshot)
        end
    end
else
    @info "issue #618: shared-interpolant thread test skipped; start Julia with more " *
        "than one thread (e.g. `julia -t auto`) to exercise it."
end
