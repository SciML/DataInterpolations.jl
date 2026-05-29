using DataInterpolations
using DataInterpolations: derivative
using Test
using Base.Threads

# Regression tests for issue #532: a `BSplineInterpolation`, `BSplineApprox` or
# (default, uncached) `QuadraticSpline` kept a scratch buffer in its `sc` field
# and overwrote it in place on every evaluation. Because a "read" (evaluation)
# mutated shared state, evaluating one shared interpolant from several threads
# silently returned wrong values. Evaluation (and differentiation) must be
# reentrant: it must not write to any state stored in the interpolant.

function thread_safety_interpolants()
    u = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
    t = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    umat = [
        0.0 1.0 0.0 1.0 0.0 1.0
        1.0 0.0 1.0 0.0 1.0 0.0
    ]            # 2 × length(t): last axis indexes t
    return [
        ("BSplineInterpolation", BSplineInterpolation(u, t, 3, :Uniform, :Average)),
        ("BSplineApprox", BSplineApprox(u, t, 3, 4, :Uniform, :Average)),
        ("BSplineInterpolation (array)", BSplineInterpolation(umat, t, 3, :Uniform, :Average)),
        ("BSplineApprox (array)", BSplineApprox(umat, t, 3, 4, :Uniform, :Average)),
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
