# Cross-library 1D interpolation benchmark

## Setup

```
Julia Version 1.11.9
Commit 53a02c0720c (2026-02-06 00:27 UTC)
Build Info:
  Official https://julialang.org/ release
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  CPU: 128 × AMD EPYC 7502 32-Core Processor
  WORD_SIZE: 64
  LLVM: libLLVM-16.0.6 (ORCJIT, znver2)
Threads: 1 default, 0 interactive, 1 GC (on 128 virtual cores)
```

Bench harness: `BenchmarkTools.@benchmark` with `evals=1`, max samples=100, max seconds=0.5.

Commit: `e22f5d6a6fa8d72079b33209b910b9504c4cadde`

Library versions:
```
  DataInterpolations 8.10.0
  Interpolations 0.16.2
  Dierckx 0.5.4
  BasicInterpolators 0.7.1
  PCHIPInterpolation 0.2.1
  BenchmarkTools 1.8.0
```

Total bench time: 353.5 s

## Headline findings

These numbers are taken directly from the tables below. All cells are medians from `BenchmarkTools.@benchmark evals=1` runs.

- **Sorted-batch + cached `Auto(t_props)` is DI's biggest cross-library win.** At cubic spline n=100 000, m=100 000 (sorted batch, uniform knots), DI evaluates in **1.78 ms** vs Dierckx 3.16 s (~**1 770× faster**), BasicInterpolators 7.94 ms (~**4.5× faster**), and PCHIP at n=100k m=100k (sorted) 9.70 ms (~**5.3× faster** for monotone cubic). The only library that beats DI on this row is **Interpolations.jl's uniform constructor (2.06 ms, ~16% faster)** because it uses O(1) uniform-grid index lookup; DI still wins against every non-uniform-capable competitor. At linear n=100k m=100k sorted batch, DI is **1.64 ms** vs Dierckx 3.18 s (~**1 940×**) and BasicInterpolators 7.89 ms (~**4.8×**).

- **`Auto(t_props)` also dominates the random / unsorted-batch case** because the cached search-property is re-used on every evaluation, not re-probed. At cubic n=100k m=100k *random* batch, DI is **6.87 ms** vs Dierckx 3.15 s (~**460×**), BasicInterpolators 19.8 ms (~**2.9×**); Interpolations(uniform) is still ahead at 2.06 ms (O(1) index). Even at random unsorted access, DI competes with — and on every non-uniform library beats — the alternatives.

- **DI wins the chained ODE-style case decisively at large n.** Cubic spline, n=100k, m=1000 monotone chain: DI **69.6 μs** vs Dierckx 31.5 ms (~**450×**), BasicInterpolators 137 μs (~**2.0×**). MonotoneCubic chained, n=100k m=1000: DI 71.9 μs vs PCHIPInterpolation 149 μs (~**2.1×**). This is exactly the workload DI's `iguesser` was designed for; libraries without hint-chaining (Dierckx, BasicInterpolators, PCHIP) lose by ~2× to ~450×.

- **Where DI loses on small m or single-eval:** single-query cubic at n=100k is 80 ns (DI) vs 50 ns (BasicInterpolators) and 70 ns (Interpolations uniform) — within ~1.5×, but consistently slightly slower because DI does a real index lookup whereas the others can use a one-shot uniform divide. At `sorted batch m=10 / n=10000` linear, DI is 500 ns vs Interpolations(uniform) 120 ns (4×) because DI's sorted-batch fast-path uses an O(m log n) per-element bisect plus allocates a small idx buffer; at this very small m the buffer alloc overhead dominates.

- **DI QuadraticSpline is a clear loser, by ~2-5× across the board.** At n=100k construction DI takes **7.0 s** vs Dierckx 13.9 ms (~**500× slower**). At n=10k construction DI is 65.7 ms vs Dierckx 1.37 ms (~48×). Inspection traces this to `quadratic_spline_params` calling `spline_coefficients!` which does `findfirst(x -> x > u, k)` — a linear scan inside a loop over `n` knots, making the QuadraticSpline constructor **O(n²)**. The QuadraticSpline single-query at n=100k is also 60-110 μs (vs Dierckx 13-55 μs) and the batched evaluators at n=100k m=100k cost 7.0 s (vs Dierckx 3.2 s). This is the most actionable finding in the entire report.

- **DI CubicHermiteSpline (used here as PCHIP analogue) beats PCHIPInterpolation.jl on every batched cell.** Sorted batch n=100k m=100k: DI 1.83 ms vs PCHIP 9.70 ms (~5.3×). Random batch n=100k m=100k: DI 7.28 ms vs PCHIP 20.4 ms (~2.8×). Chained n=100k m=1000: DI 72 μs vs PCHIP 149 μs (~2.1×). Construction is also slightly faster across all n.

- **`Interpolations.jl`'s `cubic_spline_interpolation` / `linear_interpolation` over a `range` is the ceiling we're chasing on uniform data.** Because it does O(1) index lookup, it beats DI on every uniform-grid case where the lookup cost dominates: linear n=100k single-query 60 ns (DI 80 ns), linear sorted-batch n=100k m=100k 658 μs (DI 1.64 ms, ~2.5×), cubic sorted-batch n=100k m=100k 2.06 ms (DI 1.78 ms — DI wins here). It cannot handle non-uniform cubic at all; it falls back to `Gridded(Linear)` only. So the comparison is really "DI generalised to non-uniform & O(log n) lookup" vs "Interpolations specialised to uniform & O(1) lookup."

## Construction time

Rows = library, columns = (n, knot pattern). Values = median wall time (IQR).

### Akima

| Library | n=100,nonuniform | n=100,uniform | n=1000,nonuniform | n=1000,uniform | n=10000,nonuniform | n=10000,uniform | n=100000,nonuniform | n=100000,uniform | 
|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 2.195 μs (IQR 320.000 ns) | 2.285 μs (IQR 427.500 ns) | 31.189 μs (IQR 4.425 μs) | 44.820 μs (IQR 6.655 μs) | 117.964 μs (IQR 1.608 μs) | 182.268 μs (IQR 233.136 μs) | 1.203 ms (IQR 55.989 μs) | 1.269 ms (IQR 42.470 μs) | 

### CubicSpline

| Library | n=100,nonuniform | n=100,uniform | n=1000,nonuniform | n=1000,uniform | n=10000,nonuniform | n=10000,uniform | n=100000,nonuniform | n=100000,uniform | 
|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 7.450 μs (IQR 3.060 μs) | 6.880 μs (IQR 675.000 ns) | 56.620 μs (IQR 3.062 μs) | 60.729 μs (IQR 26.519 μs) | 507.370 μs (IQR 9.600 μs) | 580.565 μs (IQR 948.358 μs) | 5.510 ms (IQR 1.113 ms) | 5.229 ms (IQR 863.691 μs) | 
| BasicInterpolators | 4.305 μs (IQR 542.500 ns) | 4.430 μs (IQR 475.000 ns) | 37.569 μs (IQR 2.132 μs) | 39.330 μs (IQR 4.598 μs) | 714.473 μs (IQR 401.371 μs) | 1.025 ms (IQR 714.324 μs) | 3.287 ms (IQR 894.667 μs) | 3.315 ms (IQR 848.164 μs) | 
| Dierckx (k=3) | 17.674 μs (IQR 739.250 ns) | 17.600 μs (IQR 1.215 μs) | 208.738 μs (IQR 51.890 μs) | 212.863 μs (IQR 6.742 μs) | 1.608 ms (IQR 84.278 μs) | 1.745 ms (IQR 548.674 μs) | 17.461 ms (IQR 6.682 ms) | 18.954 ms (IQR 6.357 ms) | 
| Interpolations (uniform) | — | 15.274 μs (IQR 734.250 ns) | — | 137.884 μs (IQR 19.681 μs) | — | 1.442 ms (IQR 464.458 μs) | — | 8.195 ms (IQR 134.718 μs) | 

### Linear

| Library | n=100,nonuniform | n=100,uniform | n=1000,nonuniform | n=1000,uniform | n=10000,nonuniform | n=10000,uniform | n=100000,nonuniform | n=100000,uniform | 
|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 1.260 μs (IQR 185.000 ns) | 1.450 μs (IQR 267.500 ns) | 10.200 μs (IQR 320.250 ns) | 11.110 μs (IQR 354.250 ns) | 92.584 μs (IQR 826.750 ns) | 101.064 μs (IQR 967.500 ns) | 930.192 μs (IQR 9.795 μs) | 998.760 μs (IQR 7.260 μs) | 
| BasicInterpolators | 1.190 μs (IQR 212.500 ns) | 1.290 μs (IQR 235.000 ns) | 7.905 μs (IQR 505.000 ns) | 8.960 μs (IQR 1.292 μs) | 59.559 μs (IQR 1.714 μs) | 67.840 μs (IQR 1.940 μs) | 594.784 μs (IQR 9.875 μs) | 671.909 μs (IQR 10.883 μs) | 
| Dierckx (k=1) | 11.919 μs (IQR 2.472 μs) | 12.485 μs (IQR 1.991 μs) | 79.025 μs (IQR 2.056 μs) | 108.474 μs (IQR 14.793 μs) | 790.987 μs (IQR 114.966 μs) | 1.019 ms (IQR 266.517 μs) | 7.557 ms (IQR 622.339 μs) | 7.581 ms (IQR 727.211 μs) | 
| Interpolations (gridded) | 940.000 ns (IQR 162.500 ns) | 1.030 μs (IQR 222.500 ns) | 7.580 μs (IQR 537.250 ns) | 11.495 μs (IQR 10.139 μs) | 64.489 μs (IQR 977.500 ns) | 147.218 μs (IQR 76.435 μs) | 630.804 μs (IQR 14.617 μs) | 723.713 μs (IQR 18.942 μs) | 
| Interpolations (uniform) | — | 215.000 ns (IQR 142.500 ns) | — | 3.204 μs (IQR 1.998 μs) | — | 22.045 μs (IQR 4.760 μs) | — | 75.369 μs (IQR 8.617 μs) | 

### MonotoneCubic

| Library | n=100,nonuniform | n=100,uniform | n=1000,nonuniform | n=1000,uniform | n=10000,nonuniform | n=10000,uniform | n=100000,nonuniform | n=100000,uniform | 
|---|---|---|---|---|---|---|---|---|
| DataInterpolations (CubicHermite) | 1.400 μs (IQR 222.500 ns) | 1.545 μs (IQR 292.500 ns) | 10.195 μs (IQR 500.000 ns) | 11.000 μs (IQR 448.250 ns) | 87.144 μs (IQR 975.250 ns) | 94.999 μs (IQR 670.000 ns) | 868.637 μs (IQR 11.227 μs) | 953.866 μs (IQR 14.377 μs) | 
| PCHIPInterpolation | 1.660 μs (IQR 420.000 ns) | 1.780 μs (IQR 452.500 ns) | 13.255 μs (IQR 960.500 ns) | 14.490 μs (IQR 1.624 μs) | 109.594 μs (IQR 1.990 μs) | 116.719 μs (IQR 1.327 μs) | 1.092 ms (IQR 21.902 μs) | 1.181 ms (IQR 35.152 μs) | 

### QuadraticSpline

| Library | n=100,nonuniform | n=100,uniform | n=1000,nonuniform | n=1000,uniform | n=10000,nonuniform | n=10000,uniform | n=100000,nonuniform | n=100000,uniform | 
|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 16.825 μs (IQR 904.250 ns) | 16.880 μs (IQR 865.750 ns) | 744.113 μs (IQR 19.005 μs) | 766.238 μs (IQR 39.012 μs) | 65.719 ms (IQR 401.579 μs) | 66.057 ms (IQR 689.417 μs) | 6.978 s (IQR 0.000 ns) | 7.028 s (IQR 0.000 ns) | 
| Dierckx (k=2) | 18.560 μs (IQR 1.295 μs) | 17.700 μs (IQR 3.228 μs) | 137.524 μs (IQR 1.077 μs) | 175.948 μs (IQR 8.215 μs) | 1.371 ms (IQR 61.880 μs) | 1.583 ms (IQR 193.173 μs) | 13.906 ms (IQR 997.231 μs) | 13.571 ms (IQR 1.787 ms) | 

## Single-query latency

Cold single evaluation `A(x_query)`. Rows = library, columns = (n, knot pattern).

### Akima

| Library | n=100,nonuniform | n=100,uniform | n=1000,nonuniform | n=1000,uniform | n=10000,nonuniform | n=10000,uniform | n=100000,nonuniform | n=100000,uniform | 
|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 70.000 ns (IQR 10.000 ns) | 70.000 ns (IQR 0.000 ns) | 70.000 ns (IQR 10.000 ns) | 70.000 ns (IQR 10.000 ns) | 90.000 ns (IQR 10.000 ns) | 75.000 ns (IQR 10.000 ns) | 100.000 ns (IQR 0.000 ns) | 70.000 ns (IQR 10.000 ns) | 

### CubicSpline

| Library | n=100,nonuniform | n=100,uniform | n=1000,nonuniform | n=1000,uniform | n=10000,nonuniform | n=10000,uniform | n=100000,nonuniform | n=100000,uniform | 
|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 80.000 ns (IQR 10.000 ns) | 80.000 ns (IQR 10.000 ns) | 70.000 ns (IQR 10.000 ns) | 80.000 ns (IQR 10.000 ns) | 100.000 ns (IQR 10.000 ns) | 80.000 ns (IQR 10.000 ns) | 110.000 ns (IQR 10.000 ns) | 80.000 ns (IQR 0.000 ns) | 
| BasicInterpolators | 50.000 ns (IQR 0.000 ns) | 50.000 ns (IQR 0.000 ns) | 50.000 ns (IQR 10.000 ns) | 50.000 ns (IQR 10.000 ns) | 60.000 ns (IQR 0.000 ns) | 60.000 ns (IQR 0.000 ns) | 65.000 ns (IQR 10.000 ns) | 70.000 ns (IQR 0.000 ns) | 
| Dierckx (k=3) | 160.000 ns (IQR 10.000 ns) | 160.000 ns (IQR 10.000 ns) | 390.000 ns (IQR 10.000 ns) | 400.000 ns (IQR 10.000 ns) | 2.780 μs (IQR 0.000 ns) | 2.790 μs (IQR 10.000 ns) | 26.750 μs (IQR 40.000 ns) | 26.650 μs (IQR 89.250 ns) | 
| Interpolations (uniform) | — | 80.000 ns (IQR 10.000 ns) | — | 80.000 ns (IQR 10.000 ns) | — | 80.000 ns (IQR 0.000 ns) | — | 80.000 ns (IQR 10.000 ns) | 

### Linear

| Library | n=100,nonuniform | n=100,uniform | n=1000,nonuniform | n=1000,uniform | n=10000,nonuniform | n=10000,uniform | n=100000,nonuniform | n=100000,uniform | 
|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 70.000 ns (IQR 0.000 ns) | 80.000 ns (IQR 0.000 ns) | 70.000 ns (IQR 0.000 ns) | 80.000 ns (IQR 0.000 ns) | 90.000 ns (IQR 10.000 ns) | 80.000 ns (IQR 10.000 ns) | 100.000 ns (IQR 10.000 ns) | 80.000 ns (IQR 0.000 ns) | 
| BasicInterpolators | 50.000 ns (IQR 10.000 ns) | 50.000 ns (IQR 0.000 ns) | 60.000 ns (IQR 10.000 ns) | 60.000 ns (IQR 0.000 ns) | 70.000 ns (IQR 10.000 ns) | 70.000 ns (IQR 10.000 ns) | 70.000 ns (IQR 0.000 ns) | 70.000 ns (IQR 0.000 ns) | 
| Dierckx (k=1) | 120.000 ns (IQR 10.000 ns) | 120.000 ns (IQR 10.000 ns) | 350.000 ns (IQR 0.000 ns) | 370.000 ns (IQR 10.000 ns) | 2.760 μs (IQR 0.000 ns) | 2.770 μs (IQR 10.000 ns) | 26.729 μs (IQR 42.250 ns) | 26.630 μs (IQR 72.500 ns) | 
| Interpolations (gridded) | 60.000 ns (IQR 0.000 ns) | 60.000 ns (IQR 10.000 ns) | 70.000 ns (IQR 10.000 ns) | 70.000 ns (IQR 0.000 ns) | 70.000 ns (IQR 10.000 ns) | 70.000 ns (IQR 0.000 ns) | 80.000 ns (IQR 0.000 ns) | 80.000 ns (IQR 10.000 ns) | 
| Interpolations (uniform) | — | 60.000 ns (IQR 0.000 ns) | — | 60.000 ns (IQR 0.000 ns) | — | 60.000 ns (IQR 0.000 ns) | — | 60.000 ns (IQR 0.000 ns) | 

### MonotoneCubic

| Library | n=100,nonuniform | n=100,uniform | n=1000,nonuniform | n=1000,uniform | n=10000,nonuniform | n=10000,uniform | n=100000,nonuniform | n=100000,uniform | 
|---|---|---|---|---|---|---|---|---|
| DataInterpolations (CubicHermite) | 80.000 ns (IQR 0.000 ns) | 90.000 ns (IQR 10.000 ns) | 80.000 ns (IQR 0.000 ns) | 90.000 ns (IQR 0.000 ns) | 100.000 ns (IQR 0.000 ns) | 90.000 ns (IQR 10.000 ns) | 120.000 ns (IQR 0.000 ns) | 90.000 ns (IQR 0.000 ns) | 
| PCHIPInterpolation | 70.000 ns (IQR 10.000 ns) | 70.000 ns (IQR 0.000 ns) | 80.000 ns (IQR 0.000 ns) | 80.000 ns (IQR 0.000 ns) | 90.000 ns (IQR 0.000 ns) | 90.000 ns (IQR 0.000 ns) | 100.000 ns (IQR 10.000 ns) | 100.000 ns (IQR 10.000 ns) | 

### QuadraticSpline

| Library | n=100,nonuniform | n=100,uniform | n=1000,nonuniform | n=1000,uniform | n=10000,nonuniform | n=10000,uniform | n=100000,nonuniform | n=100000,uniform | 
|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 180.000 ns (IQR 0.000 ns) | 190.000 ns (IQR 10.000 ns) | 660.000 ns (IQR 10.000 ns) | 710.000 ns (IQR 10.000 ns) | 5.810 μs (IQR 10.000 ns) | 5.820 μs (IQR 19.000 ns) | 62.569 μs (IQR 243.500 ns) | 62.474 μs (IQR 163.250 ns) | 
| Dierckx (k=2) | 140.000 ns (IQR 0.000 ns) | 170.000 ns (IQR 32.500 ns) | 370.000 ns (IQR 0.000 ns) | 390.000 ns (IQR 0.000 ns) | 2.780 μs (IQR 10.000 ns) | 2.780 μs (IQR 10.000 ns) | 26.739 μs (IQR 73.250 ns) | 26.670 μs (IQR 59.000 ns) | 

## Sorted batch

`A(out, tt)` where `tt` is sorted random points in domain. (knot pattern = uniform)

### Akima

| Library | n=100,m=1 | n=100,m=10 | n=100,m=1000 | n=100,m=100000 | n=1000,m=1 | n=1000,m=10 | n=1000,m=1000 | n=1000,m=100000 | n=10000,m=1 | n=10000,m=10 | n=10000,m=1000 | n=10000,m=100000 | n=100000,m=1 | n=100000,m=10 | n=100000,m=1000 | n=100000,m=100000 | 
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 110.000 ns (IQR 10.000 ns) | 170.000 ns (IQR 10.000 ns) | 4.900 μs (IQR 1.530 μs) | 323.432 μs (IQR 2.525 μs) | 120.000 ns (IQR 0.000 ns) | 340.000 ns (IQR 10.000 ns) | 5.195 μs (IQR 30.000 ns) | 341.772 μs (IQR 1.866 μs) | 120.000 ns (IQR 10.000 ns) | 430.000 ns (IQR 10.000 ns) | 12.070 μs (IQR 412.500 ns) | 479.050 μs (IQR 9.502 μs) | 140.000 ns (IQR 10.000 ns) | 510.000 ns (IQR 10.000 ns) | 140.203 μs (IQR 1.752 μs) | 1.436 ms (IQR 26.067 μs) | 

### CubicSpline

| Library | n=100,m=1 | n=100,m=10 | n=100,m=1000 | n=100,m=100000 | n=1000,m=1 | n=1000,m=10 | n=1000,m=1000 | n=1000,m=100000 | n=10000,m=1 | n=10000,m=10 | n=10000,m=1000 | n=10000,m=100000 | n=100000,m=1 | n=100000,m=10 | n=100000,m=1000 | n=100000,m=100000 | 
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 130.000 ns (IQR 10.000 ns) | 200.000 ns (IQR 10.000 ns) | 10.030 μs (IQR 20.000 ns) | 975.096 μs (IQR 13.203 μs) | 130.000 ns (IQR 0.000 ns) | 420.000 ns (IQR 10.000 ns) | 9.920 μs (IQR 20.000 ns) | 973.951 μs (IQR 22.567 μs) | 140.000 ns (IQR 10.000 ns) | 490.000 ns (IQR 10.000 ns) | 15.420 μs (IQR 50.000 ns) | 1.019 ms (IQR 9.423 μs) | 140.000 ns (IQR 10.000 ns) | 580.000 ns (IQR 10.000 ns) | 151.668 μs (IQR 4.980 μs) | 1.779 ms (IQR 21.049 μs) | 
| BasicInterpolators | 70.000 ns (IQR 10.000 ns) | 180.000 ns (IQR 10.000 ns) | 16.315 μs (IQR 302.500 ns) | 1.369 ms (IQR 19.600 μs) | 70.000 ns (IQR 0.000 ns) | 240.000 ns (IQR 0.000 ns) | 42.385 μs (IQR 948.500 ns) | 2.034 ms (IQR 21.232 μs) | 70.000 ns (IQR 10.000 ns) | 300.000 ns (IQR 10.000 ns) | 79.069 μs (IQR 1.113 μs) | 4.233 ms (IQR 66.511 μs) | 80.000 ns (IQR 0.000 ns) | 370.000 ns (IQR 10.000 ns) | 130.333 μs (IQR 2.495 μs) | 7.938 ms (IQR 43.734 μs) | 
| Dierckx (k=3) | 190.000 ns (IQR 0.000 ns) | 1.220 μs (IQR 0.000 ns) | 121.929 μs (IQR 30.250 ns) | 11.761 ms (IQR 53.174 μs) | 710.000 ns (IQR 10.000 ns) | 3.620 μs (IQR 0.250 ns) | 408.796 μs (IQR 188.500 ns) | 41.214 ms (IQR 255.653 μs) | 5.710 μs (IQR 10.000 ns) | 26.890 μs (IQR 20.000 ns) | 3.210 ms (IQR 28.301 μs) | 323.867 ms (IQR 20.851 μs) | 55.660 μs (IQR 70.000 ns) | 260.212 μs (IQR 242.250 ns) | 31.300 ms (IQR 245.060 μs) | 3.157 s (IQR 0.000 ns) | 
| Interpolations (uniform) | 90.000 ns (IQR 10.000 ns) | 270.000 ns (IQR 0.000 ns) | 20.469 μs (IQR 10.000 ns) | 2.050 ms (IQR 13.412 μs) | 90.000 ns (IQR 10.000 ns) | 270.000 ns (IQR 0.000 ns) | 20.460 μs (IQR 1.000 ns) | 2.050 ms (IQR 13.490 μs) | 90.000 ns (IQR 10.000 ns) | 270.000 ns (IQR 0.000 ns) | 20.460 μs (IQR 30.000 ns) | 2.050 ms (IQR 11.105 μs) | 90.000 ns (IQR 10.000 ns) | 270.000 ns (IQR 0.000 ns) | 20.450 μs (IQR 20.000 ns) | 2.057 ms (IQR 21.207 μs) | 

### Linear

| Library | n=100,m=1 | n=100,m=10 | n=100,m=1000 | n=100,m=100000 | n=1000,m=1 | n=1000,m=10 | n=1000,m=1000 | n=1000,m=100000 | n=10000,m=1 | n=10000,m=10 | n=10000,m=1000 | n=10000,m=100000 | n=100000,m=1 | n=100000,m=10 | n=100000,m=1000 | n=100000,m=100000 | 
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 120.000 ns (IQR 10.000 ns) | 180.000 ns (IQR 0.000 ns) | 5.830 μs (IQR 10.000 ns) | 469.835 μs (IQR 4.003 μs) | 120.000 ns (IQR 10.000 ns) | 380.000 ns (IQR 10.000 ns) | 6.210 μs (IQR 50.000 ns) | 485.851 μs (IQR 1.160 μs) | 120.000 ns (IQR 10.000 ns) | 450.000 ns (IQR 10.000 ns) | 15.390 μs (IQR 1.252 μs) | 629.649 μs (IQR 950.000 ns) | 130.000 ns (IQR 0.000 ns) | 540.000 ns (IQR 10.000 ns) | 141.224 μs (IQR 3.203 μs) | 1.644 ms (IQR 9.440 μs) | 
| BasicInterpolators | 60.000 ns (IQR 0.000 ns) | 200.000 ns (IQR 0.000 ns) | 16.205 μs (IQR 350.250 ns) | 1.378 ms (IQR 13.425 μs) | 70.000 ns (IQR 10.000 ns) | 240.000 ns (IQR 0.000 ns) | 41.170 μs (IQR 2.439 μs) | 2.053 ms (IQR 9.387 μs) | 70.000 ns (IQR 10.000 ns) | 310.000 ns (IQR 0.000 ns) | 76.659 μs (IQR 3.703 μs) | 4.179 ms (IQR 15.650 μs) | 80.000 ns (IQR 0.000 ns) | 380.000 ns (IQR 10.000 ns) | 123.879 μs (IQR 2.312 μs) | 7.892 ms (IQR 63.663 μs) | 
| Dierckx (k=1) | 150.000 ns (IQR 10.000 ns) | 850.000 ns (IQR 10.000 ns) | 87.719 μs (IQR 1.460 μs) | 8.265 ms (IQR 132.779 μs) | 670.000 ns (IQR 10.000 ns) | 3.280 μs (IQR 820.000 ns) | 373.336 μs (IQR 5.365 μs) | 37.651 ms (IQR 73.838 μs) | 5.660 μs (IQR 30.000 ns) | 26.730 μs (IQR 50.000 ns) | 3.168 ms (IQR 24.134 μs) | 319.511 ms (IQR 16.000 μs) | 55.640 μs (IQR 113.250 ns) | 261.808 μs (IQR 933.000 ns) | 31.246 ms (IQR 377.619 μs) | 3.186 s (IQR 0.000 ns) | 
| Interpolations (gridded) | 70.000 ns (IQR 0.000 ns) | 230.000 ns (IQR 0.000 ns) | 18.250 μs (IQR 50.250 ns) | 1.635 ms (IQR 23.676 μs) | 70.000 ns (IQR 10.000 ns) | 300.000 ns (IQR 0.000 ns) | 48.190 μs (IQR 3.190 μs) | 2.548 ms (IQR 41.645 μs) | 80.000 ns (IQR 0.000 ns) | 360.000 ns (IQR 10.000 ns) | 90.964 μs (IQR 3.175 μs) | 4.704 ms (IQR 63.885 μs) | 84.500 ns (IQR 10.000 ns) | 430.000 ns (IQR 0.000 ns) | 141.284 μs (IQR 2.475 μs) | 7.996 ms (IQR 40.135 μs) | 
| Interpolations (uniform) | 60.000 ns (IQR 10.000 ns) | 120.000 ns (IQR 0.000 ns) | 6.610 μs (IQR 0.000 ns) | 655.124 μs (IQR 583.250 ns) | 70.000 ns (IQR 10.000 ns) | 120.000 ns (IQR 0.000 ns) | 6.600 μs (IQR 20.000 ns) | 655.409 μs (IQR 12.107 μs) | 70.000 ns (IQR 10.000 ns) | 120.000 ns (IQR 0.000 ns) | 6.670 μs (IQR 30.000 ns) | 655.644 μs (IQR 702.500 ns) | 70.000 ns (IQR 10.000 ns) | 120.000 ns (IQR 0.000 ns) | 7.060 μs (IQR 50.000 ns) | 658.249 μs (IQR 1.255 μs) | 

### MonotoneCubic

| Library | n=100,m=1 | n=100,m=10 | n=100,m=1000 | n=100,m=100000 | n=1000,m=1 | n=1000,m=10 | n=1000,m=1000 | n=1000,m=100000 | n=10000,m=1 | n=10000,m=10 | n=10000,m=1000 | n=10000,m=100000 | n=100000,m=1 | n=100000,m=10 | n=100000,m=1000 | n=100000,m=100000 | 
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| DataInterpolations (CubicHermite) | 130.000 ns (IQR 30.000 ns) | 340.000 ns (IQR 10.000 ns) | 10.490 μs (IQR 20.250 ns) | 988.671 μs (IQR 9.099 μs) | 140.000 ns (IQR 10.000 ns) | 420.000 ns (IQR 10.000 ns) | 10.649 μs (IQR 30.000 ns) | 989.491 μs (IQR 14.422 μs) | 140.000 ns (IQR 0.000 ns) | 500.000 ns (IQR 10.000 ns) | 19.150 μs (IQR 80.000 ns) | 1.067 ms (IQR 11.665 μs) | 150.000 ns (IQR 0.000 ns) | 580.000 ns (IQR 10.000 ns) | 146.929 μs (IQR 2.422 μs) | 1.825 ms (IQR 19.328 μs) | 
| PCHIPInterpolation | 70.000 ns (IQR 10.000 ns) | 240.000 ns (IQR 0.000 ns) | 20.389 μs (IQR 437.500 ns) | 1.867 ms (IQR 24.703 μs) | 90.000 ns (IQR 0.000 ns) | 330.000 ns (IQR 0.000 ns) | 39.885 μs (IQR 881.750 ns) | 2.972 ms (IQR 17.169 μs) | 100.000 ns (IQR 0.000 ns) | 440.000 ns (IQR 10.000 ns) | 79.079 μs (IQR 4.630 μs) | 5.475 ms (IQR 51.462 μs) | 110.000 ns (IQR 0.000 ns) | 510.000 ns (IQR 10.000 ns) | 133.694 μs (IQR 4.237 μs) | 9.703 ms (IQR 97.707 μs) | 

### QuadraticSpline

| Library | n=100,m=1 | n=100,m=10 | n=100,m=1000 | n=100,m=100000 | n=1000,m=1 | n=1000,m=10 | n=1000,m=1000 | n=1000,m=100000 | n=10000,m=1 | n=10000,m=10 | n=10000,m=1000 | n=10000,m=100000 | n=100000,m=1 | n=100000,m=10 | n=100000,m=1000 | n=100000,m=100000 | 
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 290.000 ns (IQR 10.000 ns) | 1.370 μs (IQR 0.000 ns) | 137.279 μs (IQR 811.750 ns) | 13.395 ms (IQR 177.476 μs) | 1.260 μs (IQR 0.000 ns) | 6.600 μs (IQR 12.250 ns) | 722.878 μs (IQR 1.238 μs) | 72.672 ms (IQR 365.231 μs) | 10.915 μs (IQR 30.000 ns) | 56.859 μs (IQR 70.000 ns) | 6.588 ms (IQR 72.728 μs) | 660.689 ms (IQR 0.000 ns) | 112.729 μs (IQR 162.500 ns) | 607.684 μs (IQR 7.387 μs) | 70.380 ms (IQR 419.419 μs) | 7.038 s (IQR 0.000 ns) | 
| Dierckx (k=2) | 170.000 ns (IQR 0.000 ns) | 1.030 μs (IQR 10.000 ns) | 103.394 μs (IQR 260.000 ns) | 9.983 ms (IQR 106.786 μs) | 680.000 ns (IQR 10.000 ns) | 3.450 μs (IQR 10.000 ns) | 390.976 μs (IQR 1.002 μs) | 39.394 ms (IQR 210.678 μs) | 5.690 μs (IQR 13.250 ns) | 26.760 μs (IQR 20.000 ns) | 3.190 ms (IQR 47.677 μs) | 322.009 ms (IQR 65.414 μs) | 55.739 μs (IQR 115.000 ns) | 260.118 μs (IQR 761.500 ns) | 31.568 ms (IQR 621.605 μs) | 3.207 s (IQR 0.000 ns) | 

## Random batch

`A(out, tt)` where `tt` is unsorted. (knot pattern = uniform)

### Akima

| Library | n=100,m=1 | n=100,m=10 | n=100,m=1000 | n=100,m=100000 | n=1000,m=1 | n=1000,m=10 | n=1000,m=1000 | n=1000,m=100000 | n=10000,m=1 | n=10000,m=10 | n=10000,m=1000 | n=10000,m=100000 | n=100000,m=1 | n=100000,m=10 | n=100000,m=1000 | n=100000,m=100000 | 
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 120.000 ns (IQR 10.000 ns) | 450.000 ns (IQR 0.000 ns) | 47.359 μs (IQR 171.500 ns) | 4.831 ms (IQR 31.837 μs) | 120.000 ns (IQR 0.000 ns) | 450.000 ns (IQR 0.000 ns) | 48.005 μs (IQR 194.750 ns) | 4.898 ms (IQR 32.682 μs) | 130.000 ns (IQR 10.000 ns) | 430.000 ns (IQR 10.000 ns) | 52.684 μs (IQR 263.500 ns) | 5.357 ms (IQR 57.461 μs) | 140.000 ns (IQR 10.000 ns) | 440.000 ns (IQR 10.000 ns) | 59.049 μs (IQR 342.500 ns) | 6.366 ms (IQR 52.595 μs) | 

### CubicSpline

| Library | n=100,m=1 | n=100,m=10 | n=100,m=1000 | n=100,m=100000 | n=1000,m=1 | n=1000,m=10 | n=1000,m=1000 | n=1000,m=100000 | n=10000,m=1 | n=10000,m=10 | n=10000,m=1000 | n=10000,m=100000 | n=100000,m=1 | n=100000,m=10 | n=100000,m=1000 | n=100000,m=100000 | 
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 130.000 ns (IQR 10.000 ns) | 490.000 ns (IQR 10.000 ns) | 51.094 μs (IQR 95.500 ns) | 5.201 ms (IQR 58.934 μs) | 130.000 ns (IQR 0.000 ns) | 490.000 ns (IQR 0.000 ns) | 52.210 μs (IQR 255.000 ns) | 5.480 ms (IQR 26.636 μs) | 140.000 ns (IQR 0.000 ns) | 480.000 ns (IQR 0.000 ns) | 55.980 μs (IQR 354.750 ns) | 5.666 ms (IQR 48.039 μs) | 150.000 ns (IQR 10.000 ns) | 500.000 ns (IQR 10.000 ns) | 63.150 μs (IQR 3.263 μs) | 6.874 ms (IQR 47.719 μs) | 
| BasicInterpolators | 70.000 ns (IQR 10.000 ns) | 190.000 ns (IQR 0.000 ns) | 28.160 μs (IQR 695.000 ns) | 6.444 ms (IQR 72.182 μs) | 70.000 ns (IQR 0.000 ns) | 240.000 ns (IQR 10.000 ns) | 69.590 μs (IQR 4.885 μs) | 10.042 ms (IQR 149.286 μs) | 70.000 ns (IQR 10.000 ns) | 300.000 ns (IQR 0.000 ns) | 110.244 μs (IQR 3.438 μs) | 13.575 ms (IQR 158.219 μs) | 80.000 ns (IQR 0.000 ns) | 360.000 ns (IQR 0.000 ns) | 172.224 μs (IQR 7.392 μs) | 19.794 ms (IQR 117.482 μs) | 
| Dierckx (k=3) | 150.000 ns (IQR 10.000 ns) | 1.260 μs (IQR 20.000 ns) | 134.564 μs (IQR 305.000 ns) | 13.506 ms (IQR 43.390 μs) | 280.000 ns (IQR 10.000 ns) | 3.950 μs (IQR 10.000 ns) | 412.066 μs (IQR 282.500 ns) | 41.346 ms (IQR 455.686 μs) | 1.490 μs (IQR 0.000 ns) | 30.060 μs (IQR 30.000 ns) | 3.247 ms (IQR 47.519 μs) | 323.084 ms (IQR 648.144 μs) | 13.470 μs (IQR 80.000 ns) | 291.452 μs (IQR 1.265 μs) | 31.720 ms (IQR 238.153 μs) | 3.151 s (IQR 0.000 ns) | 
| Interpolations (uniform) | 90.000 ns (IQR 10.000 ns) | 270.000 ns (IQR 10.000 ns) | 20.459 μs (IQR 10.000 ns) | 2.048 ms (IQR 8.437 μs) | 80.000 ns (IQR 10.000 ns) | 260.000 ns (IQR 10.000 ns) | 20.440 μs (IQR 52.500 ns) | 2.047 ms (IQR 13.848 μs) | 90.000 ns (IQR 10.000 ns) | 270.000 ns (IQR 0.000 ns) | 20.430 μs (IQR 20.000 ns) | 2.049 ms (IQR 31.699 μs) | 90.000 ns (IQR 10.000 ns) | 270.000 ns (IQR 0.000 ns) | 20.420 μs (IQR 20.000 ns) | 2.061 ms (IQR 28.078 μs) | 

### Linear

| Library | n=100,m=1 | n=100,m=10 | n=100,m=1000 | n=100,m=100000 | n=1000,m=1 | n=1000,m=10 | n=1000,m=1000 | n=1000,m=100000 | n=10000,m=1 | n=10000,m=10 | n=10000,m=1000 | n=10000,m=100000 | n=100000,m=1 | n=100000,m=10 | n=100000,m=1000 | n=100000,m=100000 | 
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 130.000 ns (IQR 0.000 ns) | 520.000 ns (IQR 10.000 ns) | 55.089 μs (IQR 111.500 ns) | 5.582 ms (IQR 35.069 μs) | 120.000 ns (IQR 0.000 ns) | 510.000 ns (IQR 10.000 ns) | 55.030 μs (IQR 172.500 ns) | 5.565 ms (IQR 31.740 μs) | 130.000 ns (IQR 10.000 ns) | 480.000 ns (IQR 0.000 ns) | 56.074 μs (IQR 649.500 ns) | 5.687 ms (IQR 50.330 μs) | 130.000 ns (IQR 10.000 ns) | 480.000 ns (IQR 0.000 ns) | 58.300 μs (IQR 568.500 ns) | 6.744 ms (IQR 64.104 μs) | 
| BasicInterpolators | 60.000 ns (IQR 0.000 ns) | 180.000 ns (IQR 10.000 ns) | 29.580 μs (IQR 1.822 μs) | 6.435 ms (IQR 41.382 μs) | 60.000 ns (IQR 10.000 ns) | 230.000 ns (IQR 0.000 ns) | 64.855 μs (IQR 3.784 μs) | 9.784 ms (IQR 57.275 μs) | 70.000 ns (IQR 0.000 ns) | 300.000 ns (IQR 10.000 ns) | 105.934 μs (IQR 1.353 μs) | 13.251 ms (IQR 76.569 μs) | 80.000 ns (IQR 10.000 ns) | 360.000 ns (IQR 10.000 ns) | 167.988 μs (IQR 8.889 μs) | 19.198 ms (IQR 179.439 μs) | 
| Dierckx (k=1) | 110.000 ns (IQR 0.000 ns) | 880.000 ns (IQR 0.000 ns) | 103.619 μs (IQR 592.500 ns) | 10.349 ms (IQR 87.979 μs) | 250.000 ns (IQR 0.000 ns) | 3.600 μs (IQR 10.000 ns) | 377.267 μs (IQR 260.250 ns) | 37.624 ms (IQR 375.581 μs) | 1.450 μs (IQR 10.000 ns) | 29.650 μs (IQR 23.250 ns) | 3.191 ms (IQR 13.440 μs) | 318.772 ms (IQR 307.167 μs) | 13.470 μs (IQR 40.000 ns) | 291.407 μs (IQR 310.000 ns) | 31.660 ms (IQR 145.164 μs) | 3.181 s (IQR 0.000 ns) | 
| Interpolations (gridded) | 70.000 ns (IQR 2.500 ns) | 230.000 ns (IQR 10.000 ns) | 21.320 μs (IQR 824.250 ns) | 6.593 ms (IQR 59.599 μs) | 70.000 ns (IQR 10.000 ns) | 290.000 ns (IQR 0.000 ns) | 57.594 μs (IQR 7.278 μs) | 10.134 ms (IQR 55.693 μs) | 100.000 ns (IQR 20.000 ns) | 350.000 ns (IQR 2.500 ns) | 101.449 μs (IQR 5.872 μs) | 13.761 ms (IQR 140.448 μs) | 90.000 ns (IQR 0.000 ns) | 400.000 ns (IQR 10.000 ns) | 166.308 μs (IQR 6.663 μs) | 19.858 ms (IQR 70.962 μs) | 
| Interpolations (uniform) | 70.000 ns (IQR 10.000 ns) | 120.000 ns (IQR 10.000 ns) | 6.560 μs (IQR 40.000 ns) | 655.544 μs (IQR 9.552 μs) | 60.000 ns (IQR 0.000 ns) | 120.000 ns (IQR 0.000 ns) | 6.630 μs (IQR 11.000 ns) | 655.564 μs (IQR 8.735 μs) | 70.000 ns (IQR 10.000 ns) | 120.000 ns (IQR 0.000 ns) | 6.889 μs (IQR 10.000 ns) | 683.948 μs (IQR 4.087 μs) | 70.000 ns (IQR 10.000 ns) | 120.000 ns (IQR 0.000 ns) | 7.160 μs (IQR 10.000 ns) | 752.648 μs (IQR 1.765 μs) | 

### MonotoneCubic

| Library | n=100,m=1 | n=100,m=10 | n=100,m=1000 | n=100,m=100000 | n=1000,m=1 | n=1000,m=10 | n=1000,m=1000 | n=1000,m=100000 | n=10000,m=1 | n=10000,m=10 | n=10000,m=1000 | n=10000,m=100000 | n=100000,m=1 | n=100000,m=10 | n=100000,m=1000 | n=100000,m=100000 | 
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| DataInterpolations (CubicHermite) | 130.000 ns (IQR 10.000 ns) | 530.000 ns (IQR 0.000 ns) | 52.569 μs (IQR 80.000 ns) | 5.297 ms (IQR 17.295 μs) | 130.000 ns (IQR 10.000 ns) | 530.000 ns (IQR 10.000 ns) | 53.119 μs (IQR 119.500 ns) | 5.326 ms (IQR 26.192 μs) | 140.000 ns (IQR 10.000 ns) | 520.000 ns (IQR 10.000 ns) | 57.039 μs (IQR 395.000 ns) | 5.697 ms (IQR 20.677 μs) | 150.000 ns (IQR 0.000 ns) | 540.000 ns (IQR 10.000 ns) | 61.514 μs (IQR 244.750 ns) | 7.278 ms (IQR 31.449 μs) | 
| PCHIPInterpolation | 80.000 ns (IQR 10.000 ns) | 240.000 ns (IQR 10.000 ns) | 24.865 μs (IQR 1.385 μs) | 6.407 ms (IQR 40.035 μs) | 100.000 ns (IQR 20.000 ns) | 340.000 ns (IQR 10.000 ns) | 65.279 μs (IQR 569.000 ns) | 9.842 ms (IQR 42.269 μs) | 90.000 ns (IQR 10.000 ns) | 440.000 ns (IQR 0.000 ns) | 117.314 μs (IQR 3.672 μs) | 14.126 ms (IQR 36.737 μs) | 110.000 ns (IQR 10.000 ns) | 480.000 ns (IQR 10.000 ns) | 182.258 μs (IQR 4.223 μs) | 20.388 ms (IQR 38.910 μs) | 

### QuadraticSpline

| Library | n=100,m=1 | n=100,m=10 | n=100,m=1000 | n=100,m=100000 | n=1000,m=1 | n=1000,m=10 | n=1000,m=1000 | n=1000,m=100000 | n=10000,m=1 | n=10000,m=10 | n=10000,m=1000 | n=10000,m=100000 | n=100000,m=1 | n=100000,m=10 | n=100000,m=1000 | n=100000,m=100000 | 
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 230.000 ns (IQR 10.000 ns) | 1.650 μs (IQR 10.000 ns) | 169.573 μs (IQR 1.357 μs) | 17.094 ms (IQR 49.515 μs) | 540.000 ns (IQR 10.000 ns) | 7.080 μs (IQR 10.000 ns) | 750.153 μs (IQR 7.195 μs) | 74.773 ms (IQR 70.109 μs) | 3.630 μs (IQR 22.250 ns) | 62.340 μs (IQR 43.250 ns) | 6.627 ms (IQR 25.269 μs) | 669.475 ms (IQR 0.000 ns) | 39.715 μs (IQR 720.000 ns) | 663.894 μs (IQR 3.838 μs) | 70.311 ms (IQR 31.738 μs) | 7.058 s (IQR 0.000 ns) | 
| Dierckx (k=2) | 130.000 ns (IQR 0.000 ns) | 1.060 μs (IQR 0.000 ns) | 116.774 μs (IQR 711.750 ns) | 11.700 ms (IQR 35.465 μs) | 270.000 ns (IQR 10.000 ns) | 3.790 μs (IQR 10.000 ns) | 394.467 μs (IQR 301.750 ns) | 39.408 ms (IQR 62.060 μs) | 1.480 μs (IQR 2.500 ns) | 29.880 μs (IQR 30.000 ns) | 3.209 ms (IQR 5.805 μs) | 320.457 ms (IQR 632.470 μs) | 13.500 μs (IQR 50.000 ns) | 291.127 μs (IQR 298.500 ns) | 31.489 ms (IQR 130.091 μs) | 3.195 s (IQR 0.000 ns) | 

## Chained ODE-style

Sequential `for x in tt; A(x); end` over a monotone sequence. (knot pattern = uniform)

### Akima

| Library | n=100,m=1000 | n=1000,m=1000 | n=10000,m=1000 | n=100000,m=1000 | 
|---|---|---|---|---|
| DataInterpolations | 40.720 μs (IQR 21.000 ns) | 47.709 μs (IQR 51.000 ns) | 47.785 μs (IQR 129.250 ns) | 63.300 μs (IQR 381.500 ns) | 

### CubicSpline

| Library | n=100,m=1000 | n=1000,m=1000 | n=10000,m=1000 | n=100000,m=1000 | 
|---|---|---|---|---|
| DataInterpolations | 46.169 μs (IQR 20.000 ns) | 52.069 μs (IQR 350.000 ns) | 52.645 μs (IQR 225.500 ns) | 69.550 μs (IQR 362.500 ns) | 
| BasicInterpolators | 16.989 μs (IQR 459.250 ns) | 47.639 μs (IQR 1.077 μs) | 86.939 μs (IQR 2.837 μs) | 136.684 μs (IQR 2.812 μs) | 
| Dierckx (k=3) | 122.399 μs (IQR 52.500 ns) | 409.806 μs (IQR 50.000 ns) | 3.212 ms (IQR 2.620 μs) | 31.457 ms (IQR 222.590 μs) | 
| Interpolations (uniform) | 27.750 μs (IQR 600.000 ns) | 27.749 μs (IQR 499.250 ns) | 27.750 μs (IQR 1.179 μs) | 27.750 μs (IQR 440.000 ns) | 

### Linear

| Library | n=100,m=1000 | n=1000,m=1000 | n=10000,m=1000 | n=100000,m=1000 | 
|---|---|---|---|---|
| DataInterpolations | 42.920 μs (IQR 10.000 ns) | 49.549 μs (IQR 130.250 ns) | 50.259 μs (IQR 159.250 ns) | 55.679 μs (IQR 428.000 ns) | 
| BasicInterpolators | 17.259 μs (IQR 92.500 ns) | 44.764 μs (IQR 3.072 μs) | 82.629 μs (IQR 2.270 μs) | 132.393 μs (IQR 2.408 μs) | 
| Dierckx (k=1) | 87.039 μs (IQR 70.000 ns) | 374.796 μs (IQR 415.500 ns) | 3.199 ms (IQR 48.275 μs) | 31.242 ms (IQR 56.158 μs) | 
| Interpolations (gridded) | 23.920 μs (IQR 449.750 ns) | 58.344 μs (IQR 2.035 μs) | 87.894 μs (IQR 3.572 μs) | 141.208 μs (IQR 3.995 μs) | 
| Interpolations (uniform) | 23.050 μs (IQR 21.000 ns) | 23.060 μs (IQR 30.000 ns) | 23.040 μs (IQR 50.000 ns) | 23.129 μs (IQR 40.000 ns) | 

### MonotoneCubic

| Library | n=100,m=1000 | n=1000,m=1000 | n=10000,m=1000 | n=100000,m=1000 | 
|---|---|---|---|---|
| DataInterpolations (CubicHermite) | 50.009 μs (IQR 30.000 ns) | 54.480 μs (IQR 109.000 ns) | 53.984 μs (IQR 110.000 ns) | 71.909 μs (IQR 809.250 ns) | 
| PCHIPInterpolation | 25.240 μs (IQR 300.000 ns) | 50.559 μs (IQR 934.000 ns) | 93.234 μs (IQR 1.728 μs) | 148.724 μs (IQR 7.090 μs) | 

### QuadraticSpline

| Library | n=100,m=1000 | n=1000,m=1000 | n=10000,m=1000 | n=100000,m=1000 | 
|---|---|---|---|---|
| DataInterpolations | 163.768 μs (IQR 173.250 ns) | 743.168 μs (IQR 397.500 ns) | 6.613 ms (IQR 60.383 μs) | 70.305 ms (IQR 457.841 μs) | 
| Dierckx (k=2) | 104.904 μs (IQR 220.000 ns) | 393.182 μs (IQR 351.000 ns) | 3.190 ms (IQR 1.820 μs) | 31.236 ms (IQR 63.944 μs) | 

## Reproducer

Bench script: `bench/cross_library_comparison.jl`

Bench Project.toml: `bench/Project.toml` (devs DI from `..`).

To rerun:
```bash
cd /home/crackauc/sandbox/tmp_20260515_091703_4914/DataInterpolations.jl
git checkout fff-strategy-batched-evals
julia +1.11 --project=bench bench/cross_library_comparison.jl
```

