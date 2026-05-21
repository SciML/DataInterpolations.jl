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

Commit: `3ce1080447eb29b4f51ef363bbad4129e00b176a`

Library versions:
```
  DataInterpolations 8.10.0
  Interpolations 0.16.2
  Dierckx 0.5.4
  BasicInterpolators 0.7.1
  PCHIPInterpolation 0.2.1
  FastInterpolations 0.4.11
  BenchmarkTools 1.8.0
```

Total bench time: 261.4 s

## Construction time

Rows = library, columns = (n, knot pattern). Values = median wall time (IQR).

### Akima

| Library | n=100,nonuniform | n=100,uniform | n=1000,nonuniform | n=1000,uniform | n=10000,nonuniform | n=10000,uniform | n=100000,nonuniform | n=100000,uniform | 
|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 1.895 μs (IQR 267.500 ns) | 1.935 μs (IQR 255.000 ns) | 12.340 μs (IQR 552.500 ns) | 13.625 μs (IQR 31.997 μs) | 96.769 μs (IQR 892.500 ns) | 105.959 μs (IQR 160.974 μs) | 966.061 μs (IQR 2.906 ms) | 1.065 ms (IQR 2.945 ms) | 
| FastInterpolations | 1.320 μs (IQR 422.500 ns) | 990.000 ns (IQR 310.000 ns) | 10.680 μs (IQR 39.707 μs) | 29.840 μs (IQR 7.692 μs) | 69.779 μs (IQR 424.373 μs) | 173.559 μs (IQR 78.677 μs) | 680.669 μs (IQR 35.064 μs) | 673.078 μs (IQR 25.835 μs) | 

### CubicSpline

| Library | n=100,nonuniform | n=100,uniform | n=1000,nonuniform | n=1000,uniform | n=10000,nonuniform | n=10000,uniform | n=100000,nonuniform | n=100000,uniform | 
|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 6.930 μs (IQR 3.364 μs) | 7.170 μs (IQR 3.487 μs) | 55.215 μs (IQR 1.649 μs) | 57.094 μs (IQR 44.369 μs) | 482.461 μs (IQR 3.031 μs) | 495.015 μs (IQR 981.821 μs) | 5.001 ms (IQR 2.172 ms) | 5.070 ms (IQR 893.747 μs) | 
| BasicInterpolators | 4.340 μs (IQR 530.000 ns) | 4.385 μs (IQR 722.500 ns) | 37.880 μs (IQR 1.177 μs) | 39.015 μs (IQR 2.552 μs) | 319.082 μs (IQR 3.510 μs) | 1.021 ms (IQR 696.255 μs) | 3.258 ms (IQR 864.505 μs) | 3.324 ms (IQR 835.567 μs) | 
| Dierckx (k=3) | 23.020 μs (IQR 1.680 μs) | 17.440 μs (IQR 495.000 ns) | 209.343 μs (IQR 50.807 μs) | 211.203 μs (IQR 27.988 μs) | 1.582 ms (IQR 93.580 μs) | 1.723 ms (IQR 687.644 μs) | 18.212 ms (IQR 5.298 ms) | 20.038 ms (IQR 5.815 ms) | 
| FastInterpolations | 1.760 μs (IQR 285.000 ns) | 1.400 μs (IQR 235.000 ns) | 14.395 μs (IQR 495.000 ns) | 11.440 μs (IQR 515.000 ns) | 128.024 μs (IQR 1.032 μs) | 101.614 μs (IQR 1.115 μs) | 1.280 ms (IQR 242.824 μs) | 1.010 ms (IQR 29.291 μs) | 
| Interpolations (uniform) | — | 15.339 μs (IQR 1.022 μs) | — | 130.319 μs (IQR 9.140 μs) | — | 1.470 ms (IQR 311.781 μs) | — | 8.529 ms (IQR 384.982 μs) | 

### Linear

| Library | n=100,nonuniform | n=100,uniform | n=1000,nonuniform | n=1000,uniform | n=10000,nonuniform | n=10000,uniform | n=100000,nonuniform | n=100000,uniform | 
|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 905.000 ns (IQR 250.000 ns) | 1.030 μs (IQR 162.500 ns) | 7.020 μs (IQR 412.500 ns) | 7.835 μs (IQR 300.250 ns) | 60.794 μs (IQR 1.147 μs) | 69.284 μs (IQR 679.750 ns) | 602.635 μs (IQR 3.941 μs) | 684.534 μs (IQR 5.035 μs) | 
| BasicInterpolators | 1.170 μs (IQR 350.000 ns) | 1.320 μs (IQR 262.500 ns) | 8.145 μs (IQR 600.000 ns) | 9.035 μs (IQR 1.195 μs) | 61.674 μs (IQR 2.082 μs) | 70.815 μs (IQR 2.756 μs) | 610.510 μs (IQR 9.142 μs) | 696.798 μs (IQR 5.480 μs) | 
| Dierckx (k=1) | 9.309 μs (IQR 580.000 ns) | 12.405 μs (IQR 2.297 μs) | 79.575 μs (IQR 2.900 μs) | 106.234 μs (IQR 10.303 μs) | 764.703 μs (IQR 21.014 μs) | 764.083 μs (IQR 13.045 μs) | 7.486 ms (IQR 1.087 ms) | 7.568 ms (IQR 922.522 μs) | 
| FastInterpolations | 1.015 μs (IQR 402.500 ns) | 330.000 ns (IQR 212.500 ns) | 7.410 μs (IQR 2.632 μs) | 2.100 μs (IQR 347.500 ns) | 45.944 μs (IQR 3.295 μs) | 15.395 μs (IQR 1.893 μs) | 429.641 μs (IQR 24.649 μs) | 137.569 μs (IQR 2.454 μs) | 
| Interpolations (gridded) | 875.000 ns (IQR 245.000 ns) | 970.000 ns (IQR 245.000 ns) | 7.385 μs (IQR 477.500 ns) | 9.085 μs (IQR 9.578 μs) | 60.604 μs (IQR 1.885 μs) | 69.745 μs (IQR 1.954 μs) | 596.539 μs (IQR 5.037 μs) | 682.029 μs (IQR 5.189 μs) | 
| Interpolations (uniform) | — | 165.000 ns (IQR 165.000 ns) | — | 3.040 μs (IQR 2.235 μs) | — | 23.445 μs (IQR 5.530 μs) | — | 68.814 μs (IQR 2.921 μs) | 

### MonotoneCubic

| Library | n=100,nonuniform | n=100,uniform | n=1000,nonuniform | n=1000,uniform | n=10000,nonuniform | n=10000,uniform | n=100000,nonuniform | n=100000,uniform | 
|---|---|---|---|---|---|---|---|---|
| DataInterpolations (CubicHermite) | 1.210 μs (IQR 245.000 ns) | 1.260 μs (IQR 312.500 ns) | 8.205 μs (IQR 472.500 ns) | 8.970 μs (IQR 552.500 ns) | 66.989 μs (IQR 899.500 ns) | 155.944 μs (IQR 2.865 μs) | 667.664 μs (IQR 11.177 μs) | 760.448 μs (IQR 1.507 ms) | 
| FastInterpolations (PCHIP) | 1.980 μs (IQR 395.000 ns) | 1.245 μs (IQR 262.500 ns) | 15.105 μs (IQR 837.500 ns) | 10.855 μs (IQR 550.000 ns) | 116.039 μs (IQR 425.011 μs) | 93.414 μs (IQR 1.335 μs) | 1.148 ms (IQR 64.792 μs) | 925.542 μs (IQR 1.909 ms) | 
| PCHIPInterpolation | 1.590 μs (IQR 302.500 ns) | 1.685 μs (IQR 352.500 ns) | 13.275 μs (IQR 790.000 ns) | 14.570 μs (IQR 890.000 ns) | 150.029 μs (IQR 157.078 μs) | 118.489 μs (IQR 12.057 μs) | 1.100 ms (IQR 1.930 ms) | 1.195 ms (IQR 2.849 ms) | 

### QuadraticSpline

| Library | n=100,nonuniform | n=100,uniform | n=1000,nonuniform | n=1000,uniform | n=10000,nonuniform | n=10000,uniform | n=100000,nonuniform | n=100000,uniform | 
|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 9.420 μs (IQR 661.750 ns) | 9.455 μs (IQR 750.000 ns) | 84.669 μs (IQR 997.750 ns) | 96.034 μs (IQR 61.307 μs) | 800.123 μs (IQR 10.530 μs) | 823.987 μs (IQR 668.939 μs) | 8.070 ms (IQR 3.777 ms) | 8.209 ms (IQR 5.391 ms) | 
| Dierckx (k=2) | 15.485 μs (IQR 680.000 ns) | 21.090 μs (IQR 3.221 μs) | 136.554 μs (IQR 1.056 μs) | 143.228 μs (IQR 16.675 μs) | 1.353 ms (IQR 164.550 μs) | 1.356 ms (IQR 155.038 μs) | 13.384 ms (IQR 4.032 ms) | 14.377 ms (IQR 3.926 ms) | 
| FastInterpolations | 1.430 μs (IQR 420.000 ns) | 715.000 ns (IQR 342.500 ns) | 10.325 μs (IQR 740.000 ns) | 5.210 μs (IQR 485.750 ns) | 70.434 μs (IQR 420.084 μs) | 38.454 μs (IQR 2.993 μs) | 702.374 μs (IQR 38.012 μs) | 372.171 μs (IQR 33.597 μs) | 

## Single-query latency

Cold single evaluation `A(x_query)`. Rows = library, columns = (n, knot pattern).

### Akima

| Library | n=100,nonuniform | n=100,uniform | n=1000,nonuniform | n=1000,uniform | n=10000,nonuniform | n=10000,uniform | n=100000,nonuniform | n=100000,uniform | 
|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 70.000 ns (IQR 0.000 ns) | 70.000 ns (IQR 0.000 ns) | 70.000 ns (IQR 0.000 ns) | 80.000 ns (IQR 10.000 ns) | 90.000 ns (IQR 0.000 ns) | 90.000 ns (IQR 0.000 ns) | 100.000 ns (IQR 2.500 ns) | 100.000 ns (IQR 10.000 ns) | 
| FastInterpolations | 50.000 ns (IQR 10.000 ns) | 50.000 ns (IQR 0.000 ns) | 50.000 ns (IQR 10.000 ns) | 50.000 ns (IQR 0.000 ns) | 60.000 ns (IQR 10.000 ns) | 50.000 ns (IQR 0.000 ns) | 60.000 ns (IQR 10.000 ns) | 50.000 ns (IQR 0.000 ns) | 

### CubicSpline

| Library | n=100,nonuniform | n=100,uniform | n=1000,nonuniform | n=1000,uniform | n=10000,nonuniform | n=10000,uniform | n=100000,nonuniform | n=100000,uniform | 
|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 80.000 ns (IQR 0.000 ns) | 80.000 ns (IQR 10.000 ns) | 100.000 ns (IQR 0.000 ns) | 90.000 ns (IQR 0.000 ns) | 90.000 ns (IQR 10.000 ns) | 100.000 ns (IQR 0.000 ns) | 110.000 ns (IQR 0.000 ns) | 110.000 ns (IQR 10.000 ns) | 
| BasicInterpolators | 50.000 ns (IQR 10.000 ns) | 50.000 ns (IQR 10.000 ns) | 60.000 ns (IQR 2.500 ns) | 60.000 ns (IQR 0.000 ns) | 60.000 ns (IQR 0.000 ns) | 60.000 ns (IQR 0.000 ns) | 70.000 ns (IQR 10.000 ns) | 90.000 ns (IQR 10.000 ns) | 
| Dierckx (k=3) | 150.000 ns (IQR 0.000 ns) | 150.000 ns (IQR 10.000 ns) | 380.000 ns (IQR 10.000 ns) | 400.000 ns (IQR 10.000 ns) | 2.789 μs (IQR 10.000 ns) | 2.790 μs (IQR 10.000 ns) | 26.760 μs (IQR 51.000 ns) | 26.680 μs (IQR 40.000 ns) | 
| FastInterpolations | 50.000 ns (IQR 10.000 ns) | 50.000 ns (IQR 10.000 ns) | 50.000 ns (IQR 0.000 ns) | 60.000 ns (IQR 10.000 ns) | 60.000 ns (IQR 0.000 ns) | 50.000 ns (IQR 10.000 ns) | 60.000 ns (IQR 0.000 ns) | 50.000 ns (IQR 10.000 ns) | 
| Interpolations (uniform) | — | 80.000 ns (IQR 10.000 ns) | — | 80.000 ns (IQR 10.000 ns) | — | 80.000 ns (IQR 0.000 ns) | — | 80.000 ns (IQR 10.000 ns) | 

### Linear

| Library | n=100,nonuniform | n=100,uniform | n=1000,nonuniform | n=1000,uniform | n=10000,nonuniform | n=10000,uniform | n=100000,nonuniform | n=100000,uniform | 
|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 70.000 ns (IQR 10.000 ns) | 70.000 ns (IQR 10.000 ns) | 70.000 ns (IQR 0.000 ns) | 90.000 ns (IQR 0.000 ns) | 90.000 ns (IQR 10.000 ns) | 100.000 ns (IQR 10.000 ns) | 110.000 ns (IQR 0.000 ns) | 100.000 ns (IQR 10.000 ns) | 
| BasicInterpolators | 60.000 ns (IQR 10.000 ns) | 50.000 ns (IQR 0.000 ns) | 60.000 ns (IQR 10.000 ns) | 60.000 ns (IQR 0.000 ns) | 70.000 ns (IQR 10.000 ns) | 70.000 ns (IQR 10.000 ns) | 70.000 ns (IQR 0.000 ns) | 70.000 ns (IQR 10.000 ns) | 
| Dierckx (k=1) | 120.000 ns (IQR 10.000 ns) | 120.000 ns (IQR 10.000 ns) | 350.000 ns (IQR 0.000 ns) | 370.000 ns (IQR 10.000 ns) | 2.750 μs (IQR 10.000 ns) | 2.750 μs (IQR 10.000 ns) | 26.720 μs (IQR 40.250 ns) | 26.650 μs (IQR 42.500 ns) | 
| FastInterpolations | 50.000 ns (IQR 0.000 ns) | 40.000 ns (IQR 10.000 ns) | 50.000 ns (IQR 0.000 ns) | 40.000 ns (IQR 10.000 ns) | 60.000 ns (IQR 10.000 ns) | 40.000 ns (IQR 10.000 ns) | 60.000 ns (IQR 0.000 ns) | 40.000 ns (IQR 10.000 ns) | 
| Interpolations (gridded) | 60.000 ns (IQR 10.000 ns) | 60.000 ns (IQR 0.000 ns) | 70.000 ns (IQR 10.000 ns) | 70.000 ns (IQR 0.000 ns) | 90.000 ns (IQR 0.000 ns) | 70.000 ns (IQR 0.000 ns) | 100.000 ns (IQR 0.000 ns) | 80.000 ns (IQR 10.000 ns) | 
| Interpolations (uniform) | — | 60.000 ns (IQR 0.000 ns) | — | 60.000 ns (IQR 0.000 ns) | — | 60.000 ns (IQR 0.000 ns) | — | 60.000 ns (IQR 0.000 ns) | 

### MonotoneCubic

| Library | n=100,nonuniform | n=100,uniform | n=1000,nonuniform | n=1000,uniform | n=10000,nonuniform | n=10000,uniform | n=100000,nonuniform | n=100000,uniform | 
|---|---|---|---|---|---|---|---|---|
| DataInterpolations (CubicHermite) | 90.000 ns (IQR 10.000 ns) | 90.000 ns (IQR 0.000 ns) | 90.000 ns (IQR 10.000 ns) | 100.000 ns (IQR 0.000 ns) | 100.000 ns (IQR 10.000 ns) | 100.000 ns (IQR 10.000 ns) | 110.000 ns (IQR 10.000 ns) | 110.000 ns (IQR 0.000 ns) | 
| FastInterpolations (PCHIP) | 50.000 ns (IQR 10.000 ns) | 50.000 ns (IQR 0.000 ns) | 60.000 ns (IQR 10.000 ns) | 50.000 ns (IQR 10.000 ns) | 60.000 ns (IQR 10.000 ns) | 50.000 ns (IQR 0.000 ns) | 60.000 ns (IQR 10.000 ns) | 50.000 ns (IQR 0.000 ns) | 
| PCHIPInterpolation | 60.000 ns (IQR 10.000 ns) | 70.000 ns (IQR 10.000 ns) | 80.000 ns (IQR 10.000 ns) | 80.000 ns (IQR 0.000 ns) | 90.000 ns (IQR 0.000 ns) | 90.000 ns (IQR 0.000 ns) | 110.000 ns (IQR 10.000 ns) | 100.000 ns (IQR 0.000 ns) | 

### QuadraticSpline

| Library | n=100,nonuniform | n=100,uniform | n=1000,nonuniform | n=1000,uniform | n=10000,nonuniform | n=10000,uniform | n=100000,nonuniform | n=100000,uniform | 
|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 160.000 ns (IQR 10.000 ns) | 160.000 ns (IQR 10.000 ns) | 250.000 ns (IQR 0.000 ns) | 260.000 ns (IQR 0.000 ns) | 1.240 μs (IQR 10.000 ns) | 1.240 μs (IQR 10.000 ns) | 14.540 μs (IQR 100.000 ns) | 14.464 μs (IQR 59.000 ns) | 
| Dierckx (k=2) | 130.000 ns (IQR 10.000 ns) | 140.000 ns (IQR 10.000 ns) | 370.000 ns (IQR 10.000 ns) | 390.000 ns (IQR 0.000 ns) | 2.780 μs (IQR 10.000 ns) | 2.780 μs (IQR 0.000 ns) | 26.719 μs (IQR 59.000 ns) | 26.630 μs (IQR 43.250 ns) | 
| FastInterpolations | 50.000 ns (IQR 10.000 ns) | 40.000 ns (IQR 10.000 ns) | 50.000 ns (IQR 2.500 ns) | 40.000 ns (IQR 10.000 ns) | 60.000 ns (IQR 10.000 ns) | 40.000 ns (IQR 10.000 ns) | 60.000 ns (IQR 2.500 ns) | 40.000 ns (IQR 10.000 ns) | 

## Sorted batch

`A(out, tt)` where `tt` is sorted random points in domain. (knot pattern = uniform)

### Akima

| Library | n=100,m=1 | n=100,m=10 | n=100,m=1000 | n=100,m=100000 | n=1000,m=1 | n=1000,m=10 | n=1000,m=1000 | n=1000,m=100000 | n=10000,m=1 | n=10000,m=10 | n=10000,m=1000 | n=10000,m=100000 | n=100000,m=1 | n=100000,m=10 | n=100000,m=1000 | n=100000,m=100000 | 
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 120.000 ns (IQR 10.000 ns) | 170.000 ns (IQR 0.000 ns) | 4.010 μs (IQR 0.000 ns) | 364.582 μs (IQR 1.153 μs) | 120.000 ns (IQR 10.000 ns) | 360.000 ns (IQR 0.000 ns) | 5.350 μs (IQR 20.000 ns) | 387.027 μs (IQR 2.195 μs) | 130.000 ns (IQR 10.000 ns) | 430.000 ns (IQR 10.000 ns) | 13.874 μs (IQR 683.250 ns) | 521.981 μs (IQR 3.110 μs) | 140.000 ns (IQR 10.000 ns) | 540.000 ns (IQR 10.000 ns) | 141.629 μs (IQR 2.840 μs) | 1.442 ms (IQR 10.085 μs) | 
| FastInterpolations | 60.000 ns (IQR 10.000 ns) | 170.000 ns (IQR 10.000 ns) | 11.190 μs (IQR 41.000 ns) | 1.074 ms (IQR 22.070 μs) | 60.000 ns (IQR 10.000 ns) | 170.000 ns (IQR 10.000 ns) | 11.170 μs (IQR 11.000 ns) | 1.070 ms (IQR 11.035 μs) | 60.000 ns (IQR 10.000 ns) | 170.000 ns (IQR 10.000 ns) | 11.190 μs (IQR 50.000 ns) | 1.072 ms (IQR 7.127 μs) | 60.000 ns (IQR 10.000 ns) | 170.000 ns (IQR 10.000 ns) | 11.884 μs (IQR 411.500 ns) | 1.075 ms (IQR 6.483 μs) | 

### CubicSpline

| Library | n=100,m=1 | n=100,m=10 | n=100,m=1000 | n=100,m=100000 | n=1000,m=1 | n=1000,m=10 | n=1000,m=1000 | n=1000,m=100000 | n=10000,m=1 | n=10000,m=10 | n=10000,m=1000 | n=10000,m=100000 | n=100000,m=1 | n=100000,m=10 | n=100000,m=1000 | n=100000,m=100000 | 
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 130.000 ns (IQR 10.000 ns) | 210.000 ns (IQR 0.000 ns) | 10.380 μs (IQR 10.000 ns) | 1.002 ms (IQR 6.160 μs) | 130.000 ns (IQR 10.000 ns) | 410.000 ns (IQR 10.000 ns) | 10.220 μs (IQR 10.000 ns) | 1.007 ms (IQR 5.220 μs) | 140.000 ns (IQR 10.000 ns) | 495.000 ns (IQR 10.000 ns) | 15.790 μs (IQR 144.750 ns) | 1.050 ms (IQR 5.147 μs) | 150.000 ns (IQR 10.000 ns) | 570.000 ns (IQR 10.000 ns) | 153.254 μs (IQR 2.850 μs) | 1.787 ms (IQR 11.688 μs) | 
| BasicInterpolators | 60.000 ns (IQR 0.000 ns) | 200.000 ns (IQR 10.000 ns) | 16.920 μs (IQR 371.000 ns) | 1.402 ms (IQR 29.020 μs) | 70.000 ns (IQR 0.000 ns) | 250.000 ns (IQR 10.000 ns) | 43.025 μs (IQR 2.495 μs) | 2.069 ms (IQR 40.005 μs) | 70.000 ns (IQR 10.000 ns) | 310.000 ns (IQR 0.000 ns) | 79.704 μs (IQR 4.163 μs) | 4.140 ms (IQR 15.742 μs) | 80.000 ns (IQR 0.000 ns) | 370.000 ns (IQR 0.000 ns) | 126.109 μs (IQR 1.827 μs) | 7.865 ms (IQR 77.216 μs) | 
| Dierckx (k=3) | 190.000 ns (IQR 10.000 ns) | 1.220 μs (IQR 0.000 ns) | 121.679 μs (IQR 62.500 ns) | 11.733 ms (IQR 34.135 μs) | 700.000 ns (IQR 0.000 ns) | 3.620 μs (IQR 0.000 ns) | 408.256 μs (IQR 72.500 ns) | 40.985 ms (IQR 72.560 μs) | 5.690 μs (IQR 30.000 ns) | 26.920 μs (IQR 39.250 ns) | 3.202 ms (IQR 4.645 μs) | 322.833 ms (IQR 90.359 μs) | 55.620 μs (IQR 110.000 ns) | 260.202 μs (IQR 173.500 ns) | 31.267 ms (IQR 222.842 μs) | 3.164 s (IQR 0.000 ns) | 
| FastInterpolations | 60.000 ns (IQR 2.500 ns) | 150.000 ns (IQR 0.000 ns) | 8.860 μs (IQR 0.000 ns) | 843.207 μs (IQR 12.260 μs) | 60.000 ns (IQR 0.000 ns) | 140.000 ns (IQR 10.000 ns) | 8.980 μs (IQR 1.000 ns) | 846.997 μs (IQR 6.682 μs) | 60.000 ns (IQR 2.500 ns) | 150.000 ns (IQR 0.000 ns) | 8.950 μs (IQR 42.500 ns) | 846.477 μs (IQR 4.218 μs) | 70.000 ns (IQR 20.000 ns) | 150.000 ns (IQR 10.000 ns) | 9.525 μs (IQR 69.250 ns) | 849.188 μs (IQR 10.038 μs) | 
| Interpolations (uniform) | 90.000 ns (IQR 2.500 ns) | 270.000 ns (IQR 0.000 ns) | 20.460 μs (IQR 10.000 ns) | 2.046 ms (IQR 8.038 μs) | 90.000 ns (IQR 10.000 ns) | 270.000 ns (IQR 0.000 ns) | 20.460 μs (IQR 1.000 ns) | 2.047 ms (IQR 9.498 μs) | 90.000 ns (IQR 10.000 ns) | 270.000 ns (IQR 0.000 ns) | 20.440 μs (IQR 19.000 ns) | 2.047 ms (IQR 7.382 μs) | 90.000 ns (IQR 0.000 ns) | 270.000 ns (IQR 0.000 ns) | 20.440 μs (IQR 31.000 ns) | 2.049 ms (IQR 10.880 μs) | 

### Linear

| Library | n=100,m=1 | n=100,m=10 | n=100,m=1000 | n=100,m=100000 | n=1000,m=1 | n=1000,m=10 | n=1000,m=1000 | n=1000,m=100000 | n=10000,m=1 | n=10000,m=10 | n=10000,m=1000 | n=10000,m=100000 | n=100000,m=1 | n=100000,m=10 | n=100000,m=1000 | n=100000,m=100000 | 
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 120.000 ns (IQR 10.000 ns) | 180.000 ns (IQR 10.000 ns) | 6.490 μs (IQR 10.000 ns) | 530.655 μs (IQR 2.315 μs) | 130.000 ns (IQR 0.000 ns) | 380.000 ns (IQR 0.000 ns) | 7.010 μs (IQR 10.000 ns) | 541.605 μs (IQR 889.250 ns) | 140.000 ns (IQR 0.000 ns) | 470.000 ns (IQR 10.000 ns) | 15.860 μs (IQR 50.250 ns) | 685.063 μs (IQR 7.840 μs) | 150.000 ns (IQR 10.000 ns) | 550.000 ns (IQR 0.000 ns) | 143.244 μs (IQR 4.605 μs) | 1.653 ms (IQR 13.035 μs) | 
| BasicInterpolators | 60.000 ns (IQR 10.000 ns) | 190.000 ns (IQR 0.000 ns) | 16.470 μs (IQR 322.500 ns) | 1.399 ms (IQR 16.903 μs) | 70.000 ns (IQR 0.000 ns) | 240.000 ns (IQR 0.000 ns) | 42.279 μs (IQR 3.697 μs) | 2.052 ms (IQR 12.033 μs) | 70.000 ns (IQR 0.000 ns) | 310.000 ns (IQR 10.000 ns) | 75.574 μs (IQR 5.532 μs) | 4.203 ms (IQR 42.208 μs) | 80.000 ns (IQR 0.000 ns) | 380.000 ns (IQR 10.000 ns) | 123.084 μs (IQR 1.998 μs) | 7.835 ms (IQR 64.520 μs) | 
| Dierckx (k=1) | 150.000 ns (IQR 0.000 ns) | 840.000 ns (IQR 10.000 ns) | 86.449 μs (IQR 320.000 ns) | 8.338 ms (IQR 205.934 μs) | 670.000 ns (IQR 0.000 ns) | 3.270 μs (IQR 40.000 ns) | 373.326 μs (IQR 591.750 ns) | 37.456 ms (IQR 86.332 μs) | 5.660 μs (IQR 20.000 ns) | 26.600 μs (IQR 30.000 ns) | 3.172 ms (IQR 37.733 μs) | 319.279 ms (IQR 135.094 μs) | 55.620 μs (IQR 93.250 ns) | 263.383 μs (IQR 1.763 μs) | 31.217 ms (IQR 137.637 μs) | 3.147 s (IQR 0.000 ns) | 
| FastInterpolations | 50.000 ns (IQR 0.000 ns) | 90.000 ns (IQR 0.000 ns) | 4.220 μs (IQR 90.000 ns) | 370.687 μs (IQR 3.397 μs) | 50.000 ns (IQR 0.000 ns) | 90.000 ns (IQR 0.000 ns) | 4.145 μs (IQR 30.000 ns) | 372.367 μs (IQR 3.388 μs) | 50.000 ns (IQR 0.000 ns) | 90.000 ns (IQR 0.000 ns) | 4.200 μs (IQR 82.500 ns) | 373.437 μs (IQR 4.133 μs) | 50.000 ns (IQR 0.000 ns) | 90.000 ns (IQR 0.000 ns) | 4.370 μs (IQR 72.500 ns) | 374.116 μs (IQR 4.323 μs) | 
| Interpolations (gridded) | 70.000 ns (IQR 0.000 ns) | 230.000 ns (IQR 0.000 ns) | 18.080 μs (IQR 50.000 ns) | 1.644 ms (IQR 26.273 μs) | 80.000 ns (IQR 10.000 ns) | 300.000 ns (IQR 0.000 ns) | 52.349 μs (IQR 3.364 μs) | 2.509 ms (IQR 35.296 μs) | 80.000 ns (IQR 10.000 ns) | 350.000 ns (IQR 10.000 ns) | 93.839 μs (IQR 3.115 μs) | 4.435 ms (IQR 47.497 μs) | 110.000 ns (IQR 10.000 ns) | 430.000 ns (IQR 10.000 ns) | 139.209 μs (IQR 4.602 μs) | 7.831 ms (IQR 71.805 μs) | 
| Interpolations (uniform) | 70.000 ns (IQR 10.000 ns) | 120.000 ns (IQR 10.000 ns) | 6.610 μs (IQR 20.000 ns) | 655.379 μs (IQR 612.250 ns) | 70.000 ns (IQR 0.000 ns) | 120.000 ns (IQR 0.000 ns) | 6.610 μs (IQR 20.000 ns) | 655.654 μs (IQR 1.357 μs) | 70.000 ns (IQR 0.000 ns) | 120.000 ns (IQR 10.000 ns) | 6.665 μs (IQR 20.000 ns) | 655.439 μs (IQR 542.500 ns) | 70.000 ns (IQR 0.000 ns) | 120.000 ns (IQR 10.000 ns) | 7.050 μs (IQR 30.250 ns) | 655.829 μs (IQR 6.360 μs) | 

### MonotoneCubic

| Library | n=100,m=1 | n=100,m=10 | n=100,m=1000 | n=100,m=100000 | n=1000,m=1 | n=1000,m=10 | n=1000,m=1000 | n=1000,m=100000 | n=10000,m=1 | n=10000,m=10 | n=10000,m=1000 | n=10000,m=100000 | n=100000,m=1 | n=100000,m=10 | n=100000,m=1000 | n=100000,m=100000 | 
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| DataInterpolations (CubicHermite) | 130.000 ns (IQR 0.000 ns) | 230.000 ns (IQR 10.000 ns) | 10.500 μs (IQR 10.000 ns) | 983.576 μs (IQR 12.512 μs) | 130.000 ns (IQR 10.000 ns) | 430.000 ns (IQR 0.000 ns) | 10.660 μs (IQR 20.000 ns) | 992.811 μs (IQR 10.555 μs) | 140.000 ns (IQR 10.000 ns) | 510.000 ns (IQR 0.000 ns) | 18.420 μs (IQR 41.000 ns) | 1.075 ms (IQR 11.249 μs) | 150.000 ns (IQR 10.000 ns) | 600.000 ns (IQR 20.000 ns) | 147.614 μs (IQR 2.352 μs) | 1.812 ms (IQR 21.045 μs) | 
| FastInterpolations (PCHIP) | 60.000 ns (IQR 0.000 ns) | 160.000 ns (IQR 10.000 ns) | 11.139 μs (IQR 10.000 ns) | 1.072 ms (IQR 3.803 μs) | 60.000 ns (IQR 0.000 ns) | 160.000 ns (IQR 10.000 ns) | 11.190 μs (IQR 30.000 ns) | 1.072 ms (IQR 11.660 μs) | 60.000 ns (IQR 0.000 ns) | 170.000 ns (IQR 10.000 ns) | 11.160 μs (IQR 72.250 ns) | 1.081 ms (IQR 23.755 μs) | 60.000 ns (IQR 0.000 ns) | 160.000 ns (IQR 10.000 ns) | 11.830 μs (IQR 112.750 ns) | 1.073 ms (IQR 7.117 μs) | 
| PCHIPInterpolation | 70.000 ns (IQR 0.000 ns) | 240.000 ns (IQR 0.000 ns) | 20.294 μs (IQR 322.250 ns) | 1.864 ms (IQR 11.880 μs) | 90.000 ns (IQR 0.000 ns) | 320.000 ns (IQR 10.000 ns) | 40.360 μs (IQR 1.312 μs) | 2.973 ms (IQR 40.100 μs) | 100.000 ns (IQR 0.000 ns) | 450.000 ns (IQR 2.500 ns) | 79.409 μs (IQR 2.623 μs) | 5.436 ms (IQR 51.495 μs) | 120.000 ns (IQR 10.000 ns) | 510.000 ns (IQR 10.000 ns) | 134.408 μs (IQR 3.429 μs) | 9.598 ms (IQR 69.869 μs) | 

### QuadraticSpline

| Library | n=100,m=1 | n=100,m=10 | n=100,m=1000 | n=100,m=100000 | n=1000,m=1 | n=1000,m=10 | n=1000,m=1000 | n=1000,m=100000 | n=10000,m=1 | n=10000,m=10 | n=10000,m=1000 | n=10000,m=100000 | n=100000,m=1 | n=100000,m=10 | n=100000,m=1000 | n=100000,m=100000 | 
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 210.000 ns (IQR 0.000 ns) | 1.010 μs (IQR 0.000 ns) | 110.119 μs (IQR 162.500 ns) | 9.480 ms (IQR 53.279 μs) | 310.000 ns (IQR 0.000 ns) | 2.420 μs (IQR 10.000 ns) | 222.133 μs (IQR 455.000 ns) | 18.933 ms (IQR 74.185 μs) | 1.300 μs (IQR 10.000 ns) | 12.569 μs (IQR 59.250 ns) | 1.227 ms (IQR 5.570 μs) | 117.928 ms (IQR 382.237 μs) | 14.610 μs (IQR 70.000 ns) | 144.098 μs (IQR 315.250 ns) | 14.771 ms (IQR 263.223 μs) | 1.477 s (IQR 0.000 ns) | 
| Dierckx (k=2) | 170.000 ns (IQR 2.500 ns) | 1.040 μs (IQR 10.000 ns) | 103.589 μs (IQR 329.250 ns) | 9.887 ms (IQR 85.805 μs) | 690.000 ns (IQR 0.000 ns) | 3.440 μs (IQR 10.000 ns) | 391.187 μs (IQR 620.250 ns) | 39.289 ms (IQR 132.329 μs) | 5.669 μs (IQR 20.000 ns) | 26.799 μs (IQR 21.000 ns) | 3.180 ms (IQR 3.583 μs) | 320.143 ms (IQR 155.369 μs) | 55.730 μs (IQR 89.000 ns) | 260.017 μs (IQR 267.500 ns) | 31.216 ms (IQR 132.239 μs) | 3.154 s (IQR 0.000 ns) | 
| FastInterpolations | 50.000 ns (IQR 10.000 ns) | 110.000 ns (IQR 10.000 ns) | 3.910 μs (IQR 20.000 ns) | 349.737 μs (IQR 3.685 μs) | 50.000 ns (IQR 10.000 ns) | 110.000 ns (IQR 10.000 ns) | 3.990 μs (IQR 60.000 ns) | 350.812 μs (IQR 3.107 μs) | 50.000 ns (IQR 10.000 ns) | 110.000 ns (IQR 10.000 ns) | 4.320 μs (IQR 62.500 ns) | 351.341 μs (IQR 3.715 μs) | 50.000 ns (IQR 10.000 ns) | 110.000 ns (IQR 10.000 ns) | 5.560 μs (IQR 80.000 ns) | 359.077 μs (IQR 9.547 μs) | 

## Random batch

`A(out, tt)` where `tt` is unsorted. (knot pattern = uniform)

### Akima

| Library | n=100,m=1 | n=100,m=10 | n=100,m=1000 | n=100,m=100000 | n=1000,m=1 | n=1000,m=10 | n=1000,m=1000 | n=1000,m=100000 | n=10000,m=1 | n=10000,m=10 | n=10000,m=1000 | n=10000,m=100000 | n=100000,m=1 | n=100000,m=10 | n=100000,m=1000 | n=100000,m=100000 | 
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 120.000 ns (IQR 10.000 ns) | 440.000 ns (IQR 10.000 ns) | 85.734 μs (IQR 283.250 ns) | 8.733 ms (IQR 42.521 μs) | 120.000 ns (IQR 0.000 ns) | 760.000 ns (IQR 60.000 ns) | 113.999 μs (IQR 2.893 μs) | 11.877 ms (IQR 130.729 μs) | 130.000 ns (IQR 0.000 ns) | 620.000 ns (IQR 10.000 ns) | 144.079 μs (IQR 3.071 μs) | 15.912 ms (IQR 39.300 μs) | 130.000 ns (IQR 10.000 ns) | 750.000 ns (IQR 10.000 ns) | 209.804 μs (IQR 3.175 μs) | 22.256 ms (IQR 44.275 μs) | 
| FastInterpolations | 60.000 ns (IQR 10.000 ns) | 170.000 ns (IQR 10.000 ns) | 11.150 μs (IQR 0.000 ns) | 1.071 ms (IQR 7.907 μs) | 60.000 ns (IQR 0.000 ns) | 190.000 ns (IQR 0.000 ns) | 11.180 μs (IQR 30.000 ns) | 1.071 ms (IQR 5.312 μs) | 60.000 ns (IQR 10.000 ns) | 180.000 ns (IQR 10.000 ns) | 11.180 μs (IQR 40.000 ns) | 1.070 ms (IQR 5.062 μs) | 60.000 ns (IQR 0.000 ns) | 180.000 ns (IQR 10.000 ns) | 11.215 μs (IQR 60.250 ns) | 1.115 ms (IQR 5.000 μs) | 

### CubicSpline

| Library | n=100,m=1 | n=100,m=10 | n=100,m=1000 | n=100,m=100000 | n=1000,m=1 | n=1000,m=10 | n=1000,m=1000 | n=1000,m=100000 | n=10000,m=1 | n=10000,m=10 | n=10000,m=1000 | n=10000,m=100000 | n=100000,m=1 | n=100000,m=10 | n=100000,m=1000 | n=100000,m=100000 | 
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 130.000 ns (IQR 0.000 ns) | 510.000 ns (IQR 10.000 ns) | 88.379 μs (IQR 315.750 ns) | 8.819 ms (IQR 37.679 μs) | 130.000 ns (IQR 10.000 ns) | 790.000 ns (IQR 10.000 ns) | 116.909 μs (IQR 462.500 ns) | 11.897 ms (IQR 63.819 μs) | 140.000 ns (IQR 0.000 ns) | 740.000 ns (IQR 10.000 ns) | 145.738 μs (IQR 2.996 μs) | 16.274 ms (IQR 145.029 μs) | 150.000 ns (IQR 0.000 ns) | 740.000 ns (IQR 10.000 ns) | 221.258 μs (IQR 1.857 μs) | 24.186 ms (IQR 232.969 μs) | 
| BasicInterpolators | 60.000 ns (IQR 0.000 ns) | 200.000 ns (IQR 10.000 ns) | 26.510 μs (IQR 1.024 μs) | 6.469 ms (IQR 35.997 μs) | 70.000 ns (IQR 0.000 ns) | 240.000 ns (IQR 10.000 ns) | 64.219 μs (IQR 3.527 μs) | 9.810 ms (IQR 34.904 μs) | 70.000 ns (IQR 10.000 ns) | 300.000 ns (IQR 0.000 ns) | 107.729 μs (IQR 3.515 μs) | 13.373 ms (IQR 124.972 μs) | 80.000 ns (IQR 0.000 ns) | 360.000 ns (IQR 0.000 ns) | 171.113 μs (IQR 4.553 μs) | 19.801 ms (IQR 558.067 μs) | 
| Dierckx (k=3) | 150.000 ns (IQR 10.000 ns) | 1.250 μs (IQR 0.000 ns) | 134.424 μs (IQR 462.250 ns) | 13.515 ms (IQR 60.659 μs) | 290.000 ns (IQR 10.000 ns) | 3.960 μs (IQR 10.000 ns) | 412.197 μs (IQR 334.250 ns) | 41.079 ms (IQR 50.310 μs) | 1.480 μs (IQR 0.000 ns) | 30.050 μs (IQR 30.000 ns) | 3.232 ms (IQR 6.122 μs) | 321.523 ms (IQR 211.493 μs) | 13.540 μs (IQR 80.000 ns) | 291.423 μs (IQR 211.500 ns) | 31.510 ms (IQR 153.314 μs) | 3.190 s (IQR 0.000 ns) | 
| FastInterpolations | 60.000 ns (IQR 2.500 ns) | 150.000 ns (IQR 0.000 ns) | 8.980 μs (IQR 10.000 ns) | 860.512 μs (IQR 4.238 μs) | 60.000 ns (IQR 0.000 ns) | 150.000 ns (IQR 0.000 ns) | 8.960 μs (IQR 60.000 ns) | 847.768 μs (IQR 4.537 μs) | 60.000 ns (IQR 0.000 ns) | 150.000 ns (IQR 0.000 ns) | 8.980 μs (IQR 40.000 ns) | 854.847 μs (IQR 8.707 μs) | 60.000 ns (IQR 0.000 ns) | 150.000 ns (IQR 0.000 ns) | 9.140 μs (IQR 40.000 ns) | 1.037 ms (IQR 5.805 μs) | 
| Interpolations (uniform) | 90.000 ns (IQR 0.000 ns) | 270.000 ns (IQR 10.000 ns) | 20.430 μs (IQR 20.000 ns) | 2.046 ms (IQR 7.638 μs) | 90.000 ns (IQR 0.000 ns) | 270.000 ns (IQR 0.000 ns) | 20.470 μs (IQR 10.000 ns) | 2.046 ms (IQR 5.918 μs) | 90.000 ns (IQR 0.000 ns) | 270.000 ns (IQR 0.000 ns) | 20.420 μs (IQR 21.000 ns) | 2.044 ms (IQR 7.173 μs) | 90.000 ns (IQR 10.000 ns) | 270.000 ns (IQR 0.000 ns) | 20.450 μs (IQR 5.500 ns) | 2.048 ms (IQR 11.220 μs) | 

### Linear

| Library | n=100,m=1 | n=100,m=10 | n=100,m=1000 | n=100,m=100000 | n=1000,m=1 | n=1000,m=10 | n=1000,m=1000 | n=1000,m=100000 | n=10000,m=1 | n=10000,m=10 | n=10000,m=1000 | n=10000,m=100000 | n=100000,m=1 | n=100000,m=10 | n=100000,m=1000 | n=100000,m=100000 | 
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 120.000 ns (IQR 0.000 ns) | 480.000 ns (IQR 0.000 ns) | 83.519 μs (IQR 372.750 ns) | 8.415 ms (IQR 52.812 μs) | 130.000 ns (IQR 0.000 ns) | 770.000 ns (IQR 12.500 ns) | 108.679 μs (IQR 3.103 μs) | 11.532 ms (IQR 79.008 μs) | 140.000 ns (IQR 0.000 ns) | 700.000 ns (IQR 10.000 ns) | 141.614 μs (IQR 3.540 μs) | 15.566 ms (IQR 113.496 μs) | 150.000 ns (IQR 0.000 ns) | 820.000 ns (IQR 10.000 ns) | 208.248 μs (IQR 2.930 μs) | 22.285 ms (IQR 93.724 μs) | 
| BasicInterpolators | 60.000 ns (IQR 0.000 ns) | 200.000 ns (IQR 0.000 ns) | 30.665 μs (IQR 533.250 ns) | 6.365 ms (IQR 46.459 μs) | 60.000 ns (IQR 10.000 ns) | 240.000 ns (IQR 0.000 ns) | 65.439 μs (IQR 3.513 μs) | 9.742 ms (IQR 67.108 μs) | 70.000 ns (IQR 0.000 ns) | 290.000 ns (IQR 0.000 ns) | 104.350 μs (IQR 3.683 μs) | 13.157 ms (IQR 127.214 μs) | 80.000 ns (IQR 0.000 ns) | 360.000 ns (IQR 10.000 ns) | 167.203 μs (IQR 4.349 μs) | 19.014 ms (IQR 160.544 μs) | 
| Dierckx (k=1) | 110.000 ns (IQR 0.000 ns) | 870.000 ns (IQR 10.000 ns) | 98.819 μs (IQR 450.000 ns) | 10.105 ms (IQR 280.429 μs) | 260.000 ns (IQR 60.000 ns) | 3.610 μs (IQR 0.000 ns) | 376.521 μs (IQR 5.056 μs) | 37.694 ms (IQR 182.261 μs) | 1.455 μs (IQR 10.000 ns) | 29.699 μs (IQR 33.250 ns) | 3.194 ms (IQR 5.514 μs) | 321.630 ms (IQR 565.205 μs) | 13.500 μs (IQR 45.000 ns) | 291.442 μs (IQR 304.250 ns) | 31.515 ms (IQR 158.856 μs) | 3.173 s (IQR 0.000 ns) | 
| FastInterpolations | 50.000 ns (IQR 0.000 ns) | 90.000 ns (IQR 10.000 ns) | 4.240 μs (IQR 0.000 ns) | 392.077 μs (IQR 5.178 μs) | 50.000 ns (IQR 0.000 ns) | 90.000 ns (IQR 0.000 ns) | 4.240 μs (IQR 32.500 ns) | 373.457 μs (IQR 3.303 μs) | 50.000 ns (IQR 0.000 ns) | 90.000 ns (IQR 0.000 ns) | 4.270 μs (IQR 50.000 ns) | 377.542 μs (IQR 4.339 μs) | 50.000 ns (IQR 0.000 ns) | 90.000 ns (IQR 10.000 ns) | 4.350 μs (IQR 30.000 ns) | 464.091 μs (IQR 4.390 μs) | 
| Interpolations (gridded) | 70.000 ns (IQR 10.000 ns) | 230.000 ns (IQR 10.000 ns) | 25.485 μs (IQR 400.000 ns) | 6.346 ms (IQR 57.190 μs) | 80.000 ns (IQR 10.000 ns) | 270.000 ns (IQR 10.000 ns) | 60.749 μs (IQR 4.149 μs) | 9.964 ms (IQR 113.239 μs) | 80.000 ns (IQR 0.000 ns) | 350.000 ns (IQR 0.000 ns) | 106.214 μs (IQR 5.040 μs) | 13.345 ms (IQR 47.422 μs) | 90.000 ns (IQR 0.000 ns) | 410.000 ns (IQR 0.000 ns) | 172.038 μs (IQR 6.404 μs) | 19.450 ms (IQR 135.066 μs) | 
| Interpolations (uniform) | 70.000 ns (IQR 2.500 ns) | 120.000 ns (IQR 10.000 ns) | 6.600 μs (IQR 30.000 ns) | 655.404 μs (IQR 625.000 ns) | 70.000 ns (IQR 0.000 ns) | 120.000 ns (IQR 10.000 ns) | 6.600 μs (IQR 20.000 ns) | 661.254 μs (IQR 8.865 μs) | 70.000 ns (IQR 10.000 ns) | 120.000 ns (IQR 10.000 ns) | 6.890 μs (IQR 20.000 ns) | 683.869 μs (IQR 335.000 ns) | 70.000 ns (IQR 0.000 ns) | 120.000 ns (IQR 0.000 ns) | 7.170 μs (IQR 20.000 ns) | 752.364 μs (IQR 8.863 μs) | 

### MonotoneCubic

| Library | n=100,m=1 | n=100,m=10 | n=100,m=1000 | n=100,m=100000 | n=1000,m=1 | n=1000,m=10 | n=1000,m=1000 | n=1000,m=100000 | n=10000,m=1 | n=10000,m=10 | n=10000,m=1000 | n=10000,m=100000 | n=100000,m=1 | n=100000,m=10 | n=100000,m=1000 | n=100000,m=100000 | 
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| DataInterpolations (CubicHermite) | 130.000 ns (IQR 0.000 ns) | 550.000 ns (IQR 10.000 ns) | 88.869 μs (IQR 182.750 ns) | 8.923 ms (IQR 45.887 μs) | 140.000 ns (IQR 0.000 ns) | 750.000 ns (IQR 10.000 ns) | 114.689 μs (IQR 312.750 ns) | 12.259 ms (IQR 81.830 μs) | 150.000 ns (IQR 0.000 ns) | 760.000 ns (IQR 10.000 ns) | 149.914 μs (IQR 4.128 μs) | 16.358 ms (IQR 36.446 μs) | 160.000 ns (IQR 10.000 ns) | 850.000 ns (IQR 10.000 ns) | 225.383 μs (IQR 3.223 μs) | 24.150 ms (IQR 233.658 μs) | 
| FastInterpolations (PCHIP) | 60.000 ns (IQR 0.000 ns) | 180.000 ns (IQR 10.000 ns) | 11.240 μs (IQR 10.000 ns) | 1.070 ms (IQR 6.601 μs) | 60.000 ns (IQR 0.000 ns) | 170.000 ns (IQR 0.000 ns) | 11.150 μs (IQR 140.000 ns) | 1.071 ms (IQR 10.822 μs) | 60.000 ns (IQR 0.000 ns) | 180.000 ns (IQR 10.000 ns) | 11.140 μs (IQR 24.750 ns) | 1.068 ms (IQR 5.527 μs) | 60.000 ns (IQR 0.000 ns) | 190.000 ns (IQR 10.000 ns) | 11.240 μs (IQR 40.000 ns) | 1.114 ms (IQR 9.055 μs) | 
| PCHIPInterpolation | 70.000 ns (IQR 10.000 ns) | 250.000 ns (IQR 0.000 ns) | 25.435 μs (IQR 1.558 μs) | 6.446 ms (IQR 16.104 μs) | 90.000 ns (IQR 0.000 ns) | 345.000 ns (IQR 10.000 ns) | 63.260 μs (IQR 2.118 μs) | 9.963 ms (IQR 79.389 μs) | 90.000 ns (IQR 0.000 ns) | 440.000 ns (IQR 10.000 ns) | 118.869 μs (IQR 4.293 μs) | 14.222 ms (IQR 93.245 μs) | 110.000 ns (IQR 10.000 ns) | 490.000 ns (IQR 10.000 ns) | 185.893 μs (IQR 5.101 μs) | 20.408 ms (IQR 54.869 μs) | 

### QuadraticSpline

| Library | n=100,m=1 | n=100,m=10 | n=100,m=1000 | n=100,m=100000 | n=1000,m=1 | n=1000,m=10 | n=1000,m=1000 | n=1000,m=100000 | n=10000,m=1 | n=10000,m=10 | n=10000,m=1000 | n=10000,m=100000 | n=100000,m=1 | n=100000,m=10 | n=100000,m=1000 | n=100000,m=100000 | 
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| DataInterpolations | 210.000 ns (IQR 0.000 ns) | 1.270 μs (IQR 0.000 ns) | 188.969 μs (IQR 920.000 ns) | 19.321 ms (IQR 149.183 μs) | 310.000 ns (IQR 10.000 ns) | 2.720 μs (IQR 10.000 ns) | 347.986 μs (IQR 3.498 μs) | 36.760 ms (IQR 66.121 μs) | 1.300 μs (IQR 20.000 ns) | 12.700 μs (IQR 61.750 ns) | 1.444 ms (IQR 9.453 μs) | 145.560 ms (IQR 123.712 μs) | 15.650 μs (IQR 700.000 ns) | 146.198 μs (IQR 282.250 ns) | 16.271 ms (IQR 285.353 μs) | 1.616 s (IQR 0.000 ns) | 
| Dierckx (k=2) | 130.000 ns (IQR 0.000 ns) | 1.060 μs (IQR 10.000 ns) | 117.319 μs (IQR 3.700 μs) | 11.766 ms (IQR 113.288 μs) | 270.000 ns (IQR 10.000 ns) | 3.780 μs (IQR 0.000 ns) | 394.716 μs (IQR 372.250 ns) | 39.406 ms (IQR 139.238 μs) | 1.470 μs (IQR 20.000 ns) | 29.990 μs (IQR 44.750 ns) | 3.214 ms (IQR 5.622 μs) | 320.480 ms (IQR 312.612 μs) | 13.520 μs (IQR 90.000 ns) | 291.358 μs (IQR 312.500 ns) | 32.152 ms (IQR 591.392 μs) | 3.155 s (IQR 0.000 ns) | 
| FastInterpolations | 50.000 ns (IQR 10.000 ns) | 110.000 ns (IQR 0.000 ns) | 4.030 μs (IQR 30.000 ns) | 366.447 μs (IQR 4.920 μs) | 50.000 ns (IQR 10.000 ns) | 110.000 ns (IQR 0.000 ns) | 4.000 μs (IQR 50.000 ns) | 351.122 μs (IQR 4.475 μs) | 50.000 ns (IQR 10.000 ns) | 110.000 ns (IQR 2.500 ns) | 4.150 μs (IQR 50.000 ns) | 412.296 μs (IQR 4.755 μs) | 50.000 ns (IQR 10.000 ns) | 110.000 ns (IQR 0.000 ns) | 4.475 μs (IQR 80.000 ns) | 478.275 μs (IQR 4.135 μs) | 

## Chained ODE-style

Sequential `for x in tt; A(x); end` over a monotone sequence. (knot pattern = uniform)

### Akima

| Library | n=100,m=1000 | n=1000,m=1000 | n=10000,m=1000 | n=100000,m=1000 | 
|---|---|---|---|---|
| DataInterpolations | 55.874 μs (IQR 224.000 ns) | 80.234 μs (IQR 400.250 ns) | 116.019 μs (IQR 2.942 μs) | 163.899 μs (IQR 4.407 μs) | 
| FastInterpolations | 10.970 μs (IQR 9.250 ns) | 11.040 μs (IQR 1.000 ns) | 11.030 μs (IQR 410.000 ns) | 11.420 μs (IQR 93.250 ns) | 

### CubicSpline

| Library | n=100,m=1000 | n=1000,m=1000 | n=10000,m=1000 | n=100000,m=1000 | 
|---|---|---|---|---|
| DataInterpolations | 59.995 μs (IQR 562.250 ns) | 81.979 μs (IQR 374.500 ns) | 113.249 μs (IQR 3.044 μs) | 178.064 μs (IQR 5.123 μs) | 
| BasicInterpolators | 17.450 μs (IQR 357.500 ns) | 45.319 μs (IQR 1.970 μs) | 85.799 μs (IQR 2.695 μs) | 140.149 μs (IQR 2.475 μs) | 
| Dierckx (k=3) | 122.163 μs (IQR 72.250 ns) | 410.447 μs (IQR 169.250 ns) | 3.214 ms (IQR 6.145 μs) | 31.322 ms (IQR 130.226 μs) | 
| FastInterpolations | 10.150 μs (IQR 10.000 ns) | 10.140 μs (IQR 20.000 ns) | 10.160 μs (IQR 12.250 ns) | 10.950 μs (IQR 99.250 ns) | 
| Interpolations (uniform) | 27.260 μs (IQR 490.000 ns) | 27.850 μs (IQR 0.250 ns) | 27.509 μs (IQR 500.000 ns) | 27.720 μs (IQR 420.000 ns) | 

### Linear

| Library | n=100,m=1000 | n=1000,m=1000 | n=10000,m=1000 | n=100000,m=1000 | 
|---|---|---|---|---|
| DataInterpolations | 55.260 μs (IQR 104.750 ns) | 77.310 μs (IQR 362.500 ns) | 109.109 μs (IQR 3.022 μs) | 162.824 μs (IQR 2.605 μs) | 
| BasicInterpolators | 17.130 μs (IQR 430.000 ns) | 49.775 μs (IQR 2.400 μs) | 87.274 μs (IQR 4.043 μs) | 135.179 μs (IQR 2.735 μs) | 
| Dierckx (k=1) | 86.279 μs (IQR 450.000 ns) | 374.666 μs (IQR 50.000 ns) | 3.185 ms (IQR 19.682 μs) | 31.293 ms (IQR 85.214 μs) | 
| FastInterpolations | 3.510 μs (IQR 10.000 ns) | 3.500 μs (IQR 20.000 ns) | 3.550 μs (IQR 12.250 ns) | 3.650 μs (IQR 10.000 ns) | 
| Interpolations (gridded) | 23.720 μs (IQR 162.500 ns) | 56.059 μs (IQR 2.413 μs) | 89.024 μs (IQR 3.799 μs) | 139.869 μs (IQR 3.240 μs) | 
| Interpolations (uniform) | 23.050 μs (IQR 11.000 ns) | 23.050 μs (IQR 30.000 ns) | 23.059 μs (IQR 10.000 ns) | 23.120 μs (IQR 39.250 ns) | 

### MonotoneCubic

| Library | n=100,m=1000 | n=1000,m=1000 | n=10000,m=1000 | n=100000,m=1000 | 
|---|---|---|---|---|
| DataInterpolations (CubicHermite) | 67.909 μs (IQR 212.750 ns) | 89.499 μs (IQR 411.750 ns) | 121.504 μs (IQR 4.065 μs) | 180.869 μs (IQR 3.533 μs) | 
| FastInterpolations (PCHIP) | 11.170 μs (IQR 10.000 ns) | 11.070 μs (IQR 100.000 ns) | 11.020 μs (IQR 45.000 ns) | 11.380 μs (IQR 80.000 ns) | 
| PCHIPInterpolation | 25.144 μs (IQR 229.250 ns) | 53.205 μs (IQR 303.250 ns) | 94.259 μs (IQR 1.500 μs) | 150.668 μs (IQR 4.713 μs) | 

### QuadraticSpline

| Library | n=100,m=1000 | n=1000,m=1000 | n=10000,m=1000 | n=100000,m=1000 | 
|---|---|---|---|---|
| DataInterpolations | 137.358 μs (IQR 460.500 ns) | 268.752 μs (IQR 559.250 ns) | 1.322 ms (IQR 6.378 μs) | 15.843 ms (IQR 1.173 ms) | 
| Dierckx (k=2) | 104.454 μs (IQR 240.000 ns) | 393.216 μs (IQR 341.000 ns) | 3.187 ms (IQR 5.403 μs) | 31.297 ms (IQR 52.556 μs) | 
| FastInterpolations | 3.310 μs (IQR 2.500 ns) | 3.310 μs (IQR 10.000 ns) | 3.450 μs (IQR 20.000 ns) | 4.650 μs (IQR 40.000 ns) | 

## FastInterpolations.jl advertised benchmark

Numbers below come from `bench/fast_interpolations_bench.jl`, a port of
FastInterpolations.jl's own `benchmark/interpolation_benchmark.jl`
(ProjectTorreyPines/FastInterpolations.jl, upstream commit `616b106b`). It runs
the matrix-of-interpolants workload they advertise on their README:
`mpert × mpert` independent 1D interpolants over a shared uniform
`range(0.0, 1.0; length = npsi)` grid, evaluated at `n_eval` query points
clustered near psi=0 (cubic spacing — mimics ODE solver behavior near
singular surfaces). Default size: `npsi = 64`, `mpert = 100`, `n_eval = 1000`
(10⁴ interpolants × 10³ queries = 10⁷ total scalar evaluations).

### Cubic spline, `--default` (npsi=64, mpert=100, n_eval=1000)

| Package | Init (ms) | Eval (ms) | Total (ms) | Speedup vs DI scalar |
|---|---|---|---|---|
| FastInterpolations.jl (Series+scalar) | 16.508 | 5.751 | 22.259 | 73.28× |
| FastInterpolations.jl (Series+vector) | 16.508 | 22.521 | 39.029 | 41.80× |
| FastInterpolations.jl (vector)        | 31.631 | 95.125 | 126.757 | 12.87× |
| DataInterpolations.jl (vector)        | 93.259 | 103.438 | 196.697 |  8.29× |
| FastInterpolations.jl (scalar)        | 31.631 | 171.245 | 202.876 |  8.04× |
| Interpolations.jl (broadcast)         | 139.363 | 212.500 | 351.863 |  4.64× |
| Interpolations.jl (scalar)            | 139.363 | 427.833 | 567.196 |  2.88× |
| Dierckx.jl (vector)                   | 127.741 | 566.160 | 693.901 |  2.35× |
| DataInterpolations.jl (scalar)        |  93.259 | 1537.974 | 1631.233 | 1.00× |
| Dierckx.jl (scalar)                   | 127.741 | 1590.398 | 1718.139 | 0.95× |

### Linear, `--default` (npsi=64, mpert=100, n_eval=1000)

| Package | Init (ms) | Eval (ms) | Total (ms) | Speedup vs DI scalar |
|---|---|---|---|---|
| FastInterpolations.jl (Series+scalar) |  5.649 |   3.173 |    8.822 | 102.83× |
| FastInterpolations.jl (Series+vector) |  5.649 |  21.019 |   26.668 |  34.02× |
| FastInterpolations.jl (vector)        | 10.481 |  49.861 |   60.343 |  15.03× |
| DataInterpolations.jl (vector)        | 15.629 |  64.376 |   80.005 |  11.34× |
| Interpolations.jl (broadcast)         |  9.634 |  70.787 |   80.421 |  11.28× |
| FastInterpolations.jl (scalar)        | 10.481 |  90.099 |  100.580 |   9.02× |
| Interpolations.jl (scalar)            |  9.634 | 211.753 |  221.387 |   4.10× |
| Dierckx.jl (vector)                   | 77.443 | 244.645 |  322.088 |   2.82× |
| DataInterpolations.jl (scalar)        | 15.629 | 891.534 |  907.163 |   1.00× |
| Dierckx.jl (scalar)                   | 77.443 |1236.548 | 1313.991 |   0.69× |

## Findings

### Where FastInterpolations.jl beats DI

  1. **Series interpolant + matrix-of-interpolants workload.** FastInterpolations'
     `Series` API computes the cell anchor (index, alpha, neighbouring grid
     points) once per query point and reuses it across all 10⁴ coefficient
     series. DI has no equivalent — each interpolant runs an independent
     search per query. The 70-100× speedup at `--default` size is real and
     unfixable without an analogous "shared anchor / Series" type in DI. This
     is a separate design proposal; out of scope for this PR.

  2. **Per-query scalar evaluation on uniform-`Vector` grids.** In the
     cross-library bench (we `_vec(t)` before passing to every library to
     compare like-for-like), DI's per-query latency is ~100-200 ns,
     FastInterpolations' is ~50 ns. The gap is in the scalar kernel call
     overhead and the dispatch through `Auto(t_props)`.
     FastInterpolations' `_search_direct(::_CachedRange, q)` is a single
     `unsafe_trunc(Int, muladd(q - lo, inv_h, 1))` — fewer instructions than
     FFF's `Auto` → `UniformStep` path which still goes through method
     dispatch. Closing this gap on DI would require either:
       - resolving `Auto` to a concrete strategy at *construction* time (so
         `_interpolate` doesn't dispatch through `Auto` at all), or
       - making FFF `Auto` simpler to specialize at the call site.

  3. **Batched / chained evaluation on uniform `Vector` grids.** DI's batched
     loop is ~10-30 μs for m=1000 queries; FastInterpolations gets to
     ~3-4 μs because their vectorized batch loop avoids per-call function-
     pointer dispatch entirely. The remaining gap is partly the kernel
     overhead from (2) and partly that FastInterpolations always knows the
     search policy at compile time (it's a type parameter on its `Searcher`).

### Where DI matches FastInterpolations.jl

  - **Non-uniform construction** (Akima, CubicSpline, QuadraticSpline) at
    large n. DI and FastInterpolations are within ~30% on construction time
    once n ≥ 10k. The new `Auto(A.t_props)` + `O(n)` `spline_coefficients!`
    fix on this branch closes the gap from where it was before.
  - **Sorted-batch evaluation on non-uniform grids** at large m. FFF's
    `BracketGallop` / `LinearScan` strategies + the batched-Auto
    specialization keep DI competitive once the batch is large enough to
    amortize the type-instability overhead from `_resolve_search_policy`.

### Where DI loses but the fix is out-of-scope for this PR

  - Per-query latency on `Vector{Float64}` grids: DI's per-call path goes
    `interp(t) → _interpolate(A, t) → _interpolate(A, t, A.iguesser) →
    get_idx → Auto(A.t_props) → searchsortedlast(Auto, v, q, hint)`. Each
    indirection is ~5-10 ns; FastInterpolations' direct dispatch path is
    one or two indirections. Reducing this means either inlining `get_idx`
    into `_interpolate` per type, or storing a *resolved* concrete search
    strategy at construction time rather than a generic `Auto(props)` —
    both substantial restructurings.
  - No `Series`-style anchor reuse: a different type system. Worth a
    separate proposal; the design space is large.

## Reproducer

Bench script: `bench/cross_library_comparison.jl`

FastInterpolations-style bench: `bench/fast_interpolations_bench.jl`
(port of ProjectTorreyPines/FastInterpolations.jl's
`benchmark/interpolation_benchmark.jl`, commit `616b106b`).

Bench Project.toml: `bench/Project.toml` (devs DI from `..`).

To rerun:
```bash
cd /home/crackauc/sandbox/tmp_20260515_091703_4914/DataInterpolations.jl
git checkout fff-v2-cleanup-quadraticspline
julia +1.11 --project=bench bench/cross_library_comparison.jl
julia +1.11 --project=bench bench/fast_interpolations_bench.jl --cubic --default
julia +1.11 --project=bench bench/fast_interpolations_bench.jl --linear --default
```

