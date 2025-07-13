"""
    LinearInterpolation(u, t; extrapolation_left::ExtrapolationType.T = ExtrapolationType.None, 
    extrapolation::ExtrapolationType.T = ExtrapolationType.None, extrapolation_right::ExtrapolationType.T = ExtrapolationType.None, 
    cache_parameters = false)

It is the method of interpolating between the data points using a linear polynomial. For any point, two data points one each side are chosen and connected with a line.
Extrapolation extends the last linear polynomial on each side.

## Arguments

  - `u`: data points.
  - `t`: time points.

## Keyword Arguments

  - `extrapolation`: The extrapolation type applied left and right of the data. Possible options
    are `ExtrapolationType.None` (default), `ExtrapolationType.Constant`, `ExtrapolationType.Linear`
    `ExtrapolationType.Extension`, `ExtrapolationType.Periodic` and `ExtrapolationType.Reflective`.
  - `extrapolation_left`: The extrapolation type applied left of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `extrapolation_right`: The extrapolation type applied right of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `cache_parameters`: precompute parameters at initialization for faster interpolation
    computations. Note: if activated, `u` and `t` should not be modified. Defaults to `false`.
  - `assume_linear_t`: boolean value to specify a faster index lookup behavior for
    evenly-distributed abscissae. Alternatively, a numerical threshold may be specified
    for a test based on the normalized standard deviation of the difference with respect
    to the straight line (see [`looks_linear`](@ref)). Defaults to 1e-2.
"""
struct LinearInterpolation{uType, tType, IType, pType, T} <: AbstractInterpolation{T}
    u::uType
    t::tType
    I::IType
    p::LinearParameterCache{pType}
    extrapolation_left::ExtrapolationType.T
    extrapolation_right::ExtrapolationType.T
    iguesser::Guesser{tType}
    cache_parameters::Bool
    linear_lookup::Bool
    function LinearInterpolation(u, t, I, p, extrapolation_left, extrapolation_right,
            cache_parameters, assume_linear_t)
        linear_lookup = seems_linear(assume_linear_t, t)
        new{typeof(u), typeof(t), typeof(I), typeof(p.slope), eltype(u)}(
            u, t, I, p, extrapolation_left, extrapolation_right,
            Guesser(t), cache_parameters, linear_lookup)
    end
end

function LinearInterpolation(
        u, t; extrapolation::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None, cache_parameters = false, assume_linear_t = 1e-2)
    extrapolation_left,
    extrapolation_right = munge_extrapolation(
        extrapolation, extrapolation_left, extrapolation_right)
    u, t = munge_data(u, t)
    p = LinearParameterCache(u, t, cache_parameters)
    A = LinearInterpolation(
        u, t, nothing, p, extrapolation_left,
        extrapolation_right, cache_parameters, assume_linear_t)
    I = cumulative_integral(A, cache_parameters)
    LinearInterpolation(
        u, t, I, p, extrapolation_left, extrapolation_right,
        cache_parameters, assume_linear_t)
end

"""
    QuadraticInterpolation(u, t, mode = :Forward; extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation::ExtrapolationType.T = ExtrapolationType.None, extrapolation_right::ExtrapolationType.T = ExtrapolationType.None, 
        cache_parameters = false)

It is the method of interpolating between the data points using quadratic polynomials. For any point, three data points nearby are taken to fit a quadratic polynomial.
Extrapolation extends the last quadratic polynomial on each side.

## Arguments

  - `u`: data points.
  - `t`: time points.
  - `mode`: `:Forward` or `:Backward`. If `:Forward`, two data points ahead of the point and one data point behind is taken for interpolation. If `:Backward`, two data points behind and one ahead is taken for interpolation.

## Keyword Arguments

  - `extrapolation`: The extrapolation type applied left and right of the data. Possible options
    are `ExtrapolationType.None` (default), `ExtrapolationType.Constant`, `ExtrapolationType.Linear`
    `ExtrapolationType.Extension`, `ExtrapolationType.Periodic` and `ExtrapolationType.Reflective`.
  - `extrapolation_left`: The extrapolation type applied left of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `extrapolation_right`: The extrapolation type applied right of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `cache_parameters`: precompute parameters at initialization for faster interpolation computations. Note: if activated, `u` and `t` should not be modified. Defaults to `false`.
  - `assume_linear_t`: boolean value to specify a faster index lookup behaviour for
    evenly-distributed abscissae. Alternatively, a numerical threshold may be specified
    for a test based on the normalized standard deviation of the difference with respect
    to the straight line (see [`looks_linear`](@ref)). Defaults to 1e-2.
"""
struct QuadraticInterpolation{uType, tType, IType, pType, T} <:
       AbstractInterpolation{T}
    u::uType
    t::tType
    I::IType
    p::QuadraticParameterCache{pType}
    mode::Symbol
    extrapolation_left::ExtrapolationType.T
    extrapolation_right::ExtrapolationType.T
    iguesser::Guesser{tType}
    cache_parameters::Bool
    linear_lookup::Bool
    function QuadraticInterpolation(
            u, t, I, p, mode, extrapolation_left,
            extrapolation_right, cache_parameters, assume_linear_t)
        mode ∈ (:Forward, :Backward) ||
            error("mode should be :Forward or :Backward for QuadraticInterpolation")
        linear_lookup = seems_linear(assume_linear_t, t)
        new{typeof(u), typeof(t), typeof(I), typeof(p.α), eltype(u)}(
            u, t, I, p, mode, extrapolation_left, extrapolation_right,
            Guesser(t), cache_parameters, linear_lookup)
    end
end

function QuadraticInterpolation(
        u, t, mode; extrapolation::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None, cache_parameters = false, assume_linear_t = 1e-2)
    extrapolation_left,
    extrapolation_right = munge_extrapolation(
        extrapolation, extrapolation_left, extrapolation_right)
    u, t = munge_data(u, t)
    linear_lookup = seems_linear(assume_linear_t, t)
    p = QuadraticParameterCache(u, t, cache_parameters, mode)
    A = QuadraticInterpolation(
        u, t, nothing, p, mode, extrapolation_left,
        extrapolation_right, cache_parameters, linear_lookup)
    I = cumulative_integral(A, cache_parameters)
    QuadraticInterpolation(u, t, I, p, mode, extrapolation_left,
        extrapolation_right, cache_parameters, linear_lookup)
end

function QuadraticInterpolation(u, t; kwargs...)
    QuadraticInterpolation(u, t, :Forward; kwargs...)
end

"""
    LagrangeInterpolation(u, t, n = length(t) - 1; extrapolation::ExtrapolationType.T = ExtrapolationType.None, 
    extrapolation_left::ExtrapolationType.T = ExtrapolationType.None, extrapolation_right::ExtrapolationType.T = ExtrapolationType.None)

It is the method of interpolation using Lagrange polynomials of (k-1)th order passing through all the data points where k is the number of data points.

## Arguments

  - `u`: data points.
  - `t`: time points.
  - `n`: order of the polynomial. Currently only (k-1)th order where k is the number of data points.

## Keyword Arguments

  - `extrapolation`: The extrapolation type applied left and right of the data. Possible options
    are `ExtrapolationType.None` (default), `ExtrapolationType.Constant`, `ExtrapolationType.Linear`
    `ExtrapolationType.Extension`, `ExtrapolationType.Periodic` and `ExtrapolationType.Reflective`.
  - `extrapolation_left`: The extrapolation type applied left of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `extrapolation_right`: The extrapolation type applied right of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
"""
struct LagrangeInterpolation{uType, tType, T, bcacheType} <:
       AbstractInterpolation{T}
    u::uType
    t::tType
    n::Int
    bcache::bcacheType
    idxs::Vector{Int}
    extrapolation_left::ExtrapolationType.T
    extrapolation_right::ExtrapolationType.T
    iguesser::Guesser{tType}
    function LagrangeInterpolation(u, t, n, extrapolation_left, extrapolation_right)
        bcache = zeros(eltype(u[1]), n + 1)
        idxs = zeros(Int, n + 1)
        fill!(bcache, NaN)
        new{typeof(u), typeof(t), eltype(u), typeof(bcache)}(u,
            t,
            n,
            bcache,
            idxs,
            extrapolation_left,
            extrapolation_right,
            Guesser(t)
        )
    end
end

function LagrangeInterpolation(
        u, t, n = length(t) - 1;
        extrapolation::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None)
    extrapolation_left,
    extrapolation_right = munge_extrapolation(
        extrapolation, extrapolation_left, extrapolation_right)
    u, t = munge_data(u, t)
    if n != length(t) - 1
        error("Currently only n=length(t) - 1 is supported")
    end
    LagrangeInterpolation(u, t, n, extrapolation_left, extrapolation_right)
end

"""
    AkimaInterpolation(u, t; extrapolation::ExtrapolationType.T = ExtrapolationType.None, extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None, cache_parameters = false)

It is a spline interpolation built from cubic polynomials. It forms a continuously differentiable function. For more details, refer: [https://en.wikipedia.org/wiki/Akima_spline](https://en.wikipedia.org/wiki/Akima_spline).
Extrapolation extends the last cubic polynomial on each side.

## Arguments

  - `u`: data points.
  - `t`: time points.

## Keyword Arguments

  - `extrapolation`: The extrapolation type applied left and right of the data. Possible options
    are `ExtrapolationType.None` (default), `ExtrapolationType.Constant`, `ExtrapolationType.Linear`
    `ExtrapolationType.Extension`, `ExtrapolationType.Periodic` and `ExtrapolationType.Reflective`.
  - `extrapolation_left`: The extrapolation type applied left of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `extrapolation_right`: The extrapolation type applied right of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `cache_parameters`: precompute parameters at initialization for faster interpolation computations. Note: if activated, `u` and `t` should not be modified. Defaults to `false`.
  - `assume_linear_t`: boolean value to specify a faster index lookup behaviour for
    evenly-distributed abscissae. Alternatively, a numerical threshold may be specified
    for a test based on the normalized standard deviation of the difference with respect
    to the straight line (see [`looks_linear`](@ref)). Defaults to 1e-2.
"""
struct AkimaInterpolation{uType, tType, IType, bType, cType, dType, T} <:
       AbstractInterpolation{T}
    u::uType
    t::tType
    I::IType
    b::bType
    c::cType
    d::dType
    extrapolation_left::ExtrapolationType.T
    extrapolation_right::ExtrapolationType.T
    iguesser::Guesser{tType}
    cache_parameters::Bool
    linear_lookup::Bool
    function AkimaInterpolation(
            u, t, I, b, c, d, extrapolation_left,
            extrapolation_right, cache_parameters, assume_linear_t)
        linear_lookup = seems_linear(assume_linear_t, t)
        new{typeof(u), typeof(t), typeof(I), typeof(b), typeof(c),
            typeof(d), eltype(u)}(u,
            t,
            I,
            b,
            c,
            d,
            extrapolation_left,
            extrapolation_right,
            Guesser(t),
            cache_parameters,
            linear_lookup
        )
    end
end

function AkimaInterpolation(
        u, t; extrapolation::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None, cache_parameters = false, assume_linear_t = 1e-2)
    extrapolation_left,
    extrapolation_right = munge_extrapolation(
        extrapolation, extrapolation_left, extrapolation_right)
    u, t = munge_data(u, t)
    linear_lookup = seems_linear(assume_linear_t, t)
    n = length(t)
    dt = diff(t)
    m = Array{eltype(u)}(undef, n + 3)
    m[3:(end - 2)] = diff(u) ./ dt
    m[2] = 2m[3] - m[4]
    m[1] = 2m[2] - m[3]
    m[end - 1] = 2m[end - 2] - m[end - 3]
    m[end] = 2m[end - 1] - m[end - 2]

    b = (m[4:end] .+ m[1:(end - 3)]) ./ 2
    dm = abs.(diff(m))
    f1 = dm[3:(n + 2)]
    f2 = dm[1:n]
    f12 = f1 + f2
    ind = findall(f12 .> 1e-9 * maximum(f12))
    b[ind] = (f1[ind] .* m[ind .+ 1] .+
              f2[ind] .* m[ind .+ 2]) ./ f12[ind]
    c = (3 .* m[3:(end - 2)] .- 2 .* b[1:(end - 1)] .- b[2:end]) ./ dt
    d = (b[1:(end - 1)] .+ b[2:end] .- 2 .* m[3:(end - 2)]) ./ dt .^ 2

    A = AkimaInterpolation(
        u, t, nothing, b, c, d, extrapolation_left,
        extrapolation_right, cache_parameters, linear_lookup)
    I = cumulative_integral(A, cache_parameters)
    AkimaInterpolation(u, t, I, b, c, d, extrapolation_left,
        extrapolation_right, cache_parameters, linear_lookup)
end

"""
    ConstantInterpolation(u, t; dir = :left, extrapolation::ExtrapolationType.T = ExtrapolationType.None, extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None, cache_parameters = false)

It is the method of interpolating using a constant polynomial. For any point, two adjacent data points are found on either side (left and right). The value at that point depends on `dir`.
If it is `:left`, then the value at the left point is chosen and if it is `:right`, the value at the right point is chosen.
Extrapolation extends the last constant polynomial at the end points on each side.

## Arguments

  - `u`: data points.
  - `t`: time points.

## Keyword Arguments

  - `dir`: indicates which value should be used for interpolation (`:left` or `:right`).
  - `extrapolation`: The extrapolation type applied left and right of the data. Possible options
    are `ExtrapolationType.None` (default), `ExtrapolationType.Constant`, `ExtrapolationType.Linear`
    `ExtrapolationType.Extension`, `ExtrapolationType.Periodic` and `ExtrapolationType.Reflective`.
  - `extrapolation_left`: The extrapolation type applied left of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `extrapolation_right`: The extrapolation type applied right of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `cache_parameters`: precompute parameters at initialization for faster interpolation computations. Note: if activated, `u` and `t` should not be modified. Defaults to `false`.
  - `assume_linear_t`: boolean value to specify a faster index lookup behaviour for
    evenly-distributed abscissae. Alternatively, a numerical threshold may be specified
    for a test based on the normalized standard deviation of the difference with respect
    to the straight line (see [`looks_linear`](@ref)). Defaults to 1e-2.
"""
struct ConstantInterpolation{uType, tType, IType, T} <: AbstractInterpolation{T}
    u::uType
    t::tType
    I::IType
    p::Nothing
    dir::Symbol # indicates if value to the $dir should be used for the interpolation
    extrapolation_left::ExtrapolationType.T
    extrapolation_right::ExtrapolationType.T
    iguesser::Guesser{tType}
    cache_parameters::Bool
    linear_lookup::Bool
    function ConstantInterpolation(
            u, t, I, dir, extrapolation_left, extrapolation_right,
            cache_parameters, assume_linear_t)
        linear_lookup = seems_linear(assume_linear_t, t)
        new{typeof(u), typeof(t), typeof(I), eltype(u)}(
            u, t, I, nothing, dir, extrapolation_left, extrapolation_right,
            Guesser(t), cache_parameters, linear_lookup)
    end
end

function ConstantInterpolation(
        u, t; dir = :left, extrapolation::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None,
        cache_parameters = false, assume_linear_t = 1e-2)
    extrapolation_left,
    extrapolation_right = munge_extrapolation(
        extrapolation, extrapolation_left, extrapolation_right)
    u, t = munge_data(u, t)
    A = ConstantInterpolation(
        u, t, nothing, dir, extrapolation_left,
        extrapolation_right, cache_parameters, assume_linear_t)
    I = cumulative_integral(A, cache_parameters)
    ConstantInterpolation(u, t, I, dir, extrapolation_left, extrapolation_right,
        cache_parameters, assume_linear_t)
end

"""
    SmoothedConstantInterpolation(u, t; d_max = Inf, extrapolate = false,
        cache_parameters = false, assume_linear_t = 1e-2)

It is a method for interpolating constantly with forward fill, with smoothing around the
value transitions to make the curve continuously differentiable while the integral never
drifts far from the integral of constant interpolation. `u[end]` is ignored,
except when using extrapolation types `Constant` or `Extension`.

## Arguments

  - `u`: data points.
  - `t`: time points.

## Keyword Arguments

  - `d_max`: Around each time point `tᵢ` there is a continuously differentiable (quadratic) transition between `uᵢ₋₁` and `uᵢ`,
    on the interval `[tᵢ - d, tᵢ + d]`. The distance `d` is determined as `d = min((tᵢ - tᵢ₋₁)/2, (tᵢ₊₁ - tᵢ)/2, d_max)`.
  - `extrapolation`: The extrapolation type applied left and right of the data. Possible options
    are `ExtrapolationType.None` (default), `ExtrapolationType.Constant`, `ExtrapolationType.Linear`
    `ExtrapolationType.Extension`, `ExtrapolationType.Periodic` (also made smooth at the boundaries) and `ExtrapolationType.Reflective`.
  - `extrapolation_left`: The extrapolation type applied left of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `extrapolation_right`: The extrapolation type applied right of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `cache_parameters`: precompute parameters at initialization for faster interpolation computations. Note: if activated, `u` and `t` should not be modified. Defaults to `false`.
  - `assume_linear_t`: boolean value to specify a faster index lookup behavior for
    evenly-distributed abscissae. Alternatively, a numerical threshold may be specified
    for a test based on the normalized standard deviation of the difference with respect
    to the straight line (see [`looks_linear`](@ref)). Defaults to 1e-2.
"""
struct SmoothedConstantInterpolation{uType, tType, IType, dType, cType, dmaxType, T} <:
       AbstractInterpolation{T}
    u::uType
    t::tType
    I::IType
    p::SmoothedConstantParameterCache{dType, cType}
    d_max::dmaxType
    extrapolation_left::ExtrapolationType.T
    extrapolation_right::ExtrapolationType.T
    iguesser::Guesser{tType}
    cache_parameters::Bool
    linear_lookup::Bool
    function SmoothedConstantInterpolation(
            u, t, I, p, d_max, extrapolation_left,
            extrapolation_right, cache_parameters, assume_linear_t)
        linear_lookup = seems_linear(assume_linear_t, t)
        new{typeof(u), typeof(t), typeof(I), typeof(p.d),
            typeof(p.c), typeof(d_max), eltype(u)}(
            u, t, I, p, d_max, extrapolation_left, extrapolation_right,
            Guesser(t), cache_parameters, linear_lookup)
    end
end

function SmoothedConstantInterpolation(
        u, t; d_max = Inf, extrapolation::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None,
        cache_parameters = false, assume_linear_t = 1e-2)
    extrapolation_left,
    extrapolation_right = munge_extrapolation(
        extrapolation, extrapolation_left, extrapolation_right)
    u, t = munge_data(u, t)
    p = SmoothedConstantParameterCache(
        u, t, cache_parameters, d_max, extrapolation_left, extrapolation_right)
    A = SmoothedConstantInterpolation(
        u, t, nothing, p, d_max, extrapolation_left,
        extrapolation_right, cache_parameters, assume_linear_t)
    I = cumulative_integral(A, cache_parameters)
    SmoothedConstantInterpolation(
        u, t, I, p, d_max, extrapolation_left,
        extrapolation_right, cache_parameters, assume_linear_t)
end

"""
    QuadraticSpline(u, t; extrapolation::ExtrapolationType.T = ExtrapolationType.None, extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None, cache_parameters = false)

It is a spline interpolation using piecewise quadratic polynomials between each pair of data points. Its first derivative is also continuous.
Extrapolation extends the last quadratic polynomial on each side.

## Arguments

  - `u`: data points.
  - `t`: time points.

## Keyword Arguments

  - `extrapolation`: The extrapolation type applied left and right of the data. Possible options
    are `ExtrapolationType.None` (default), `ExtrapolationType.Constant`, `ExtrapolationType.Linear`
    `ExtrapolationType.Extension`, `ExtrapolationType.Periodic` and `ExtrapolationType.Reflective`.
  - `extrapolation_left`: The extrapolation type applied left of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `extrapolation_right`: The extrapolation type applied right of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `cache_parameters`: precompute parameters at initialization for faster interpolation computations. Note: if activated, `u` and `t` should not be modified. Defaults to `false`.
  - `assume_linear_t`: boolean value to specify a faster index lookup behaviour for
    evenly-distributed abscissae. Alternatively, a numerical threshold may be specified
    for a test based on the normalized standard deviation of the difference with respect
    to the straight line (see [`looks_linear`](@ref)). Defaults to 1e-2.
"""
struct QuadraticSpline{uType, tType, IType, pType, kType, cType, scType, T} <:
       AbstractInterpolation{T}
    u::uType
    t::tType
    I::IType
    p::QuadraticSplineParameterCache{pType}
    k::kType # knot vector
    c::cType # B-spline control points
    sc::scType # Spline coefficients (preallocated memory)
    extrapolation_left::ExtrapolationType.T
    extrapolation_right::ExtrapolationType.T
    iguesser::Guesser{tType}
    cache_parameters::Bool
    linear_lookup::Bool
    function QuadraticSpline(
            u, t, I, p, k, c, sc, extrapolation_left,
            extrapolation_right, cache_parameters, assume_linear_t)
        linear_lookup = seems_linear(assume_linear_t, t)
        new{typeof(u), typeof(t), typeof(I), typeof(p.α), typeof(k),
            typeof(c), typeof(sc), eltype(u)}(u,
            t,
            I,
            p,
            k,
            c,
            sc,
            extrapolation_left,
            extrapolation_right,
            Guesser(t),
            cache_parameters,
            linear_lookup
        )
    end
end

function QuadraticSpline(
        u::uType, t; extrapolation::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None,
        cache_parameters = false, assume_linear_t = 1e-2) where {uType <:
                                                                 AbstractVector{<:Number}}
    extrapolation_left,
    extrapolation_right = munge_extrapolation(
        extrapolation, extrapolation_left, extrapolation_right)
    u, t = munge_data(u, t)

    n = length(t)
    dtype_sc = typeof(t[1] / t[1])
    sc = zeros(dtype_sc, n)
    k, A = quadratic_spline_params(t, sc)
    c = A \ u

    p = QuadraticSplineParameterCache(u, t, k, c, sc, cache_parameters)
    A = QuadraticSpline(
        u, t, nothing, p, k, c, sc, extrapolation_left,
        extrapolation_right, cache_parameters, assume_linear_t)
    I = cumulative_integral(A, cache_parameters)
    QuadraticSpline(u, t, I, p, k, c, sc, extrapolation_left,
        extrapolation_right, cache_parameters, assume_linear_t)
end

function QuadraticSpline(
        u::uType, t; extrapolation::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None, cache_parameters = false,
        assume_linear_t = 1e-2) where {uType <:
                                       AbstractVector}
    extrapolation_left,
    extrapolation_right = munge_extrapolation(
        extrapolation, extrapolation_left, extrapolation_right)
    u, t = munge_data(u, t)

    n = length(t)
    dtype_sc = typeof(t[1] / t[1])
    sc = zeros(dtype_sc, n)
    k, A = quadratic_spline_params(t, sc)

    eltype_c_prototype = one(dtype_sc) * first(u)
    c = [similar(eltype_c_prototype) for _ in 1:n]

    # Assuming u contains arrays of equal shape
    for j in eachindex(eltype_c_prototype)
        c_dim = A \ [u_[j] for u_ in u]
        for (i, c_dim_) in enumerate(c_dim)
            c[i][j] = c_dim_
        end
    end

    p = QuadraticSplineParameterCache(u, t, k, c, sc, cache_parameters)
    A = QuadraticSpline(
        u, t, nothing, p, k, c, sc, extrapolation_left,
        extrapolation_right, cache_parameters, assume_linear_t)
    I = cumulative_integral(A, cache_parameters)
    QuadraticSpline(u, t, I, p, k, c, sc, extrapolation_left,
        extrapolation_right, cache_parameters, assume_linear_t)
end

"""
    CubicSpline(u, t; extrapolation::ExtrapolationType.T = ExtrapolationType.None, extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None, cache_parameters = false)

It is a spline interpolation using piecewise cubic polynomials between each pair of data points. Its first and second derivative is also continuous.
Second derivative on both ends are zero, which are also called "natural" boundary conditions. Extrapolation extends the last cubic polynomial on each side.

## Arguments

  - `u`: data points.
  - `t`: time points.

## Keyword Arguments

  - `extrapolation`: The extrapolation type applied left and right of the data. Possible options
    are `ExtrapolationType.None` (default), `ExtrapolationType.Constant`, `ExtrapolationType.Linear`
    `ExtrapolationType.Extension`, `ExtrapolationType.Periodic` and `ExtrapolationType.Reflective`.
  - `extrapolation_left`: The extrapolation type applied left of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `extrapolation_right`: The extrapolation type applied right of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `cache_parameters`: precompute parameters at initialization for faster interpolation computations. Note: if activated, `u` and `t` should not be modified. Defaults to `false`.
  - `assume_linear_t`: boolean value to specify a faster index lookup behaviour for
    evenly-distributed abscissae. Alternatively, a numerical threshold may be specified
    for a test based on the normalized standard deviation of the difference with respect
    to the straight line (see [`looks_linear`](@ref)). Defaults to 1e-2.
"""
struct CubicSpline{uType, tType, IType, pType, hType, zType, T} <:
       AbstractInterpolation{T}
    u::uType
    t::tType
    I::IType
    p::CubicSplineParameterCache{pType}
    h::hType
    z::zType
    extrapolation_left::ExtrapolationType.T
    extrapolation_right::ExtrapolationType.T
    iguesser::Guesser{tType}
    cache_parameters::Bool
    linear_lookup::Bool
    function CubicSpline(u, t, I, p, h, z, extrapolation_left,
            extrapolation_right, cache_parameters, assume_linear_t)
        linear_lookup = seems_linear(assume_linear_t, t)
        new{typeof(u), typeof(t), typeof(I), typeof(p.c₁),
            typeof(h), typeof(z), eltype(u)}(
            u,
            t,
            I,
            p,
            h,
            z,
            extrapolation_left,
            extrapolation_right,
            Guesser(t),
            cache_parameters,
            linear_lookup
        )
    end
end

function CubicSpline(u::AbstractVector{<:Number},
        t; extrapolation::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None, cache_parameters = false,
        assume_linear_t = 1e-2)
    extrapolation_left,
    extrapolation_right = munge_extrapolation(
        extrapolation, extrapolation_left, extrapolation_right)
    u, t = munge_data(u, t)
    n = length(t) - 1
    h = vcat(0, map(k -> t[k + 1] - t[k], 1:(length(t) - 1)), 0)
    dl = vcat(h[2:n], zero(eltype(h)))
    d_tmp = 2 .* (h[1:(n + 1)] .+ h[2:(n + 2)])
    du = vcat(zero(eltype(h)), h[3:(n + 1)])
    tA = Tridiagonal(dl, d_tmp, du)

    # zero for element type of d, which we don't know yet
    typed_zero = zero(6(u[begin + 2] - u[begin + 1]) / h[begin + 2] -
                      6(u[begin + 1] - u[begin]) / h[begin + 1])

    d = map(
        i -> i == 1 || i == n + 1 ? typed_zero :
             6(u[i + 1] - u[i]) / h[i + 1] - 6(u[i] - u[i - 1]) / h[i],
        1:(n + 1))
    z = tA \ d
    linear_lookup = seems_linear(assume_linear_t, t)
    p = CubicSplineParameterCache(u, h, z, cache_parameters)
    A = CubicSpline(
        u, t, nothing, p, h[1:(n + 1)], z, extrapolation_left,
        extrapolation_right, cache_parameters, linear_lookup)
    I = cumulative_integral(A, cache_parameters)
    CubicSpline(u, t, I, p, h[1:(n + 1)], z, extrapolation_left,
        extrapolation_right, cache_parameters, linear_lookup)
end

function CubicSpline(u::AbstractArray{T, N},
        t;
        extrapolation::ExtrapolationType.T = ExtrapolationType.None, extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None, cache_parameters = false,
        assume_linear_t = 1e-2) where {T, N}
    extrapolation_left,
    extrapolation_right = munge_extrapolation(
        extrapolation, extrapolation_left, extrapolation_right)
    u, t = munge_data(u, t)
    n = length(t) - 1
    h = vcat(0, map(k -> t[k + 1] - t[k], 1:(length(t) - 1)), 0)
    dl = vcat(h[2:n], zero(eltype(h)))
    d_tmp = 2 .* (h[1:(n + 1)] .+ h[2:(n + 2)])
    du = vcat(zero(eltype(h)), h[3:(n + 1)])
    tA = Tridiagonal(dl, d_tmp, du)

    # zero for element type of d, which we don't know yet
    ax = axes(u)[1:(end - 1)]
    typed_zero = zero(6(u[ax..., begin + 2] - u[ax..., begin + 1]) / h[begin + 2] -
                      6(u[ax..., begin + 1] - u[ax..., begin]) / h[begin + 1])

    h_ = reshape(h, ones(Int64, N - 1)..., :)
    ax_h = axes(h_)[1:(end - 1)]
    d = 6 * ((u[ax..., 3:(n + 1)] - u[ax..., 2:n]) ./ h_[ax_h..., 3:(n + 1)]) -
        6 * ((u[ax..., 2:n] - u[ax..., 1:(n - 1)]) ./ h_[ax_h..., 2:n])
    d = cat(typed_zero, d, typed_zero; dims = ndims(d))
    d_reshaped = reshape(d, prod(size(d)[1:(end - 1)]), :)
    z = (tA \ d_reshaped')'
    z = reshape(z, size(u)...)
    linear_lookup = seems_linear(assume_linear_t, t)
    p = CubicSplineParameterCache(u, h, z, cache_parameters)
    A = CubicSpline(
        u, t, nothing, p, h[1:(n + 1)], z, extrapolation_left,
        extrapolation_right, cache_parameters, linear_lookup)
    I = cumulative_integral(A, cache_parameters)
    CubicSpline(u, t, I, p, h[1:(n + 1)], z, extrapolation_left,
        extrapolation_right, cache_parameters, linear_lookup)
end

function CubicSpline(
        u::AbstractVector, t; extrapolation::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None, cache_parameters = false,
        assume_linear_t = 1e-2)
    extrapolation_left,
    extrapolation_right = munge_extrapolation(
        extrapolation, extrapolation_left, extrapolation_right)
    u, t = munge_data(u, t)
    n = length(t) - 1
    h = vcat(0, map(k -> t[k + 1] - t[k], 1:(length(t) - 1)), 0)
    dl = vcat(h[2:n], zero(eltype(h)))
    d_tmp = 2 .* (h[1:(n + 1)] .+ h[2:(n + 2)])
    du = vcat(zero(eltype(h)), h[3:(n + 1)])
    tA = Tridiagonal(dl, d_tmp, du)
    d_ = map(
        i -> i == 1 || i == n + 1 ? zeros(eltype(t), size(u[1])) :
             6(u[i + 1] - u[i]) / h[i + 1] - 6(u[i] - u[i - 1]) / h[i],
        1:(n + 1))
    d = transpose(reshape(reduce(hcat, d_), :, n + 1))
    z_ = reshape(transpose(tA \ d), size(u[1])..., :)
    z = [z_s for z_s in eachslice(z_, dims = ndims(z_))]

    p = CubicSplineParameterCache(u, h, z, cache_parameters)
    A = CubicSpline(
        u, t, nothing, p, h[1:(n + 1)], z, extrapolation_left,
        extrapolation_right, cache_parameters, assume_linear_t)
    I = cumulative_integral(A, cache_parameters)
    CubicSpline(u, t, I, p, h[1:(n + 1)], z, extrapolation_left,
        extrapolation_right, cache_parameters, assume_linear_t)
end

"""
    BSplineInterpolation(u, t, d, pVecType, knotVecType; extrapolation::ExtrapolationType.T = ExtrapolationType.None, extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None)

It is a curve defined by the linear combination of `n` basis functions of degree `d` where `n` is the number of data points. For more information, refer [https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-curve.html](https://pages.mtu.edu/%7Eshene/COURSES/cs3621/NOTES/spline/B-spline/bspline-curve.html).
Extrapolation is a constant polynomial of the end points on each side.

## Arguments

  - `u`: data points.
  - `t`: time points.
  - `d`: degree of the piecewise polynomial.
  - `pVecType`: symbol to parameters vector, `:Uniform` for uniform spaced parameters and `:ArcLen` for parameters generated by chord length method.
  - `knotVecType`: symbol to knot vector, `:Uniform` for uniform knot vector, `:Average` for average spaced knot vector.

## Keyword Arguments

  - `extrapolation`: The extrapolation type applied left and right of the data. Possible options
    are `ExtrapolationType.None` (default), `ExtrapolationType.Constant`, `ExtrapolationType.Linear`
    `ExtrapolationType.Extension`, `ExtrapolationType.Periodic` and `ExtrapolationType.Reflective`.
  - `extrapolation_left`: The extrapolation type applied left of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `extrapolation_right`: The extrapolation type applied right of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `assume_linear_t`: boolean value to specify a faster index lookup behavior for
    evenly-distributed abscissae. Alternatively, a numerical threshold may be specified
    for a test based on the normalized standard deviation of the difference with respect
    to the straight line (see [`looks_linear`](@ref)). Defaults to 1e-2.
"""
struct BSplineInterpolation{uType, tType, pType, kType, cType, scType, T} <:
       AbstractInterpolation{T}
    u::uType
    t::tType
    d::Int    # degree
    p::pType  # params vector
    k::kType  # knot vector
    c::cType  # control points
    sc::scType  # Spline coefficients (preallocated memory)
    pVecType::Symbol
    knotVecType::Symbol
    extrapolation_left::ExtrapolationType.T
    extrapolation_right::ExtrapolationType.T
    iguesser::Guesser{tType}
    linear_lookup::Bool
    function BSplineInterpolation(u,
            t,
            d,
            p,
            k,
            c,
            sc,
            pVecType,
            knotVecType,
            extrapolation_left,
            extrapolation_right,
            assume_linear_t)
        linear_lookup = seems_linear(assume_linear_t, t)
        new{typeof(u), typeof(t), typeof(p), typeof(k), typeof(c), typeof(sc), eltype(u)}(
            u,
            t,
            d,
            p,
            k,
            c,
            sc,
            pVecType,
            knotVecType,
            extrapolation_left,
            extrapolation_right,
            Guesser(t),
            linear_lookup
        )
    end
end

function BSplineInterpolation(
        u::AbstractVector, t, d, pVecType, knotVecType;
        extrapolation::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None, assume_linear_t = 1e-2)
    extrapolation_left,
    extrapolation_right = munge_extrapolation(
        extrapolation, extrapolation_left, extrapolation_right)
    u, t = munge_data(u, t)
    n = length(t)
    n < d + 1 && error("BSplineInterpolation needs at least d + 1, i.e. $(d+1) points.")
    
    # Initialize parameter vector
    param_vec = zero(t)
    param_vec[1] = zero(eltype(t))
    param_vec[end] = one(eltype(t))
    
    # Initialize knot vector
    knot_vec = zeros(eltype(t), n + d + 1)

    # Compute parameter vector based on type
    if pVecType == :Uniform
        for i in 2:(n - 1)
            param_vec[i] = param_vec[1] + (i - 1) * (param_vec[end] - param_vec[1]) / (n - 1)
        end
    elseif pVecType == :ArcLen
        # Only compute arc lengths when needed, using diff and broadcasting
        distances = sqrt.((diff(t)).^2 .+ (diff(u)).^2)
        cumulative_distances = cumsum(distances)
        total_distance = cumulative_distances[end]
        for i in 2:(n - 1)
            param_vec[i] = param_vec[1] + cumulative_distances[i - 1] / total_distance * (param_vec[end] - param_vec[1])
        end
    end

    # Set boundary knots using vectorized assignment
    knot_vec[1:d+1] .= param_vec[1]
    knot_vec[end-d:end] .= param_vec[end]

    # Compute cumulative parameter sum for internal knots using cumsum
    param_cumsum = cumsum(param_vec[2:end-1])

    if knotVecType == :Uniform
        # uniformly spaced knot vector
        # this method is not recommended because, if it is used with the chord length method for global interpolation,
        # the system of linear equations would be singular.
        for i in (d + 2):n
            knot_vec[i] = knot_vec[1] + (i - d - 1) // (n - d) * (knot_vec[end] - knot_vec[1])
        end
    elseif knotVecType == :Average
        # average spaced knot vector
        index = 1
        if d + 2 <= n
            knot_vec[d + 2] = 1 // d * param_cumsum[d]
        end
        for i in (d + 3):n
            knot_vec[i] = 1 // d * (param_cumsum[index + d] - param_cumsum[index])
            index += 1
        end
    end
    
    # control points
    spline_coeffs = zeros(eltype(t), n, n)
    spline_coefficients!(spline_coeffs, d, knot_vec, param_vec)
    control_points = vec(spline_coeffs \ u[:, :])
    spline_coeffs = zeros(eltype(t), n)
    BSplineInterpolation(
        u, t, d, param_vec, knot_vec, control_points, spline_coeffs, pVecType, knotVecType,
        extrapolation_left, extrapolation_right, assume_linear_t)
end

function BSplineInterpolation(
        u::AbstractArray, t, d, pVecType, knotVecType;
        extrapolation::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None,
        assume_linear_t = 1e-2)
    extrapolation_left,
    extrapolation_right = munge_extrapolation(
        extrapolation, extrapolation_left, extrapolation_right)
    u, t = munge_data(u, t)
    n = length(t)
    n < d + 1 && error("BSplineInterpolation needs at least d + 1, i.e. $(d+1) points.")
    
    # Initialize parameter vector
    param_vec = zero(t)
    param_vec[1] = zero(eltype(t))
    param_vec[end] = one(eltype(t))
    
    # Initialize knot vector
    knot_vec = zeros(eltype(t), n + d + 1)

    # Compute parameter vector based on type
    array_axes = axes(u)[1:(end - 1)]

    if pVecType == :Uniform
        for i in 2:(n - 1)
            param_vec[i] = param_vec[1] + (i - 1) * (param_vec[end] - param_vec[1]) / (n - 1)
        end
    elseif pVecType == :ArcLen
        # Only compute arc lengths when needed, using diff and broadcasting
        time_diffs = diff(t)
        spatial_diffs = [sqrt(sum((u[array_axes..., i+1] - u[array_axes..., i]).^2)) for i in 1:(n-1)]
        distances = sqrt.(time_diffs.^2 .+ spatial_diffs.^2)
        cumulative_distances = cumsum(distances)
        total_distance = cumulative_distances[end]
        for i in 2:(n - 1)
            param_vec[i] = param_vec[1] + cumulative_distances[i - 1] / total_distance * (param_vec[end] - param_vec[1])
        end
    end

    # Set boundary knots using vectorized assignment
    knot_vec[1:d+1] .= param_vec[1]
    knot_vec[end-d:end] .= param_vec[end]

    # Compute cumulative parameter sum for internal knots using cumsum
    param_cumsum = cumsum(param_vec[2:end-1])

    if knotVecType == :Uniform
        # uniformly spaced knot vector
        # this method is not recommended because, if it is used with the chord length method for global interpolation,
        # the system of linear equations would be singular.
        for i in (d + 2):n
            knot_vec[i] = knot_vec[1] + (i - d - 1) // (n - d) * (knot_vec[end] - knot_vec[1])
        end
    elseif knotVecType == :Average
        # average spaced knot vector
        index = 1
        if d + 2 <= n
            knot_vec[d + 2] = 1 // d * param_cumsum[d]
        end
        for i in (d + 3):n
            knot_vec[i] = 1 // d * (param_cumsum[index + d] - param_cumsum[index])
            index += 1
        end
    end
    
    # control points
    spline_coeffs = zeros(eltype(t), n, n)
    spline_coefficients!(spline_coeffs, d, knot_vec, param_vec)
    control_points = (spline_coeffs \ reshape(u, prod(size(u)[1:(end - 1)]), :)')'
    control_points = reshape(control_points, size(u)...)
    spline_coeffs = zeros(eltype(t), n)
    BSplineInterpolation(
        u, t, d, param_vec, knot_vec, control_points, spline_coeffs, pVecType, knotVecType,
        extrapolation_left, extrapolation_right, assume_linear_t)
end

"""
    BSplineApprox(u, t, d, h, pVecType, knotVecType; extrapolation::ExtrapolationType.T = ExtrapolationType.None, extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None)

It is a regression based B-spline. The argument choices are the same as the `BSplineInterpolation`, with the additional parameter `h < length(t)` which is the number of control points to use, with smaller `h` indicating more smoothing.
For more information, refer [http://www.cad.zju.edu.cn/home/zhx/GM/009/00-bsia.pdf](http://www.cad.zju.edu.cn/home/zhx/GM/009/00-bsia.pdf).
Extrapolation is a constant polynomial of the end points on each side.

## Arguments

  - `u`: data points.
  - `t`: time points.
  - `d`: degree of the piecewise polynomial.
  - `h`: number of control points to use.
  - `pVecType`: symbol to parameters vector, `:Uniform` for uniform spaced parameters and `:ArcLen` for parameters generated by chord length method.
  - `knotVecType`: symbol to knot vector, `:Uniform` for uniform knot vector, `:Average` for average spaced knot vector.

## Keyword Arguments

  - `extrapolation`: The extrapolation type applied left and right of the data. Possible options
    are `ExtrapolationType.None` (default), `ExtrapolationType.Constant`, `ExtrapolationType.Linear`
    `ExtrapolationType.Extension`, `ExtrapolationType.Periodic` and `ExtrapolationType.Reflective`.
  - `extrapolation_left`: The extrapolation type applied left of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `extrapolation_right`: The extrapolation type applied right of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `assume_linear_t`: boolean value to specify a faster index lookup behaviour for
    evenly-distributed abscissae. Alternatively, a numerical threshold may be specified
    for a test based on the normalized standard deviation of the difference with respect
    to the straight line (see [`looks_linear`](@ref)). Defaults to 1e-2.
"""
struct BSplineApprox{uType, tType, pType, kType, cType, scType, T} <:
       AbstractInterpolation{T}
    u::uType
    t::tType
    d::Int    # degree
    h::Int    # number of control points (n => h >= d >= 1)
    p::pType  # params vector
    k::kType  # knot vector
    c::cType  # control points
    sc::scType  # Spline coefficients (preallocated memory)
    pVecType::Symbol
    knotVecType::Symbol
    extrapolation_left::ExtrapolationType.T
    extrapolation_right::ExtrapolationType.T
    iguesser::Guesser{tType}
    linear_lookup::Bool
    function BSplineApprox(u,
            t,
            d,
            h,
            p,
            k,
            c,
            sc,
            pVecType,
            knotVecType,
            extrapolation_left,
            extrapolation_right,
            assume_linear_t
    )
        linear_lookup = seems_linear(assume_linear_t, t)
        new{typeof(u), typeof(t), typeof(p), typeof(k), typeof(c), typeof(sc), eltype(u)}(
            u,
            t,
            d,
            h,
            p,
            k,
            c,
            sc,
            pVecType,
            knotVecType,
            extrapolation_left,
            extrapolation_right,
            Guesser(t),
            linear_lookup
        )
    end
end

function BSplineApprox(
        u::AbstractVector, t, d, h, pVecType, knotVecType;
        extrapolation::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None, assume_linear_t = 1e-2)
    extrapolation_left,
    extrapolation_right = munge_extrapolation(
        extrapolation, extrapolation_left, extrapolation_right)
    u, t = munge_data(u, t)
    n = length(t)
    h < d + 1 && error("BSplineApprox needs at least d + 1, i.e. $(d+1) control points.")
    
    # Initialize parameter vector
    param_vec = zero(t)
    param_vec[1] = zero(eltype(t))
    param_vec[end] = one(eltype(t))
    
    # Initialize knot vector
    knot_vec = zeros(eltype(t), h + d + 1)

    # Compute parameter vector based on type
    if pVecType == :Uniform
        for i in 2:(n - 1)
            param_vec[i] = param_vec[1] + (i - 1) * (param_vec[end] - param_vec[1]) / (n - 1)
        end
    elseif pVecType == :ArcLen
        # Only compute arc lengths when needed, using diff and broadcasting
        distances = sqrt.((diff(t)).^2 .+ (diff(u)).^2)
        cumulative_distances = cumsum(distances)
        total_distance = cumulative_distances[end]
        for i in 2:(n - 1)
            param_vec[i] = param_vec[1] + cumulative_distances[i - 1] / total_distance * (param_vec[end] - param_vec[1])
        end
    end

    # Set boundary knots using vectorized assignment
    knot_vec[1:d+1] .= param_vec[1]
    knot_vec[end-d:end] .= param_vec[end]

    # Compute cumulative parameter sum for internal knots using cumsum
    param_cumsum = cumsum(param_vec[2:end-1])

    if knotVecType == :Uniform
        # uniformly spaced knot vector
        # this method is not recommended because, if it is used with the chord length method for global interpolation,
        # the system of linear equations would be singular.
        for i in (d + 2):h
            knot_vec[i] = knot_vec[1] + (i - d - 1) // (h - d) * (knot_vec[end] - knot_vec[1])
        end
    elseif knotVecType == :Average
        # average spaced knot vector using improved distribution
        # The goal is to fill k[d+2:h] with h-d-1 values from param_vec[1] to param_vec[end] 
        # in a similar distribution to how the vector param_vec is distributed
        num_internal_knots = h - d - 1
        if num_internal_knots > 0
            # Use LinearInterpolation for better distribution across parameter domain
            param_interp = LinearInterpolation(param_vec, range(0, 1, length = n))
            if num_internal_knots == 1
                # Special case for single knot to avoid range issues
                knot_vec[d + 2] = param_interp(0.5)
            else
                internal_positions = range(0, 1, length = num_internal_knots)
                knot_vec[(d+2):h] .= param_interp.(internal_positions)
            end
        end
    end
    
    # control points
    control_points = zeros(eltype(u), h)
    control_points[1] = u[1]
    control_points[end] = u[end]
    data_residual = zeros(eltype(u), n)
    spline_coeffs = zeros(eltype(t), n, h)
    for i in 1:n
        spline_coefficients!(view(spline_coeffs, i, :), d, knot_vec, param_vec[i])
    end
    for k in 2:(n - 1)
        data_residual[k] = u[k] - spline_coeffs[k, 1] * u[1] - spline_coeffs[k, h] * u[end]
    end
    coeff_matrix = Matrix{eltype(u)}(undef, h - 2, 1)
    for i in 2:(h - 1)
        residual_sum = 0.0
        for k in 2:(n - 1)
            residual_sum += spline_coeffs[k, i] * data_residual[k]
        end
        coeff_matrix[i - 1] = residual_sum
    end
    spline_coeffs = spline_coeffs[2:(end - 1), 2:(h - 1)]
    system_matrix = transpose(spline_coeffs) * spline_coeffs
    solution = system_matrix \ coeff_matrix
    control_points[2:(end - 1)] .= vec(solution)
    spline_coeffs = zeros(eltype(t), h)
    BSplineApprox(
        u, t, d, h, param_vec, knot_vec, control_points, spline_coeffs, pVecType, knotVecType,
        extrapolation_left, extrapolation_right, assume_linear_t)
end

function BSplineApprox(
        u::AbstractArray{T, N}, t, d, h, pVecType, knotVecType;
        extrapolation::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None,
        assume_linear_t = 1e-2) where {T, N}
    extrapolation_left,
    extrapolation_right = munge_extrapolation(
        extrapolation, extrapolation_left, extrapolation_right)
    u, t = munge_data(u, t)
    n = length(t)
    h < d + 1 && error("BSplineApprox needs at least d + 1, i.e. $(d+1) control points.")
    
    # Initialize parameter vector
    param_vec = zero(t)
    param_vec[1] = zero(eltype(t))
    param_vec[end] = one(eltype(t))
    
    # Initialize knot vector
    knot_vec = zeros(eltype(t), h + d + 1)

    # Compute parameter vector based on type
    array_axes = axes(u)[1:(end - 1)]

    if pVecType == :Uniform
        for i in 2:(n - 1)
            param_vec[i] = param_vec[1] + (i - 1) * (param_vec[end] - param_vec[1]) / (n - 1)
        end
    elseif pVecType == :ArcLen
        # Only compute arc lengths when needed, using diff and broadcasting
        time_diffs = diff(t)
        spatial_diffs = [sqrt(sum((u[array_axes..., i+1] - u[array_axes..., i]).^2)) for i in 1:(n-1)]
        distances = sqrt.(time_diffs.^2 .+ spatial_diffs.^2)
        cumulative_distances = cumsum(distances)
        total_distance = cumulative_distances[end]
        for i in 2:(n - 1)
            param_vec[i] = param_vec[1] + cumulative_distances[i - 1] / total_distance * (param_vec[end] - param_vec[1])
        end
    end

    # Set boundary knots using vectorized assignment
    knot_vec[1:d+1] .= param_vec[1]
    knot_vec[end-d:end] .= param_vec[end]

    # Compute cumulative parameter sum for internal knots using cumsum
    param_cumsum = cumsum(param_vec[2:end-1])

    if knotVecType == :Uniform
        # uniformly spaced knot vector
        # this method is not recommended because, if it is used with the chord length method for global interpolation,
        # the system of linear equations would be singular.
        for i in (d + 2):h
            knot_vec[i] = knot_vec[1] + (i - d - 1) // (h - d) * (knot_vec[end] - knot_vec[1])
        end
    elseif knotVecType == :Average
        # average spaced knot vector using improved distribution
        # The goal is to fill k[d+2:h] with h-d-1 values from param_vec[1] to param_vec[end] 
        # in a similar distribution to how the vector param_vec is distributed
        num_internal_knots = h - d - 1
        if num_internal_knots > 0
            # Use LinearInterpolation for better distribution across parameter domain
            param_interp = LinearInterpolation(param_vec, range(0, 1, length = n))
            if num_internal_knots == 1
                # Special case for single knot to avoid range issues
                knot_vec[d + 2] = param_interp(0.5)
            else
                internal_positions = range(0, 1, length = num_internal_knots)
                knot_vec[(d+2):h] .= param_interp.(internal_positions)
            end
        end
    end
    
    # control points
    control_points = zeros(eltype(u), size(u)[1:(end - 1)]..., h)
    control_points[array_axes..., 1] = u[array_axes..., 1]
    control_points[array_axes..., end] = u[array_axes..., end]
    data_residual = zeros(eltype(u), size(u)[1:(end - 1)]..., n)
    spline_coeffs = zeros(eltype(t), n, h)
    for i in 1:n
        spline_coefficients!(view(spline_coeffs, i, :), d, knot_vec, param_vec[i])
    end
    for k in 2:(n - 1)
        data_residual[array_axes...,
            k] = u[array_axes..., k] - spline_coeffs[k, 1] * u[array_axes..., 1] -
                 spline_coeffs[k, h] * u[array_axes..., end]
    end
    coeff_matrix = Array{T, N}(undef, size(u)[1:(end - 1)]..., h - 2)
    for i in 2:(h - 1)
        residual_sum = zeros(eltype(spline_coeffs), size(u)[1:(end - 1)]...)
        for k in 2:(n - 1)
            residual_sum = residual_sum + spline_coeffs[k, i] * data_residual[array_axes..., k]
        end
        coeff_matrix[array_axes..., i - 1] = residual_sum
    end
    spline_coeffs = spline_coeffs[2:(end - 1), 2:(h - 1)]
    system_matrix = transpose(spline_coeffs) * spline_coeffs
    coeff_matrix = reshape(coeff_matrix, prod(size(u)[1:(end - 1)]), :)
    solution = (system_matrix \ coeff_matrix')'
    solution = reshape(solution, size(u)[1:(end - 1)]..., :)
    control_points[array_axes..., 2:(end - 1)] = solution
    spline_coeffs = zeros(eltype(t), h)
    BSplineApprox(
        u, t, d, h, param_vec, knot_vec, control_points, spline_coeffs, pVecType, knotVecType,
        extrapolation_left, extrapolation_right, assume_linear_t)
end
"""
    CubicHermiteSpline(du, u, t; extrapolation::ExtrapolationType.T = ExtrapolationType.None, extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None, cache_parameters = false)

It is a Cubic Hermite interpolation, which is a piece-wise third degree polynomial such that the value and the first derivative are equal to given values in the data points.

## Arguments

  - `du`: the derivative at the data points.
  - `u`: data points.
  - `t`: time points.

## Keyword Arguments

  - `extrapolation`: The extrapolation type applied left and right of the data. Possible options
    are `ExtrapolationType.None` (default), `ExtrapolationType.Constant`, `ExtrapolationType.Linear`
    `ExtrapolationType.Extension`, `ExtrapolationType.Periodic` and `ExtrapolationType.Reflective`.
  - `extrapolation_left`: The extrapolation type applied left of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `extrapolation_right`: The extrapolation type applied right of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `cache_parameters`: precompute parameters at initialization for faster interpolation computations. Note: if activated, `u` and `t` should not be modified. Defaults to `false`.
  - `assume_linear_t`: boolean value to specify a faster index lookup behaviour for
    evenly-distributed abscissae. Alternatively, a numerical threshold may be specified
    for a test based on the normalized standard deviation of the difference with respect
    to the straight line (see [`looks_linear`](@ref)). Defaults to 1e-2.
"""
struct CubicHermiteSpline{uType, tType, IType, duType, pType, T} <:
       AbstractInterpolation{T}
    du::duType
    u::uType
    t::tType
    I::IType
    p::CubicHermiteParameterCache{pType}
    extrapolation_left::ExtrapolationType.T
    extrapolation_right::ExtrapolationType.T
    iguesser::Guesser{tType}
    cache_parameters::Bool
    linear_lookup::Bool
    function CubicHermiteSpline(
            du, u, t, I, p, extrapolation_left, extrapolation_right,
            cache_parameters, assume_linear_t)
        linear_lookup = seems_linear(assume_linear_t, t)
        new{typeof(u), typeof(t), typeof(I), typeof(du), typeof(p.c₁), eltype(u)}(
            du, u, t, I, p, extrapolation_left, extrapolation_right,
            Guesser(t), cache_parameters, linear_lookup)
    end
end

function CubicHermiteSpline(
        du, u, t; extrapolation::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None, cache_parameters = false, assume_linear_t = 1e-2)
    @assert length(u)==length(du) "Length of `u` is not equal to length of `du`."
    extrapolation_left,
    extrapolation_right = munge_extrapolation(
        extrapolation, extrapolation_left, extrapolation_right)
    u, t = munge_data(u, t)
    linear_lookup = seems_linear(assume_linear_t, t)
    p = CubicHermiteParameterCache(du, u, t, cache_parameters)
    A = CubicHermiteSpline(
        du, u, t, nothing, p, extrapolation_left,
        extrapolation_right, cache_parameters, linear_lookup)
    I = cumulative_integral(A, cache_parameters)
    CubicHermiteSpline(du, u, t, I, p, extrapolation_left,
        extrapolation_right, cache_parameters, linear_lookup)
end

"""
    PCHIPInterpolation(u, t; extrapolation::ExtrapolationType.T = ExtrapolationType.None, extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None)

It is a PCHIP Interpolation, which is a type of [`CubicHermiteSpline`](@ref) where the derivative values `du` are derived from the input data
in such a way that the interpolation never overshoots the data. See [here](https://www.mathworks.com/content/dam/mathworks/mathworks-dot-com/moler/interp.pdf),
section 3.4 for more details.

## Arguments

  - `u`: data points.
  - `t`: time points.

## Keyword Arguments

  - `extrapolation`: The extrapolation type applied left and right of the data. Possible options
    are `ExtrapolationType.None` (default), `ExtrapolationType.Constant`, `ExtrapolationType.Linear`
    `ExtrapolationType.Extension`, `ExtrapolationType.Periodic` and `ExtrapolationType.Reflective`.
  - `extrapolation_left`: The extrapolation type applied left of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `extrapolation_right`: The extrapolation type applied right of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `cache_parameters`: precompute parameters at initialization for faster interpolation computations. Note: if activated, `u` and `t` should not be modified. Defaults to `false`.
  - `assume_linear_t`: boolean value to specify a faster index lookup behaviour for
    evenly-distributed abscissae. Alternatively, a numerical threshold may be specified
    for a test based on the normalized standard deviation of the difference with respect
    to the straight line (see [`looks_linear`](@ref)). Defaults to 1e-2.
"""
function PCHIPInterpolation(u, t; kwargs...)
    u, t = munge_data(u, t)
    du = du_PCHIP(u, t)
    CubicHermiteSpline(du, u, t; kwargs...)
end

"""
    QuinticHermiteSpline(ddu, du, u, t; extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None)

It is a Quintic Hermite interpolation, which is a piece-wise fifth degree polynomial such that the value and the first and second derivative are equal to given values in the data points.

## Arguments

  - `ddu`: the second derivative at the data points.
  - `du`: the derivative at the data points.
  - `u`: data points.
  - `t`: time points.

## Keyword Arguments

  - `extrapolation`: The extrapolation type applied left and right of the data. Possible options
    are `ExtrapolationType.None` (default), `ExtrapolationType.Constant`, `ExtrapolationType.Linear`
    `ExtrapolationType.Extension`, `ExtrapolationType.Periodic` and `ExtrapolationType.Reflective`.
  - `extrapolation_left`: The extrapolation type applied left of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `extrapolation_right`: The extrapolation type applied right of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `cache_parameters`: precompute parameters at initialization for faster interpolation computations. Note: if activated, `u` and `t` should not be modified. Defaults to `false`.
  - `assume_linear_t`: boolean value to specify a faster index lookup behaviour for
    evenly-distributed abscissae. Alternatively, a numerical threshold may be specified
    for a test based on the normalized standard deviation of the difference with respect
    to the straight line (see [`looks_linear`](@ref)). Defaults to 1e-2.
"""
struct QuinticHermiteSpline{uType, tType, IType, duType, dduType, pType, T} <:
       AbstractInterpolation{T}
    ddu::dduType
    du::duType
    u::uType
    t::tType
    I::IType
    p::QuinticHermiteParameterCache{pType}
    extrapolation_left::ExtrapolationType.T
    extrapolation_right::ExtrapolationType.T
    iguesser::Guesser{tType}
    cache_parameters::Bool
    linear_lookup::Bool
    function QuinticHermiteSpline(
            ddu, du, u, t, I, p, extrapolation_left,
            extrapolation_right, cache_parameters, assume_linear_t)
        linear_lookup = seems_linear(assume_linear_t, t)
        new{typeof(u), typeof(t), typeof(I), typeof(du),
            typeof(ddu), typeof(p.c₁), eltype(u)}(
            ddu, du, u, t, I, p, extrapolation_left, extrapolation_right,
            Guesser(t), cache_parameters, linear_lookup)
    end
end

function QuinticHermiteSpline(
        ddu, du, u, t; extrapolation::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None,
        cache_parameters = false, assume_linear_t = 1e-2)
    @assert length(u)==length(du)==length(ddu) "Length of `u` is not equal to length of `du` or `ddu`."
    extrapolation_left,
    extrapolation_right = munge_extrapolation(
        extrapolation, extrapolation_left, extrapolation_right)
    u, t = munge_data(u, t)
    linear_lookup = seems_linear(assume_linear_t, t)
    p = QuinticHermiteParameterCache(ddu, du, u, t, cache_parameters)
    A = QuinticHermiteSpline(
        ddu, du, u, t, nothing, p, extrapolation_left,
        extrapolation_right, cache_parameters, linear_lookup)
    I = cumulative_integral(A, cache_parameters)
    QuinticHermiteSpline(
        ddu, du, u, t, I, p, extrapolation_left,
        extrapolation_right, cache_parameters, linear_lookup)
end

struct SmoothArcLengthInterpolation{
    uType, tType, IType, P, D, S <: Union{AbstractInterpolation, Nothing}, T} <:
       AbstractInterpolation{T}
    u::uType
    t::tType
    d::Matrix{D}
    shape_itp::S
    Δt_circle_segment::Vector{P}
    Δt_line_segment::Vector{P}
    center::Matrix{P}
    radius::Vector{P}
    dir_1::Matrix{P}
    dir_2::Matrix{P}
    # short_side_left[i] = true means that the line segment comes after the circle segment
    short_side_left::Vector{Bool}
    I::IType
    p::Nothing
    extrapolation_left::ExtrapolationType.T
    extrapolation_right::ExtrapolationType.T
    iguesser::Guesser{tType}
    cache_parameters::Bool
    linear_lookup::Bool
    out::Vector{P}
    derivative::Vector{P}
    in_place::Bool
    function SmoothArcLengthInterpolation(
            u, t, d, shape_itp, Δt_circle_segment, Δt_line_segment,
            center, radius, dir_1, dir_2, short_side_left,
            I, extrapolation_left, extrapolation_right,
            assume_linear_t, out, derivative, in_place)
        linear_lookup = seems_linear(assume_linear_t, t)
        new{typeof(u), typeof(t), typeof(I), eltype(radius),
            eltype(d), typeof(shape_itp), eltype(u)}(
            u, t, d, shape_itp, Δt_circle_segment, Δt_line_segment,
            center, radius, dir_1, dir_2, short_side_left,
            I, nothing, extrapolation_left, extrapolation_right,
            Guesser(t), false, linear_lookup, out, derivative, in_place
        )
    end
end

"""
     SmoothArcLengthInterpolation(
        u::AbstractMatrix{U};
        t::Union{AbstractVector, Nothing} = nothing,
        interpolation_type::Type{<:AbstractInterpolation} = QuadraticSpline,
        kwargs...) where {U}

Interpolate in a C¹ smooth way trough the data with unit speed by approximating
an interpolation (the shape interpolation) with line segments and circle segments.

## Arguments

  - `u`: The data to be interpolated in matrix form; (ndim, ndata).

NOTE: With this method it is not possible to pass keyword arguments to the constructor of the shape interpolation.
If you want to do this, construct the shape interpolation yourself and use the
`SmoothArcLengthInterpolation(shape_itp::AbstractInterpolation; kwargs...)` method.

## Keyword Arguments

  - `t`: The time points of the shape interpolation. By default given by the cumulative sum of the Euclidean
    distances between the points `u`.
  - `interpolation_type`: The type of the shape interpolation. Defaults to `QuadraticSpline`. Note that
    for the `SmoothArcLengthInterpolation` to be C¹ smooth, the `interpolation_type` must be C¹ smooth as well.
  - `m`: The number of points at which the shape interpolation is evaluated in each interval between time points.
    The `SmoothArcLengthInterpolation` converges to the shape interpolation (in shape) as m → ∞.
  - `extrapolation`: The extrapolation type applied left and right of the data. Possible options
    are `ExtrapolationType.None` (default), `ExtrapolationType.Constant`, `ExtrapolationType.Linear`
    `ExtrapolationType.Extension`, `ExtrapolationType.Periodic` and `ExtrapolationType.Reflective`.
  - `extrapolation_left`: The extrapolation type applied left of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `extrapolation_right`: The extrapolation type applied right of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `assume_linear_t`: boolean value to specify a faster index lookup behaviour for
    evenly-distributed abscissae. Alternatively, a numerical threshold may be specified
    for a test based on the normalized standard deviation of the difference with respect
    to the straight line (see [`looks_linear`](@ref)). Defaults to 1e-2.
"""
function SmoothArcLengthInterpolation(
        u::AbstractMatrix{U};
        t::Union{AbstractVector, Nothing} = nothing,
        interpolation_type::Type{<:AbstractInterpolation} = QuadraticSpline,
        kwargs...
) where {U}
    if isnothing(t)
        # Compute default t based on point distances
        N, n = size(u)
        t = Vector{U}(undef, n)
        t[1] = zero(U)
        Δu = Vector{U}(undef, N)
        for i in 2:n
            @. Δu = u[:, i] - u[:, i - 1]
            t[i] = t[i - 1] + norm(Δu)
        end
    end
    shape_itp = interpolation_type(collect.(eachcol(u)), t)
    SmoothArcLengthInterpolation(shape_itp; kwargs...)
end

"""
    function SmoothArcLengthInterpolation(
            shape_itp::AbstractInterpolation;
            m::Integer = 2,
            kwargs...)

Approximate the `shape_itp` with a C¹ unit speed interpolation using line segments and circle segments.

## Arguments

  - `shape_itp`: The interpolation to be approximated. Note that
    for the `SmoothArcLengthInterpolation` to be C¹ smooth, the `shape_itp` must be C¹ smooth as well.

## Keyword Arguments

  - `m`: The number of points at which the shape interpolation is evaluated in each interval between time points.
    The `SmoothArcLengthInterpolation` converges to the shape interpolation (in shape) as m → ∞.
  - `extrapolation`: The extrapolation type applied left and right of the data. Possible options
    are `ExtrapolationType.None` (default), `ExtrapolationType.Constant`, `ExtrapolationType.Linear`
    `ExtrapolationType.Extension`, `ExtrapolationType.Periodic` and `ExtrapolationType.Reflective`.
  - `extrapolation_left`: The extrapolation type applied left of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `extrapolation_right`: The extrapolation type applied right of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `assume_linear_t`: boolean value to specify a faster index lookup behaviour for
    evenly-distributed abscissae. Alternatively, a numerical threshold may be specified
    for a test based on the normalized standard deviation of the difference with respect
    to the straight line (see [`looks_linear`](@ref)). Defaults to 1e-2.
"""
function SmoothArcLengthInterpolation(
        shape_itp::AbstractInterpolation;
        m::Integer = 2,
        kwargs...
)
    (; u, t) = shape_itp
    T = promote_type(eltype(eltype(u)), eltype(t))

    # Resp. the output dimensionality and the number of data points in the original interpolation
    N = length(first(u))
    n = length(u)

    # Number of points defining the tangent curve of shape_itp
    n_tilde = m * (n - 1) + 1

    # The evaluations of shape_itp
    u_tilde = Matrix{T}(undef, N, n_tilde)
    d_tilde = Matrix{T}(undef, N, n_tilde)

    j = 1

    for i in 1:(n - 1)
        for t_eval in range(t[i], t[i + 1], length = m + 1)[1:(end - 1)]
            u_tilde[:, j] .= shape_itp(t_eval)
            d_tilde[:, j] .= derivative(shape_itp, t_eval)
            normalize!(view(d_tilde, :, j))
            j += 1
        end
    end

    u_tilde[:, end] .= shape_itp(last(t))
    d_tilde[:, end] .= derivative(shape_itp, last(t))
    normalize!(view(d_tilde, :, n_tilde))

    return SmoothArcLengthInterpolation(u_tilde, d_tilde; shape_itp, kwargs...)
end

"""
    function SmoothArcLengthInterpolation(
        u::AbstractMatrix,
        d::AbstractMatrix
        [, make_intersections::Val{<:Bool}];
        shape_itp::Union{AbstractInterpolation, Nothing} = nothing,
        extrapolation::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None,
        cache_parameters::Bool = false,
        assume_linear_t = 1e-2,
        in_place::Bool = true)

Make a C¹ smooth unit speed interpolation through the given data with the given tangents using line
segments and circle segments.

## Arguments

  - `u`: The data to be interpolated in matrix form; (ndim, ndata).
  - `d`: The tangents to the curve in the points `u`.
  - `make_intersections`: Whether additional (point, tangent) pairs have to be added in between the provided
    data to ensure that the consecutive (tangent) lines intersect. Defaults to `Val(true)`.

## Keyword Arguments

  - `shape_itp`: The interpolation that is being approximated, if one exists. Note that this
    interpolation is not being used; it is just passed along to keep track of where the shape
    of the `SmoothArcLengthInterpolation` originated.
  - `extrapolation`: The extrapolation type applied left and right of the data. Possible options
    are `ExtrapolationType.None` (default), `ExtrapolationType.Constant`, `ExtrapolationType.Linear`
    `ExtrapolationType.Extension`, `ExtrapolationType.Periodic` and `ExtrapolationType.Reflective`.
  - `extrapolation_left`: The extrapolation type applied left of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `extrapolation_right`: The extrapolation type applied right of the data. See `extrapolation` for
    the possible options. This keyword is ignored if `extrapolation != Extrapolation.none`.
  - `assume_linear_t`: boolean value to specify a faster index lookup behaviour for
    evenly-distributed abscissae. Alternatively, a numerical threshold may be specified
    for a test based on the normalized standard deviation of the difference with respect
    to the straight line (see [`looks_linear`](@ref)). Defaults to 1e-2.
"""
function SmoothArcLengthInterpolation(
        u::AbstractMatrix,
        d::AbstractMatrix;
        kwargs...)
    SmoothArcLengthInterpolation(u, d, Val{true}(); kwargs...)
end

function SmoothArcLengthInterpolation(
        u::AbstractMatrix,
        d::AbstractMatrix,
        make_intersections::Val{true};
        kwargs...)
    N, n = size(u)

    # Number of points in the augmented tangent curve
    n_hat = 2 * n - 1

    # The data defining the augmented tangent curve
    T = promote_type(eltype(eltype(u)), eltype(d))
    u_hat = Matrix{T}(undef, N, n_hat)
    d_hat = Matrix{T}(undef, N, n_hat)

    k = 1
    Δu = Vector{T}(undef, N)
    uⱼ_close_left = Vector{T}(undef, N)
    uⱼ_close_right = Vector{T}(undef, N)
    uⱼ_int = Vector{T}(undef, N)
    uⱼ_int_left = Vector{T}(undef, N)
    uⱼ_int_right = Vector{T}(undef, N)

    for j in 1:(n - 1)
        u_hat[:, k] .= u[:, j]
        d_hat[:, k] .= d[:, j]

        uⱼ, uⱼ₊₁, dⱼ, dⱼ₊₁, d_inner = smooth_arc_length_params_1!(Δu, u, d, j)

        inner_1 = dot(Δu, dⱼ)
        inner_2 = dot(Δu, dⱼ₊₁)
        denom = 1 - d_inner^2
        dⱼ_coef = (inner_1 - d_inner * inner_2) / denom
        dⱼ₊₁_coef = (d_inner * inner_1 - inner_2) / denom

        if !((dⱼ_coef >= 0) && (dⱼ₊₁_coef <= 0))
            error("Some consecutive tangent lines do not converge, consider increasing m.")
        end

        @. uⱼ_close_left = uⱼ + dⱼ_coef * dⱼ
        @. uⱼ_close_right = uⱼ₊₁ +
                            dⱼ₊₁_coef * dⱼ₊₁
        @. uⱼ_int = (uⱼ_close_left + uⱼ_close_right) / 2

        # compute δ_star
        δⱼ, _, _ = smooth_arc_length_params_2(uⱼ_int, uⱼ, uⱼ₊₁)
        δⱼ_star = δⱼ * (2 - sqrt(2 + 2 * d_inner)) / (1 - d_inner)

        # Compute the points whose connecting line defines the tangent curve augmenting point
        @. uⱼ_int_left = uⱼ_close_left - δⱼ_star * dⱼ
        @. uⱼ_int_right = uⱼ_close_right + δⱼ_star * dⱼ₊₁

        # Compute tangent curve augmenting point
        uⱼ_plus_half = view(u_hat, :, k + 1)
        dⱼ_plus_half = view(d_hat, :, k + 1)

        @. uⱼ_plus_half = (uⱼ_int_left + uⱼ_int_right) / 2
        @. dⱼ_plus_half = uⱼ_int_right - uⱼ_int_left
        normalize!(dⱼ_plus_half)

        k += 2
    end

    u_hat[:, end] .= u[:, end]
    d_hat[:, end] .= d[:, end]

    return SmoothArcLengthInterpolation(u_hat, d_hat, Val{false}(); kwargs...)
end

function SmoothArcLengthInterpolation(
        u::AbstractMatrix,
        d::AbstractMatrix,
        ::Val{false};
        shape_itp::Union{AbstractInterpolation, Nothing} = nothing,
        extrapolation::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None,
        cache_parameters::Bool = false,
        assume_linear_t = 1e-2,
        in_place::Bool = true)
    N = size(u, 1)
    n_circle_arcs = size(u, 2) - 1

    P = promote_type(eltype(u), eltype(d))
    t = zeros(P, n_circle_arcs + 1)
    Δt_circle_segment = zeros(P, n_circle_arcs)
    Δt_line_segment = zeros(P, n_circle_arcs)
    center = Matrix{P}(undef, N, n_circle_arcs)
    radius = Vector{P}(undef, n_circle_arcs)
    dir_1 = Matrix{P}(undef, N, n_circle_arcs)
    dir_2 = Matrix{P}(undef, N, n_circle_arcs)
    short_side_left = zeros(Bool, n_circle_arcs)

    # Intermediate results
    Δu = Vector{P}(undef, N)
    u_int = Vector{P}(undef, N)

    # Compute circle segments and line segments
    for j in 1:n_circle_arcs
        uⱼ, uⱼ₊₁, dⱼ, dⱼ₊₁, d_inner = smooth_arc_length_params_1!(Δu, u, d, j)

        dⱼ_coef = (dot(Δu, dⱼ) - d_inner * dot(Δu, dⱼ₊₁)) / (1 - d_inner^2)
        @. u_int = uⱼ + dⱼ_coef * dⱼ

        δⱼ, short_side_left_, Δt_line_seg = smooth_arc_length_params_2(u_int, uⱼ, uⱼ₊₁)
        short_side_left[j] = short_side_left_

        Rⱼ = δⱼ * sqrt((1 + d_inner) / (1 - d_inner))
        radius[j] = Rⱼ
        cⱼ = view(center, :, j)
        v₁ = view(dir_1, :, j)
        v₂ = view(dir_2, :, j)

        @. cⱼ = u_int + δⱼ * (dⱼ₊₁ - dⱼ) / (1 - d_inner)
        @. v₁ = -δⱼ * (dⱼ₊₁ - d_inner * dⱼ) / (1 - d_inner)
        @. v₂ = Rⱼ * dⱼ

        Δt_circle_seg = 2Rⱼ * atan(δⱼ, Rⱼ)
        Δt_circle_segment[j] = Δt_circle_seg
        Δt_line_segment[j] = Δt_line_seg

        t[j + 1] = t[j] + Δt_circle_seg + Δt_line_seg
    end

    extrapolation_left,
    extrapolation_right = munge_extrapolation(
        extrapolation, extrapolation_left, extrapolation_right)
    linear_lookup = seems_linear(assume_linear_t, t)

    out = Vector{P}(undef, N)
    derivative = Vector{P}(undef, N)

    return SmoothArcLengthInterpolation(
        u, t, d, shape_itp, Δt_circle_segment, Δt_line_segment,
        center, radius, dir_1, dir_2, short_side_left,
        nothing, extrapolation_left, extrapolation_right, linear_lookup, out, derivative, in_place)
end
