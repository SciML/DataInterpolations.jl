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
    extrapolation_left, extrapolation_right = munge_extrapolation(
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
    extrapolation_left, extrapolation_right = munge_extrapolation(
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
    extrapolation_left, extrapolation_right = munge_extrapolation(
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
    extrapolation_left, extrapolation_right = munge_extrapolation(
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
    extrapolation_left, extrapolation_right = munge_extrapolation(
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
    extrapolation_left, extrapolation_right = munge_extrapolation(
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
    extrapolation_left, extrapolation_right = munge_extrapolation(
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
    extrapolation_left, extrapolation_right = munge_extrapolation(
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
    extrapolation_left, extrapolation_right = munge_extrapolation(
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
    extrapolation_left, extrapolation_right = munge_extrapolation(
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
    extrapolation_left, extrapolation_right = munge_extrapolation(
        extrapolation, extrapolation_left, extrapolation_right)
    u, t = munge_data(u, t)
    n = length(t)
    n < d + 1 && error("BSplineInterpolation needs at least d + 1, i.e. $(d+1) points.")
    s = zero(eltype(u))
    p = zero(t)
    k = zeros(eltype(t), n + d + 1)
    l = zeros(eltype(u), n - 1)
    p[1] = zero(eltype(t))
    p[end] = one(eltype(t))

    for i in 2:n
        s += √((t[i] - t[i - 1])^2 + (u[i] - u[i - 1])^2)
        l[i - 1] = s
    end
    if pVecType == :Uniform
        for i in 2:(n - 1)
            p[i] = p[1] + (i - 1) * (p[end] - p[1]) / (n - 1)
        end
    elseif pVecType == :ArcLen
        for i in 2:(n - 1)
            p[i] = p[1] + l[i - 1] / s * (p[end] - p[1])
        end
    end

    lidx = 1
    ridx = length(k)
    while lidx <= (d + 1) && ridx >= (length(k) - d)
        k[lidx] = p[1]
        k[ridx] = p[end]
        lidx += 1
        ridx -= 1
    end

    ps = zeros(eltype(t), n - 2)
    s = zero(eltype(t))
    for i in 2:(n - 1)
        s += p[i]
        ps[i - 1] = s
    end

    if knotVecType == :Uniform
        # uniformly spaced knot vector
        # this method is not recommended because, if it is used with the chord length method for global interpolation,
        # the system of linear equations would be singular.
        for i in (d + 2):n
            k[i] = k[1] + (i - d - 1) // (n - d) * (k[end] - k[1])
        end
    elseif knotVecType == :Average
        # average spaced knot vector
        idx = 1
        if d + 2 <= n
            k[d + 2] = 1 // d * ps[d]
        end
        for i in (d + 3):n
            k[i] = 1 // d * (ps[idx + d] - ps[idx])
            idx += 1
        end
    end
    # control points
    sc = zeros(eltype(t), n, n)
    spline_coefficients!(sc, d, k, p)
    c = vec(sc \ u[:, :])
    sc = zeros(eltype(t), n)
    BSplineInterpolation(
        u, t, d, p, k, c, sc, pVecType, knotVecType,
        extrapolation_left, extrapolation_right, assume_linear_t)
end

function BSplineInterpolation(
        u::AbstractArray, t, d, pVecType, knotVecType;
        extrapolation::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None,
        assume_linear_t = 1e-2)
    extrapolation_left, extrapolation_right = munge_extrapolation(
        extrapolation, extrapolation_left, extrapolation_right)
    u, t = munge_data(u, t)
    n = length(t)
    n < d + 1 && error("BSplineInterpolation needs at least d + 1, i.e. $(d+1) points.")
    s = zero(eltype(u))
    p = zero(t)
    k = zeros(eltype(t), n + d + 1)
    l = zeros(eltype(u), n - 1)
    p[1] = zero(eltype(t))
    p[end] = one(eltype(t))

    ax_u = axes(u)[1:(end - 1)]

    for i in 2:n
        s += √((t[i] - t[i - 1])^2 + sum((u[ax_u..., i] - u[ax_u..., i - 1]) .^ 2))
        l[i - 1] = s
    end
    if pVecType == :Uniform
        for i in 2:(n - 1)
            p[i] = p[1] + (i - 1) * (p[end] - p[1]) / (n - 1)
        end
    elseif pVecType == :ArcLen
        for i in 2:(n - 1)
            p[i] = p[1] + l[i - 1] / s * (p[end] - p[1])
        end
    end

    lidx = 1
    ridx = length(k)
    while lidx <= (d + 1) && ridx >= (length(k) - d)
        k[lidx] = p[1]
        k[ridx] = p[end]
        lidx += 1
        ridx -= 1
    end

    ps = zeros(eltype(t), n - 2)
    s = zero(eltype(t))
    for i in 2:(n - 1)
        s += p[i]
        ps[i - 1] = s
    end

    if knotVecType == :Uniform
        # uniformly spaced knot vector
        # this method is not recommended because, if it is used with the chord length method for global interpolation,
        # the system of linear equations would be singular.
        for i in (d + 2):n
            k[i] = k[1] + (i - d - 1) // (n - d) * (k[end] - k[1])
        end
    elseif knotVecType == :Average
        # average spaced knot vector
        idx = 1
        if d + 2 <= n
            k[d + 2] = 1 // d * ps[d]
        end
        for i in (d + 3):n
            k[i] = 1 // d * (ps[idx + d] - ps[idx])
            idx += 1
        end
    end
    # control points
    sc = zeros(eltype(t), n, n)
    spline_coefficients!(sc, d, k, p)
    c = (sc \ reshape(u, prod(size(u)[1:(end - 1)]), :)')'
    c = reshape(c, size(u)...)
    sc = zeros(eltype(t), n)
    BSplineInterpolation(
        u, t, d, p, k, c, sc, pVecType, knotVecType,
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
    extrapolation_left, extrapolation_right = munge_extrapolation(
        extrapolation, extrapolation_left, extrapolation_right)
    u, t = munge_data(u, t)
    n = length(t)
    h < d + 1 && error("BSplineApprox needs at least d + 1, i.e. $(d+1) control points.")
    s = zero(eltype(u))
    p = zero(t)
    k = zeros(eltype(t), h + d + 1)
    l = zeros(eltype(u), n - 1)
    p[1] = zero(eltype(t))
    p[end] = one(eltype(t))

    for i in 2:n
        s += √((t[i] - t[i - 1])^2 + (u[i] - u[i - 1])^2)
        l[i - 1] = s
    end
    if pVecType == :Uniform
        for i in 2:(n - 1)
            p[i] = p[1] + (i - 1) * (p[end] - p[1]) / (n - 1)
        end
    elseif pVecType == :ArcLen
        for i in 2:(n - 1)
            p[i] = p[1] + l[i - 1] / s * (p[end] - p[1])
        end
    end

    lidx = 1
    ridx = length(k)
    while lidx <= (d + 1) && ridx >= (length(k) - d)
        k[lidx] = p[1]
        k[ridx] = p[end]
        lidx += 1
        ridx -= 1
    end

    ps = zeros(eltype(t), n - 2)
    s = zero(eltype(t))
    for i in 2:(n - 1)
        s += p[i]
        ps[i - 1] = s
    end

    if knotVecType == :Uniform
        # uniformly spaced knot vector
        # this method is not recommended because, if it is used with the chord length method for global interpolation,
        # the system of linear equations would be singular.
        for i in (d + 2):h
            k[i] = k[1] + (i - d - 1) // (h - d) * (k[end] - k[1])
        end
    elseif knotVecType == :Average
        # NOTE: verify that average method can be applied when size of k is less than size of p
        # average spaced knot vector
        idx = 1
        if d + 2 <= h
            k[d + 2] = 1 // d * ps[d]
        end
        for i in (d + 3):h
            k[i] = 1 // d * (ps[idx + d] - ps[idx])
            idx += 1
        end
    end
    # control points
    c = zeros(eltype(u), h)
    c[1] = u[1]
    c[end] = u[end]
    q = zeros(eltype(u), n)
    sc = zeros(eltype(t), n, h)
    for i in 1:n
        spline_coefficients!(view(sc, i, :), d, k, p[i])
    end
    for k in 2:(n - 1)
        q[k] = u[k] - sc[k, 1] * u[1] - sc[k, h] * u[end]
    end
    Q = Matrix{eltype(u)}(undef, h - 2, 1)
    for i in 2:(h - 1)
        s = 0.0
        for k in 2:(n - 1)
            s += sc[k, i] * q[k]
        end
        Q[i - 1] = s
    end
    sc = sc[2:(end - 1), 2:(h - 1)]
    M = transpose(sc) * sc
    P = M \ Q
    c[2:(end - 1)] .= vec(P)
    sc = zeros(eltype(t), h)
    BSplineApprox(
        u, t, d, h, p, k, c, sc, pVecType, knotVecType,
        extrapolation_left, extrapolation_right, assume_linear_t)
end

function BSplineApprox(
        u::AbstractArray{T, N}, t, d, h, pVecType, knotVecType;
        extrapolation::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None,
        assume_linear_t = 1e-2) where {T, N}
    extrapolation_left, extrapolation_right = munge_extrapolation(
        extrapolation, extrapolation_left, extrapolation_right)
    u, t = munge_data(u, t)
    n = length(t)
    h < d + 1 && error("BSplineApprox needs at least d + 1, i.e. $(d+1) control points.")
    s = zero(eltype(u))
    p = zero(t)
    k = zeros(eltype(t), h + d + 1)
    l = zeros(eltype(u), n - 1)
    p[1] = zero(eltype(t))
    p[end] = one(eltype(t))

    ax_u = axes(u)[1:(end - 1)]

    for i in 2:n
        s += √((t[i] - t[i - 1])^2 + sum((u[ax_u..., i] - u[ax_u..., i - 1]) .^ 2))
        l[i - 1] = s
    end
    if pVecType == :Uniform
        for i in 2:(n - 1)
            p[i] = p[1] + (i - 1) * (p[end] - p[1]) / (n - 1)
        end
    elseif pVecType == :ArcLen
        for i in 2:(n - 1)
            p[i] = p[1] + l[i - 1] / s * (p[end] - p[1])
        end
    end

    lidx = 1
    ridx = length(k)
    while lidx <= (d + 1) && ridx >= (length(k) - d)
        k[lidx] = p[1]
        k[ridx] = p[end]
        lidx += 1
        ridx -= 1
    end

    ps = zeros(eltype(t), n - 2)
    s = zero(eltype(t))
    for i in 2:(n - 1)
        s += p[i]
        ps[i - 1] = s
    end

    if knotVecType == :Uniform
        # uniformly spaced knot vector
        # this method is not recommended because, if it is used with the chord length method for global interpolation,
        # the system of linear equations would be singular.
        for i in (d + 2):h
            k[i] = k[1] + (i - d - 1) // (h - d) * (k[end] - k[1])
        end
    elseif knotVecType == :Average
        # NOTE: verify that average method can be applied when size of k is less than size of p
        # average spaced knot vector
        idx = 1
        if d + 2 <= h
            k[d + 2] = 1 // d * ps[d]
        end
        for i in (d + 3):h
            k[i] = 1 // d * (ps[idx + d] - ps[idx])
            idx += 1
        end
    end
    # control points
    c = zeros(eltype(u), size(u)[1:(end - 1)]..., h)
    c[ax_u..., 1] = u[ax_u..., 1]
    c[ax_u..., end] = u[ax_u..., end]
    q = zeros(eltype(u), size(u)[1:(end - 1)]..., n)
    sc = zeros(eltype(t), n, h)
    for i in 1:n
        spline_coefficients!(view(sc, i, :), d, k, p[i])
    end
    for k in 2:(n - 1)
        q[ax_u..., k] = u[ax_u..., k] - sc[k, 1] * u[ax_u..., 1] -
                        sc[k, h] * u[ax_u..., end]
    end
    Q = Array{T, N}(undef, size(u)[1:(end - 1)]..., h - 2)
    for i in 2:(h - 1)
        s = zeros(eltype(sc), size(u)[1:(end - 1)]...)
        for k in 2:(n - 1)
            s = s + sc[k, i] * q[ax_u..., k]
        end
        Q[ax_u..., i - 1] = s
    end
    sc = sc[2:(end - 1), 2:(h - 1)]
    M = transpose(sc) * sc
    Q = reshape(Q, prod(size(u)[1:(end - 1)]), :)
    P = (M \ Q')'
    P = reshape(P, size(u)[1:(end - 1)]..., :)
    c[ax_u..., 2:(end - 1)] = P
    sc = zeros(eltype(t), h)
    BSplineApprox(
        u, t, d, h, p, k, c, sc, pVecType, knotVecType,
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
    extrapolation_left, extrapolation_right = munge_extrapolation(
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
    extrapolation_left, extrapolation_right = munge_extrapolation(
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

"""
    Signature...

Interpolate in a C¹ smooth way trough the data with unit speed.

## Arguments

...

## Keyword Arguments

...
"""
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

function SmoothArcLengthInterpolation(
        shape_itp::AbstractInterpolation; m::Integer = 2,
        kwargs...
)
    (; u, t) = shape_itp
    T = promote_type(eltype(eltype(u)), eltype(t))

    # Resp. the output dimensionality and the number of data points in the original interpolation
    N = length(first(u))
    n = length(u)

    # Number of points defining the tangent curve of itp
    n_tilde = m * (n - 1) + 1

    # The evaluations of itp
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

    # Number of points in the augmented tangent curve of itp
    n_hat = 2 * n_tilde - 1

    # The data defining the augmented tangent curve of itp
    u_hat = Matrix{T}(undef, N, n_hat)
    d_hat = Matrix{T}(undef, N, n_hat)

    k = 1
    Δu_tilde = Vector{T}(undef, N)
    u_tilde_j_close_left = Vector{T}(undef, N)
    u_tilde_j_close_right = Vector{T}(undef, N)
    u_tilde_j_int = Vector{T}(undef, N)
    u_tilde_j_int_left = Vector{T}(undef, N)
    u_tilde_j_int_right = Vector{T}(undef, N)

    for j in 1:(n_tilde - 1)
        u_hat[:, k] .= u_tilde[:, j]
        d_hat[:, k] .= d_tilde[:, j]

        u_tilde_j = view(u_tilde, :, j)
        d_tilde_j = view(d_tilde, :, j)
        u_tilde_j_plus_1 = view(u_tilde, :, j + 1)
        d_tilde_j_plus_1 = view(d_tilde, :, j + 1)
        d_tilde_inner = dot(d_tilde_j, d_tilde_j_plus_1)

        @. Δu_tilde = u_tilde_j_plus_1 - u_tilde_j

        inner_1 = dot(Δu_tilde, d_tilde_j)
        inner_2 = dot(Δu_tilde, d_tilde_j_plus_1)
        denom = 1 - d_tilde_inner^2
        d_tilde_j_coef = (inner_1 - d_tilde_inner * inner_2) / denom
        d_tilde_j_plus_1_coef = (d_tilde_inner * inner_1 - inner_2) / denom

        if !((d_tilde_j_coef >= 0) && (d_tilde_j_plus_1_coef <= 0))
            error("Some consecutive tangent lines do not converge, consider increasing m.")
        end

        @. u_tilde_j_close_left = u_tilde_j + d_tilde_j_coef * d_tilde_j
        @. u_tilde_j_close_right = u_tilde_j_plus_1 +
                                   d_tilde_j_plus_1_coef * d_tilde_j_plus_1
        @. u_tilde_j_int = (u_tilde_j_close_left + u_tilde_j_close_right) / 2

        # compute δ_star
        δ_j = min(
            euclidean(u_tilde_j_int, u_tilde_j),
            euclidean(u_tilde_j_int, u_tilde_j_plus_1)
        )
        δ_j_star = δ_j * (2 - sqrt(2 + 2 * d_tilde_inner)) / (1 - d_tilde_inner)

        # Compute the points whose connecting line defines the tangent curve augmenting point
        @. u_tilde_j_int_left = u_tilde_j_close_left - δ_j_star * d_tilde_j
        @. u_tilde_j_int_right = u_tilde_j_close_right + δ_j_star * d_tilde_j_plus_1

        # Compute tangent curve augmenting point
        u_tilde_j_plus_half = view(u_hat, :, k + 1)
        d_tilde_j_plus_half = view(d_hat, :, k + 1)

        @. u_tilde_j_plus_half = (u_tilde_j_int_left + u_tilde_j_int_right) / 2
        @. d_tilde_j_plus_half = u_tilde_j_int_right - u_tilde_j_int_left
        normalize!(d_tilde_j_plus_half)

        k += 2
    end

    u_hat[:, end] .= u_tilde[:, end]
    d_hat[:, end] .= d_tilde[:, end]

    return SmoothArcLengthInterpolation(u_hat, d_hat; shape_itp, kwargs...)
end

function SmoothArcLengthInterpolation(
        u::AbstractMatrix,
        d::AbstractMatrix;
        shape_itp::Union{AbstractInterpolation, Nothing} = nothing,
        extrapolation::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_left::ExtrapolationType.T = ExtrapolationType.None,
        extrapolation_right::ExtrapolationType.T = ExtrapolationType.None,
        cache_parameters::Bool = false,
        assume_linear_t = 1e-2,
        in_place::Bool = true)
    # Note: this method assumes that consecutive tangent lines are coplanar,
    # which is generally not the case for >2 dimensional points. Use the other constructors
    # to make sure this is satisfied.

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
        uⱼ = view(u, :, j)
        uⱼ₊₁ = view(u, :, j + 1)
        @. Δu = uⱼ₊₁ - uⱼ

        dⱼ = view(d, :, j)
        dⱼ₊₁ = view(d, :, j + 1)
        d_inner = dot(dⱼ, dⱼ₊₁)

        dⱼ_coef = (dot(Δu, dⱼ) - d_inner * dot(Δu, dⱼ₊₁)) / (1 - d_inner^2)
        @. u_int = uⱼ + dⱼ_coef * dⱼ

        dist₁ = euclidean(u_int, uⱼ)
        dist₂ = euclidean(u_int, uⱼ₊₁)

        δⱼ = if dist₁ < dist₂
            short_side_left[j] = true
            dist₁
        else
            dist₂
        end

        Rⱼ = δⱼ * sqrt((1 + d_inner) / (1 - d_inner))
        radius[j] = Rⱼ
        cⱼ = view(center, :, j)
        v₁ = view(dir_1, :, j)
        v₂ = view(dir_2, :, j)

        @. cⱼ = u_int + δⱼ * (dⱼ₊₁ - dⱼ) / (1 - d_inner)
        @. v₁ = -δⱼ * (dⱼ₊₁ - d_inner * dⱼ) / (1 - d_inner)
        @. v₂ = Rⱼ * dⱼ

        Δt_circle_seg = 2Rⱼ * atan(δⱼ, Rⱼ)
        Δt_line_seg = abs(dist₂ - dist₁)
        Δt_circle_segment[j] = Δt_circle_seg
        Δt_line_segment[j] = Δt_line_seg

        t[j + 1] = t[j] + Δt_circle_seg + Δt_line_seg
    end

    extrapolation_left, extrapolation_right = munge_extrapolation(
        extrapolation, extrapolation_left, extrapolation_right)
    linear_lookup = seems_linear(assume_linear_t, t)

    out = Vector{P}(undef, N)
    derivative = Vector{P}(undef, N)

    return SmoothArcLengthInterpolation(
        u, t, d, shape_itp, Δt_circle_segment, Δt_line_segment,
        center, radius, dir_1, dir_2, short_side_left,
        nothing, extrapolation_left, extrapolation_right, linear_lookup, out, derivative, in_place)
end
