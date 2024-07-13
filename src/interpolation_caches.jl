"""
    LinearInterpolation(u, t; extrapolate = false)

It is the method of interpolating between the data points using a linear polynomial. For any point, two data points one each side are chosen and connected with a line.
Extrapolation extends the last linear polynomial on each side.

## Arguments

  - `u`: data points.
  - `t`: time points.

## Keyword Arguments

  - `extrapolate`: boolean value to allow extrapolation. Defaults to `false`.
  - `safetycopy`: boolean value to make a copy of `u` and `t`. Defaults to `true`.
"""
struct LinearInterpolation{uType, tType, IType, pType, T} <: AbstractInterpolation{T}
    u::uType
    t::tType
    I::IType
    p::LinearParameterCache{pType}
    extrapolate::Bool
    idx_prev::Base.RefValue{Int}
    safetycopy::Bool
    function LinearInterpolation(u, t, I, p, extrapolate, safetycopy)
        new{typeof(u), typeof(t), typeof(I), typeof(p.slope), eltype(u)}(
            u, t, I, p, extrapolate, Ref(1), safetycopy)
    end
end

function LinearInterpolation(u, t; extrapolate = false, safetycopy = true)
    u, t = munge_data(u, t, safetycopy)
    p = LinearParameterCache(u, t)
    A = LinearInterpolation(u, t, nothing, p, extrapolate, safetycopy)
    I = cumulative_integral(A, A.u)
    LinearInterpolation(u, t, I, p, extrapolate, safetycopy)
end

"""
    QuadraticInterpolation(u, t, mode = :Forward; extrapolate = false)

It is the method of interpolating between the data points using quadratic polynomials. For any point, three data points nearby are taken to fit a quadratic polynomial.
Extrapolation extends the last quadratic polynomial on each side.

## Arguments

  - `u`: data points.
  - `t`: time points.
  - `mode`: `:Forward` or `:Backward`. If `:Forward`, two data points ahead of the point and one data point behind is taken for interpolation. If `:Backward`, two data points behind and one ahead is taken for interpolation.

## Keyword Arguments

  - `extrapolate`: boolean value to allow extrapolation. Defaults to `false`.
  - `safetycopy`: boolean value to make a copy of `u` and `t`. Defaults to `true`.
"""
struct QuadraticInterpolation{uType, tType, IType, pType, T} <: AbstractInterpolation{T}
    u::uType
    t::tType
    I::IType
    p::QuadraticParameterCache{pType}
    mode::Symbol
    extrapolate::Bool
    idx_prev::Base.RefValue{Int}
    safetycopy::Bool
    function QuadraticInterpolation(u, t, I, p, mode, extrapolate, safetycopy)
        mode ∈ (:Forward, :Backward) ||
            error("mode should be :Forward or :Backward for QuadraticInterpolation")
        new{typeof(u), typeof(t), typeof(I), typeof(p.l₀), eltype(u)}(
            u, t, I, p, mode, extrapolate, Ref(1), safetycopy)
    end
end

function QuadraticInterpolation(u, t, mode; extrapolate = false, safetycopy = true)
    u, t = munge_data(u, t, safetycopy)
    p = QuadraticParameterCache(u, t)
    A = QuadraticInterpolation(u, t, nothing, p, mode, extrapolate, safetycopy)
    I = cumulative_integral(A, A.u)
    QuadraticInterpolation(u, t, I, p, mode, extrapolate, safetycopy)
end

function QuadraticInterpolation(u, t; extrapolate = false, safetycopy = true)
    QuadraticInterpolation(u, t, :Forward; extrapolate, safetycopy)
end

"""
    LagrangeInterpolation(u, t, n = length(t) - 1; extrapolate = false)

It is the method of interpolation using Lagrange polynomials of (k-1)th order passing through all the data points where k is the number of data points.

## Arguments

  - `u`: data points.
  - `t`: time points.
  - `n`: order of the polynomial. Currently only (k-1)th order where k is the number of data points.

## Keyword Arguments

  - `extrapolate`: boolean value to allow extrapolation. Defaults to `false`.
  - `safetycopy`: boolean value to make a copy of `u` and `t`. Defaults to `true`.
"""
struct LagrangeInterpolation{uType, tType, T, bcacheType} <:
       AbstractInterpolation{T}
    u::uType
    t::tType
    n::Int
    bcache::bcacheType
    idxs::Vector{Int}
    extrapolate::Bool
    idx_prev::Base.RefValue{Int}
    safetycopy::Bool
    function LagrangeInterpolation(u, t, n, extrapolate, safetycopy)
        bcache = zeros(eltype(u[1]), n + 1)
        idxs = zeros(Int, n + 1)
        fill!(bcache, NaN)
        new{typeof(u), typeof(t), eltype(u), typeof(bcache)}(u,
            t,
            n,
            bcache,
            idxs,
            extrapolate,
            Ref(1),
            safetycopy
        )
    end
end

function LagrangeInterpolation(
        u, t, n = length(t) - 1; extrapolate = false, safetycopy = true)
    u, t = munge_data(u, t, safetycopy)
    if n != length(t) - 1
        error("Currently only n=length(t) - 1 is supported")
    end
    LagrangeInterpolation(u, t, n, extrapolate, safetycopy)
end

"""
    AkimaInterpolation(u, t; extrapolate = false)

It is a spline interpolation built from cubic polynomials. It forms a continuously differentiable function. For more details, refer: https://en.wikipedia.org/wiki/Akima_spline.
Extrapolation extends the last cubic polynomial on each side.

## Arguments

  - `u`: data points.
  - `t`: time points.

## Keyword Arguments

  - `extrapolate`: boolean value to allow extrapolation. Defaults to `false`.
  - `safetycopy`: boolean value to make a copy of `u` and `t`. Defaults to `true`.
"""
struct AkimaInterpolation{uType, tType, IType, bType, cType, dType, T} <:
       AbstractInterpolation{T}
    u::uType
    t::tType
    I::IType
    b::bType
    c::cType
    d::dType
    extrapolate::Bool
    idx_prev::Base.RefValue{Int}
    safetycopy::Bool
    function AkimaInterpolation(u, t, I, b, c, d, extrapolate, safetycopy)
        new{typeof(u), typeof(t), typeof(I), typeof(b), typeof(c),
            typeof(d), eltype(u)}(u,
            t,
            I,
            b,
            c,
            d,
            extrapolate,
            Ref(1),
            safetycopy
        )
    end
end

function AkimaInterpolation(u, t; extrapolate = false, safetycopy = true)
    u, t = munge_data(u, t, safetycopy)
    n = length(t)
    dt = diff(t)
    m = Array{eltype(u)}(undef, n + 3)
    m[3:(end - 2)] = diff(u) ./ dt
    m[2] = 2m[3] - m[4]
    m[1] = 2m[2] - m[3]
    m[end - 1] = 2m[end - 2] - m[end - 3]
    m[end] = 2m[end - 1] - m[end - 2]

    b = 0.5 .* (m[4:end] .+ m[1:(end - 3)])
    dm = abs.(diff(m))
    f1 = dm[3:(n + 2)]
    f2 = dm[1:n]
    f12 = f1 + f2
    ind = findall(f12 .> 1e-9 * maximum(f12))
    b[ind] = (f1[ind] .* m[ind .+ 1] .+
              f2[ind] .* m[ind .+ 2]) ./ f12[ind]
    c = (3.0 .* m[3:(end - 2)] .- 2.0 .* b[1:(end - 1)] .- b[2:end]) ./ dt
    d = (b[1:(end - 1)] .+ b[2:end] .- 2.0 .* m[3:(end - 2)]) ./ dt .^ 2

    A = AkimaInterpolation(u, t, nothing, b, c, d, extrapolate, safetycopy)
    I = cumulative_integral(A, A.u)
    AkimaInterpolation(u, t, I, b, c, d, extrapolate, safetycopy)
end

"""
    ConstantInterpolation(u, t; dir = :left, extrapolate = false)

It is the method of interpolating using a constant polynomial. For any point, two adjacent data points are found on either side (left and right). The value at that point depends on `dir`.
If it is `:left`, then the value at the left point is chosen and if it is `:right`, the value at the right point is chosen.
Extrapolation extends the last constant polynomial at the end points on each side.

## Arguments

  - `u`: data points.
  - `t`: time points.

## Keyword Arguments

  - `dir`: indicates which value should be used for interpolation (`:left` or `:right`).
  - `extrapolate`: boolean value to allow extrapolation. Defaults to `false`.
  - `safetycopy`: boolean value to make a copy of `u` and `t`. Defaults to `true`.
"""
struct ConstantInterpolation{uType, tType, IType, T} <: AbstractInterpolation{T}
    u::uType
    t::tType
    I::IType
    p::Nothing
    dir::Symbol # indicates if value to the $dir should be used for the interpolation
    extrapolate::Bool
    idx_prev::Base.RefValue{Int}
    safetycopy::Bool
    function ConstantInterpolation(u, t, I, dir, extrapolate, safetycopy)
        new{typeof(u), typeof(t), typeof(I), eltype(u)}(
            u, t, I, nothing, dir, extrapolate, Ref(1), safetycopy)
    end
end

function ConstantInterpolation(u, t; dir = :left, extrapolate = false, safetycopy = true)
    u, t = munge_data(u, t, safetycopy)
    A = ConstantInterpolation(u, t, nothing, dir, extrapolate, safetycopy)
    I = cumulative_integral(A, A.u)
    ConstantInterpolation(u, t, I, dir, extrapolate, safetycopy)
end

"""
    QuadraticSpline(u, t; extrapolate = false)

It is a spline interpolation using piecewise quadratic polynomials between each pair of data points. Its first derivative is also continuous.
Extrapolation extends the last quadratic polynomial on each side.

## Arguments

  - `u`: data points.
  - `t`: time points.

## Keyword Arguments

  - `extrapolate`: boolean value to allow extrapolation. Defaults to `false`.
  - `safetycopy`: boolean value to make a copy of `u` and `t`. Defaults to `true`.
"""
struct QuadraticSpline{uType, tType, IType, pType, tAType, dType, zType, T} <:
       AbstractInterpolation{T}
    u::uType
    t::tType
    I::IType
    p::QuadraticSplineParameterCache{tAType, dType, zType, pType}
    extrapolate::Bool
    idx_prev::Base.RefValue{Int}
    safetycopy::Bool
    function QuadraticSpline(u, t, I, p, extrapolate, safetycopy)
        new{typeof(u), typeof(t), typeof(I), typeof(p.σ), typeof(p.tA),
            typeof(p.d), typeof(p.z), eltype(u)}(u,
            t,
            I,
            p,
            extrapolate,
            Ref(1),
            safetycopy
        )
    end
end

function QuadraticSpline(
        u, t; extrapolate = false,
        safetycopy = true)
    u, t = munge_data(u, t, safetycopy)
    p = QuadraticSplineParameterCache(u, t)
    A = QuadraticSpline(u, t, nothing, p, extrapolate, safetycopy)
    I = cumulative_integral(A, A.u)
    QuadraticSpline(u, t, I, p, extrapolate, safetycopy)
end

"""
    CubicSpline(u, t; extrapolate = false)

It is a spline interpolation using piecewise cubic polynomials between each pair of data points. Its first and second derivative is also continuous.
Second derivative on both ends are zero, which are also called "natural" boundary conditions. Extrapolation extends the last cubic polynomial on each side.

## Arguments

  - `u`: data points.
  - `t`: time points.

## Keyword Arguments

  - `extrapolate`: boolean value to allow extrapolation. Defaults to `false`.
  - `safetycopy`: boolean value to make a copy of `u` and `t`. Defaults to `true`.
"""
struct CubicSpline{uType, tType, IType, pType, hType, zType, T} <: AbstractInterpolation{T}
    u::uType
    t::tType
    I::IType
    p::CubicSplineParameterCache{pType}
    h::hType
    z::zType
    extrapolate::Bool
    idx_prev::Base.RefValue{Int}
    safetycopy::Bool
    function CubicSpline(u, t, I, p, h, z, extrapolate, safetycopy)
        new{typeof(u), typeof(t), typeof(I), typeof(p.c₁), typeof(h), typeof(z), eltype(u)}(
            u,
            t,
            I,
            p,
            h,
            z,
            extrapolate,
            Ref(1),
            safetycopy
        )
    end
end

function CubicSpline(u::uType,
        t;
        extrapolate = false, safetycopy = true) where {uType <: AbstractVector{<:Number}}
    u, t = munge_data(u, t, safetycopy)
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
    p = CubicSplineParameterCache(u, h, z)
    A = CubicSpline(u, t, nothing, p, h[1:(n + 1)], z, extrapolate, safetycopy)
    I = cumulative_integral(A, A.u)
    CubicSpline(u, t, I, p, h[1:(n + 1)], z, extrapolate, safetycopy)
end

function CubicSpline(
        u::uType, t; extrapolate = false, safetycopy = true) where {uType <: AbstractVector}
    u, t = munge_data(u, t, safetycopy)
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
    p = CubicSplineParameterCache(u, h, z)
    A = CubicSpline(u, t, nothing, p, h[1:(n + 1)], z, extrapolate, safetycopy)
    I = cumulative_integral(A, A.u)
    CubicSpline(u, t, I, p, h[1:(n + 1)], z, extrapolate, safetycopy)
end

"""
    BSplineInterpolation(u, t, d, pVecType, knotVecType; extrapolate = false)

It is a curve defined by the linear combination of `n` basis functions of degree `d` where `n` is the number of data points. For more information, refer https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-curve.html.
Extrapolation is a constant polynomial of the end points on each side.

## Arguments

  - `u`: data points.
  - `t`: time points.
  - `d`: degree of the piecewise polynomial.
  - `pVecType`: symbol to parameters vector, `:Uniform` for uniform spaced parameters and `:ArcLen` for parameters generated by chord length method.
  - `knotVecType`: symbol to knot vector, `:Uniform` for uniform knot vector, `:Average` for average spaced knot vector.

## Keyword Arguments

  - `extrapolate`: boolean value to allow extrapolation. Defaults to `false`.
  - `safetycopy`: boolean value to make a copy of `u` and `t`. Defaults to `true`.
"""
struct BSplineInterpolation{uType, tType, pType, kType, cType, NType, T} <:
       AbstractInterpolation{T}
    u::uType
    t::tType
    d::Int    # degree
    p::pType  # params vector
    k::kType  # knot vector
    c::cType  # control points
    N::NType  # Spline coefficients (preallocated memory)
    pVecType::Symbol
    knotVecType::Symbol
    extrapolate::Bool
    idx_prev::Base.RefValue{Int}
    safetycopy::Bool
    function BSplineInterpolation(u,
            t,
            d,
            p,
            k,
            c,
            N,
            pVecType,
            knotVecType,
            extrapolate,
            safetycopy)
        new{typeof(u), typeof(t), typeof(p), typeof(k), typeof(c), typeof(N), eltype(u)}(u,
            t,
            d,
            p,
            k,
            c,
            N,
            pVecType,
            knotVecType,
            extrapolate,
            Ref(1),
            safetycopy
        )
    end
end

function BSplineInterpolation(
        u, t, d, pVecType, knotVecType; extrapolate = false, safetycopy = true)
    u, t = munge_data(u, t, safetycopy)
    n = length(t)
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
    N = zeros(eltype(t), n, n)
    spline_coefficients!(N, d, k, p)
    c = vec(N \ u[:, :])
    N = zeros(eltype(t), n)
    BSplineInterpolation(
        u, t, d, p, k, c, N, pVecType, knotVecType, extrapolate, safetycopy)
end

"""
    BSplineApprox(u, t, d, h, pVecType, knotVecType; extrapolate = false)

It is a regression based B-spline. The argument choices are the same as the `BSplineInterpolation`, with the additional parameter `h < length(t)` which is the number of control points to use, with smaller `h` indicating more smoothing.
For more information, refer http://www.cad.zju.edu.cn/home/zhx/GM/009/00-bsia.pdf.
Extrapolation is a constant polynomial of the end points on each side.

## Arguments

  - `u`: data points.
  - `t`: time points.
  - `d`: degree of the piecewise polynomial.
  - `h`: number of control points to use.
  - `pVecType`: symbol to parameters vector, `:Uniform` for uniform spaced parameters and `:ArcLen` for parameters generated by chord length method.
  - `knotVecType`: symbol to knot vector, `:Uniform` for uniform knot vector, `:Average` for average spaced knot vector.

## Keyword Arguments

  - `extrapolate`: boolean value to allow extrapolation. Defaults to `false`.
  - `safetycopy`: boolean value to make a copy of `u` and `t`. Defaults to `true`.
"""
struct BSplineApprox{uType, tType, pType, kType, cType, NType, T} <:
       AbstractInterpolation{T}
    u::uType
    t::tType
    d::Int    # degree
    h::Int    # number of control points (n => h >= d >= 1)
    p::pType  # params vector
    k::kType  # knot vector
    c::cType  # control points
    N::NType  # Spline coefficients (preallocated memory)
    pVecType::Symbol
    knotVecType::Symbol
    extrapolate::Bool
    idx_prev::Base.RefValue{Int}
    safetycopy::Bool
    function BSplineApprox(u,
            t,
            d,
            h,
            p,
            k,
            c,
            N,
            pVecType,
            knotVecType,
            extrapolate,
            safetycopy
    )
        new{typeof(u), typeof(t), typeof(p), typeof(k), typeof(c), typeof(N), eltype(u)}(u,
            t,
            d,
            h,
            p,
            k,
            c,
            N,
            pVecType,
            knotVecType,
            extrapolate,
            Ref(1),
            safetycopy::Bool
        )
    end
end

function BSplineApprox(
        u, t, d, h, pVecType, knotVecType; extrapolate = false, safetycopy = true)
    u, t = munge_data(u, t, safetycopy)
    n = length(t)
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
    N = zeros(eltype(t), n, h)
    for i in 1:n
        spline_coefficients!(view(N, i, :), d, k, p[i])
    end
    for k in 2:(n - 1)
        q[k] = u[k] - N[k, 1] * u[1] - N[k, h] * u[end]
    end
    Q = Matrix{eltype(u)}(undef, h - 2, 1)
    for i in 2:(h - 1)
        s = 0.0
        for k in 2:(n - 1)
            s += N[k, i] * q[k]
        end
        Q[i - 1] = s
    end
    N = N[2:(end - 1), 2:(h - 1)]
    M = transpose(N) * N
    P = M \ Q
    c[2:(end - 1)] .= vec(P)
    N = zeros(eltype(t), h)
    BSplineApprox(u, t, d, h, p, k, c, N, pVecType, knotVecType, extrapolate, safetycopy)
end

"""
    CubicHermiteSpline(du, u, t; extrapolate = false)

It is a Cubic Hermite interpolation, which is a piece-wise third degree polynomial such that the value and the first derivative are equal to given values in the data points.

## Arguments

  - `du`: the derivative at the data points.
  - `u`: data points.
  - `t`: time points.

## Keyword Arguments

  - `extrapolate`: boolean value to allow extrapolation. Defaults to `false`.
  - `safetycopy`: boolean value to make a copy of `u` and `t`. Defaults to `true`.
"""
struct CubicHermiteSpline{uType, tType, IType, duType, pType, T} <: AbstractInterpolation{T}
    du::duType
    u::uType
    t::tType
    I::IType
    p::CubicHermiteParameterCache{pType}
    extrapolate::Bool
    idx_prev::Base.RefValue{Int}
    safetycopy::Bool
    function CubicHermiteSpline(du, u, t, I, p, extrapolate, safetycopy)
        new{typeof(u), typeof(t), typeof(I), typeof(du), typeof(p.c₁), eltype(u)}(
            du, u, t, I, p, extrapolate, Ref(1), safetycopy)
    end
end

function CubicHermiteSpline(du, u, t; extrapolate = false, safetycopy = true)
    @assert length(u)==length(du) "Length of `u` is not equal to length of `du`."
    u, t = munge_data(u, t, safetycopy)
    p = CubicHermiteParameterCache(du, u, t)
    A = CubicHermiteSpline(du, u, t, nothing, p, extrapolate, safetycopy)
    I = cumulative_integral(A, A.u)
    CubicHermiteSpline(du, u, t, I, p, extrapolate, safetycopy)
end

"""
    QuinticHermiteSpline(ddu, du, u, t; extrapolate = false)

It is a Quintic Hermite interpolation, which is a piece-wise fifth degree polynomial such that the value and the first and second derivative are equal to given values in the data points.

## Arguments

  - `ddu`: the second derivative at the data points.
  - `du`: the derivative at the data points.
  - `u`: data points.
  - `t`: time points.

## Keyword Arguments

  - `extrapolate`: boolean value to allow extrapolation. Defaults to `false`.
  - `safetycopy`: boolean value to make a copy of `u` and `t`. Defaults to `true`.
"""
struct QuinticHermiteSpline{uType, tType, IType, duType, dduType, pType, T} <:
       AbstractInterpolation{T}
    ddu::dduType
    du::duType
    u::uType
    t::tType
    I::IType
    p::QuinticHermiteParameterCache{pType}
    extrapolate::Bool
    idx_prev::Base.RefValue{Int}
    safetycopy::Bool
    function QuinticHermiteSpline(ddu, du, u, t, I, p, extrapolate, safetycopy)
        new{typeof(u), typeof(t), typeof(I), typeof(du),
            typeof(ddu), typeof(p.c₁), eltype(u)}(
            ddu, du, u, t, I, p, extrapolate, Ref(1), safetycopy)
    end
end

function QuinticHermiteSpline(ddu, du, u, t; extrapolate = false, safetycopy = true)
    @assert length(u)==length(du)==length(ddu) "Length of `u` is not equal to length of `du` or `ddu`."
    u, t = munge_data(u, t, safetycopy)
    p = QuinticHermiteParameterCache(ddu, du, u, t)
    A = QuinticHermiteSpline(ddu, du, u, t, nothing, p, extrapolate, safetycopy)
    I = cumulative_integral(A, A.u)
    QuinticHermiteSpline(ddu, du, u, t, I, p, extrapolate, safetycopy)
end
