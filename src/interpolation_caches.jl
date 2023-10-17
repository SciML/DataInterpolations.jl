### Linear Interpolation
struct LinearInterpolation{uType, tType, FT, T} <: AbstractInterpolation{FT, T}
    u::uType
    t::tType
    extrapolate::Bool
    function LinearInterpolation{FT}(u, t, extrapolate) where {FT}
        new{typeof(u), typeof(t), FT, eltype(u)}(u, t, extrapolate)
    end
end

function LinearInterpolation(u, t; extrapolate = false)
    u, t = munge_data(u, t)
    LinearInterpolation{true}(u, t, extrapolate)
end

### Quadratic Interpolation
struct QuadraticInterpolation{uType, tType, FT, T} <: AbstractInterpolation{FT, T}
    u::uType
    t::tType
    mode::Symbol
    extrapolate::Bool
    function QuadraticInterpolation{FT}(u, t, mode, extrapolate) where {FT}
        mode ∈ (:Forward, :Backward) ||
            error("mode should be :Forward or :Backward for QuadraticInterpolation")
        new{typeof(u), typeof(t), FT, eltype(u)}(u, t, mode, extrapolate)
    end
end

function QuadraticInterpolation(u, t, mode; extrapolate = false)
    u, t = munge_data(u, t)
    QuadraticInterpolation{true}(u, t, mode, extrapolate)
end

function QuadraticInterpolation(u, t; extrapolate = false)
    QuadraticInterpolation(u, t, :Forward; extrapolate)
end

### Lagrange Interpolation
struct LagrangeInterpolation{uType, tType, FT, T, bcacheType} <:
       AbstractInterpolation{FT, T}
    u::uType
    t::tType
    n::Int
    bcache::bcacheType
    extrapolate::Bool
    function LagrangeInterpolation{FT}(u, t, n, extrapolate) where {FT}
        bcache = zeros(eltype(u[1]), n + 1)
        fill!(bcache, NaN)
        new{typeof(u), typeof(t), FT, eltype(u), typeof(bcache)}(u,
            t,
            n,
            bcache,
            extrapolate)
    end
end

function LagrangeInterpolation(u, t, n = nothing; extrapolate = false)
    u, t = munge_data(u, t)
    if isnothing(n)
        n = length(t) - 1 # degree
    end
    if n != length(t) - 1
        error("Currently only n=length(t) - 1 is supported")
    end
    LagrangeInterpolation{true}(u, t, n, extrapolate)
end

### Akima Interpolation
struct AkimaInterpolation{uType, tType, bType, cType, dType, FT, T} <:
       AbstractInterpolation{FT, T}
    u::uType
    t::tType
    b::bType
    c::cType
    d::dType
    extrapolate::Bool
    function AkimaInterpolation{FT}(u, t, b, c, d, extrapolate) where {FT}
        new{typeof(u), typeof(t), typeof(b), typeof(c),
            typeof(d), FT, eltype(u)}(u,
            t,
            b,
            c,
            d,
            extrapolate)
    end
end

function AkimaInterpolation(u, t; extrapolate = false)
    u, t = munge_data(u, t)
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

    AkimaInterpolation{true}(u, t, b, c, d, extrapolate)
end

### ConstantInterpolation Interpolation
struct ConstantInterpolation{uType, tType, dirType, FT, T} <: AbstractInterpolation{FT, T}
    u::uType
    t::tType
    dir::Symbol # indicates if value to the $dir should be used for the interpolation
    extrapolate::Bool
    function ConstantInterpolation{FT}(u, t, dir, extrapolate) where {FT}
        new{typeof(u), typeof(t), typeof(dir), FT, eltype(u)}(u, t, dir, extrapolate)
    end
end

function ConstantInterpolation(u, t; dir = :left, extrapolate = false)
    u, t = munge_data(u, t)
    ConstantInterpolation{true}(u, t, dir, extrapolate)
end

Base.@deprecate_binding ZeroSpline ConstantInterpolation

### QuadraticSpline Interpolation
struct QuadraticSpline{uType, tType, tAType, dType, zType, FT, T} <:
       AbstractInterpolation{FT, T}
    u::uType
    t::tType
    tA::tAType
    d::dType
    z::zType
    extrapolate::Bool
    function QuadraticSpline{FT}(u, t, tA, d, z, extrapolate) where {FT}
        new{typeof(u), typeof(t), typeof(tA),
            typeof(d), typeof(z), FT, eltype(u)}(u,
            t,
            tA,
            d,
            z,
            extrapolate)
    end
end

function QuadraticSpline(u::uType,
    t;
    extrapolate = false) where {uType <: AbstractVector{<:Number}}
    u, t = munge_data(u, t)
    s = length(t)
    dl = ones(eltype(t), s - 1)
    d_tmp = ones(eltype(t), s)
    du = zeros(eltype(t), s - 1)
    tA = Tridiagonal(dl, d_tmp, du)

    # zero for element type of d, which we don't know yet
    typed_zero = zero(2 // 1 * (u[begin + 1] - u[begin]) / (t[begin + 1] - t[begin]))

    d = map(i -> i == 1 ? typed_zero : 2 // 1 * (u[i] - u[i - 1]) / (t[i] - t[i - 1]), 1:s)
    z = tA \ d
    QuadraticSpline{true}(u, t, tA, d, z, extrapolate)
end

function QuadraticSpline(u::uType, t; extrapolate = false) where {uType <: AbstractVector}
    u, t = munge_data(u, t)
    s = length(t)
    dl = ones(eltype(t), s - 1)
    d_tmp = ones(eltype(t), s)
    du = zeros(eltype(t), s - 1)
    tA = Tridiagonal(dl, d_tmp, du)
    d_ = map(i -> i == 1 ? zeros(eltype(t), size(u[1])) :
                  2 // 1 * (u[i] - u[i - 1]) / (t[i] - t[i - 1]),
        1:s)
    d = transpose(reshape(reduce(hcat, d_), :, s))
    z_ = reshape(transpose(tA \ d), size(u[1])..., :)
    z = [z_s for z_s in eachslice(z_, dims = ndims(z_))]
    QuadraticSpline{true}(u, t, tA, d, z, extrapolate)
end

# Cubic Spline Interpolation
struct CubicSpline{uType, tType, hType, zType, FT, T} <: AbstractInterpolation{FT, T}
    u::uType
    t::tType
    h::hType
    z::zType
    extrapolate::Bool
    function CubicSpline{FT}(u, t, h, z, extrapolate) where {FT}
        new{typeof(u), typeof(t), typeof(h), typeof(z), FT, eltype(u)}(u,
            t,
            h,
            z,
            extrapolate)
    end
end

function CubicSpline(u::uType,
    t;
    extrapolate = false) where {uType <: AbstractVector{<:Number}}
    u, t = munge_data(u, t)
    n = length(t) - 1
    h = vcat(0, map(k -> t[k + 1] - t[k], 1:(length(t) - 1)), 0)
    dl = h[2:(n + 1)]
    d_tmp = 2 .* (h[1:(n + 1)] .+ h[2:(n + 2)])
    du = h[2:(n + 1)]
    tA = Tridiagonal(dl, d_tmp, du)

    # zero for element type of d, which we don't know yet
    typed_zero = zero(6(u[begin + 2] - u[begin + 1]) / h[begin + 2] -
                      6(u[begin + 1] - u[begin]) / h[begin + 1])

    d = map(i -> i == 1 || i == n + 1 ? typed_zero :
                 6(u[i + 1] - u[i]) / h[i + 1] - 6(u[i] - u[i - 1]) / h[i],
        1:(n + 1))
    z = tA \ d
    CubicSpline{true}(u, t, h[1:(n + 1)], z, extrapolate)
end

function CubicSpline(u::uType, t; extrapolate = false) where {uType <: AbstractVector}
    u, t = munge_data(u, t)
    n = length(t) - 1
    h = vcat(0, map(k -> t[k + 1] - t[k], 1:(length(t) - 1)), 0)
    dl = h[2:(n + 1)]
    d_tmp = 2 .* (h[1:(n + 1)] .+ h[2:(n + 2)])
    du = h[2:(n + 1)]
    tA = Tridiagonal(dl, d_tmp, du)
    d_ = map(i -> i == 1 || i == n + 1 ? zeros(eltype(t), size(u[1])) :
                  6(u[i + 1] - u[i]) / h[i + 1] - 6(u[i] - u[i - 1]) / h[i],
        1:(n + 1))
    d = transpose(reshape(reduce(hcat, d_), :, n + 1))
    z_ = reshape(transpose(tA \ d), size(u[1])..., :)
    z = [z_s for z_s in eachslice(z_, dims = ndims(z_))]
    CubicSpline{true}(u, t, h[1:(n + 1)], z, extrapolate)
end

### BSpline Curve Interpolation
struct BSplineInterpolation{uType, tType, pType, kType, cType, FT, T} <:
       AbstractInterpolation{FT, T}
    u::uType
    t::tType
    d::Int    # degree
    p::pType  # params vector
    k::kType  # knot vector
    c::cType  # control points
    pVecType::Symbol
    knotVecType::Symbol
    extrapolate::Bool
    function BSplineInterpolation{FT}(u,
        t,
        d,
        p,
        k,
        c,
        pVecType,
        knotVecType,
        extrapolate) where {FT}
        new{typeof(u), typeof(t), typeof(p), typeof(k), typeof(c), FT, eltype(u)}(u,
            t,
            d,
            p,
            k,
            c,
            pVecType,
            knotVecType,
            extrapolate)
    end
end

function BSplineInterpolation(u, t, d, pVecType, knotVecType; extrapolate = false)
    u, t = munge_data(u, t)
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
    N = spline_coefficients(n, d, k, p)
    c = vec(N \ u[:, :])
    BSplineInterpolation{true}(u, t, d, p, k, c, pVecType, knotVecType, extrapolate)
end

### BSpline Curve Approx
struct BSplineApprox{uType, tType, pType, kType, cType, FT, T} <:
       AbstractInterpolation{FT, T}
    u::uType
    t::tType
    d::Int    # degree
    h::Int    # number of control points (n => h >= d >= 1)
    p::pType  # params vector
    k::kType  # knot vector
    c::cType  # control points
    pVecType::Symbol
    knotVecType::Symbol
    extrapolate::Bool
    function BSplineApprox{FT}(u,
        t,
        d,
        h,
        p,
        k,
        c,
        pVecType,
        knotVecType,
        extrapolate) where {FT}
        new{typeof(u), typeof(t), typeof(p), typeof(k), typeof(c), FT, eltype(u)}(u,
            t,
            d,
            h,
            p,
            k,
            c,
            pVecType,
            knotVecType,
            extrapolate)
    end
end

function BSplineApprox(u, t, d, h, pVecType, knotVecType; extrapolate = false)
    u, t = munge_data(u, t)
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
        N[i, :] .= spline_coefficients(h, d, k, p[i])
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
    BSplineApprox{true}(u, t, d, h, p, k, c, pVecType, knotVecType, extrapolate)
end
