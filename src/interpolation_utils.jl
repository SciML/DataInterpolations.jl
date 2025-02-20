function findRequiredIdxs!(A::LagrangeInterpolation, t, idx)
    n = length(A.t) - 1
    i_min, idx_min, idx_max = if t == A.t[idx]
        A.idxs[1] = idx
        2, idx, idx
    else
        1, idx + 1, idx
    end
    for i in i_min:(n + 1)
        if idx_min == 1
            A.idxs[i:end] .= range(idx_max + 1, idx_max + (n + 2 - i))
            break
        elseif idx_max == length(A.t)
            A.idxs[i:end] .= (idx_min - 1):-1:(idx_min - (n + 2 - i))
            break
        else
            left_diff = abs(t - A.t[idx_min - 1])
            right_diff = abs(t - A.t[idx_max + 1])
            left_expand = left_diff <= right_diff
        end
        if left_expand
            idx_min -= 1
            A.idxs[i] = idx_min
        else
            idx_max += 1
            A.idxs[i] = idx_max
        end
    end
    return idx
end

function spline_coefficients!(N, d, k, u::Number)
    N .= zero(u)
    if u == k[1]
        N[1] = one(u)
        return 1:1
    elseif u == k[end]
        N[end] = one(u)
        return length(N):length(N)
    else
        i = findfirst(x -> x > u, k) - 1
        N[i] = one(u)
        for deg in 1:d
            N[i - deg] = (k[i + 1] - u) / (k[i + 1] - k[i - deg + 1]) * N[i - deg + 1]
            for j in (i - deg + 1):(i - 1)
                N[j] = (u - k[j]) / (k[j + deg] - k[j]) * N[j] +
                       (k[j + deg + 1] - u) / (k[j + deg + 1] - k[j + 1]) * N[j + 1]
            end
            N[i] = (u - k[i]) / (k[i + deg] - k[i]) * N[i]
        end
        return (i - d):i
    end
end

function spline_coefficients!(N, d, k, u::AbstractVector)
    for i in 1:size(N)[2]
        spline_coefficients!(view(N, i, :), d, k, u[i])
    end
    return nothing
end

# Get the number of dimensions `ndims(interp(x))` of the interpolation `interp` evaluated at a single input `x`
# It is derived from the set of values `u` at the interpolation nodes
output_ndims(::AbstractVector) = 0 # each value is a scalar
output_ndims(::AbstractVector{<:AbstractArray{<:Any,N}}) where {N} = N # each value is an array but values are not stacked
output_ndims(::AbstractArray{<:Any,N}) where {N} = N - 1 # each value is an array but multiple values are stacked

function quadratic_spline_params(t::AbstractVector, sc::AbstractVector)

    # Create knot vector
    # Don't use x[end-1] as knot to match number of degrees of freedom with data
    k = zeros(length(t) + 3)
    k[1:3] .= t[1]
    k[(end - 2):end] .= t[end]
    k[4:(end - 3)] .= t[2:(end - 2)]

    # Create linear system Ac = u, where:
    # - A consists of basis function evaulations in t
    # - c are 1D control points 
    n = length(t)
    dtype_sc = typeof(t[1] / t[1])

    diag = Vector{dtype_sc}(undef, n)
    diag_hi = Vector{dtype_sc}(undef, n - 1)
    diag_lo = Vector{dtype_sc}(undef, n - 1)

    for (i, tᵢ) in enumerate(t)
        spline_coefficients!(sc, 2, k, tᵢ)
        diag[i] = sc[i]
        (i > 1) && (diag_lo[i - 1] = sc[i - 1])
        (i < n) && (diag_hi[i] = sc[i + 1])
    end

    A = Tridiagonal(diag_lo, diag, diag_hi)

    return k, A
end

# helper function for data manipulation
function munge_data(u::AbstractVector, t::AbstractVector)
    Tu = nonmissingtype(eltype(u))
    Tt = nonmissingtype(eltype(t))
    if Tu === eltype(u) && Tt === eltype(t)
        return u, t
    end

    @assert length(t) == length(u)
    non_missing_mask = map((ui, ti) -> !ismissing(ui) && !ismissing(ti), u, t)
    u = convert(AbstractVector{Tu}, u[non_missing_mask])
    t = convert(AbstractVector{Tt}, t[non_missing_mask])

    return u, t
end

function munge_data(U::AbstractMatrix, t::AbstractVector)
    TU = nonmissingtype(eltype(U))
    Tt = nonmissingtype(eltype(t))
    if TU === eltype(U) && Tt === eltype(t)
        return U, t
    end

    @assert length(t) == size(U, 2)
    non_missing_mask = map(
        (uis, ti) -> !any(ismissing, uis) && !ismissing(ti), eachcol(U), t)
    U = convert(AbstractMatrix{TU}, U[:, non_missing_mask])
    t = convert(AbstractVector{Tt}, t[non_missing_mask])

    return U, t
end

function munge_data(U::AbstractArray{T, N}, t) where {T, N}
    TU = nonmissingtype(eltype(U))
    Tt = nonmissingtype(eltype(t))
    if TU === eltype(U) && Tt === eltype(t)
        return U, t
    end

    @assert length(t) == size(U, N)
    non_missing_mask = map(
        (uis, ti) -> !any(ismissing, uis) && !ismissing(ti), eachslice(U; dims = N), t)
    U = convert(AbstractArray{TU, N}, copy(selectdim(U, N, non_missing_mask)))
    t = convert(AbstractVector{Tt}, t[non_missing_mask])

    return U, t
end

seems_linear(assume_linear_t::Bool, _) = assume_linear_t
seems_linear(assume_linear_t::Number, t) = looks_linear(t; threshold = assume_linear_t)

"""
    looks_linear(t; threshold = 1e-2)

Determine if the abscissae `t` are regularly distributed, taking the standard deviation of
the difference between the array of abscissae with respect to the straight line linking
its first and last elements, normalized by the range of `t`. If this standard deviation is
below the given `threshold`, the vector looks linear (return true). Internal function -
interface may change.
"""
function looks_linear(t; threshold = 1e-2)
    length(t) <= 2 && return true
    t_0, t_f = first(t), last(t)
    t_span = t_f - t_0
    tspan_over_N = t_span * length(t)^(-1)
    norm_var = sum(
        (t_i - t_0 - i * tspan_over_N)^2 for (i, t_i) in enumerate(t)
    ) / (length(t) * t_span^2)
    norm_var < threshold^2
end

function get_idx(A::AbstractInterpolation, t, iguess::Union{<:Integer, Guesser}; lb = 1,
        ub_shift = -1, idx_shift = 0, side = :last)
    tvec = A.t
    ub = length(tvec) + ub_shift
    return if side == :last
        clamp(searchsortedlastcorrelated(tvec, t, iguess) + idx_shift, lb, ub)
    elseif side == :first
        clamp(searchsortedfirstcorrelated(tvec, t, iguess) + idx_shift, lb, ub)
    else
        error("side must be :first or :last")
    end
end

cumulative_integral(::AbstractInterpolation, ::Bool) = nothing
function cumulative_integral(A::AbstractInterpolation{<:Number}, cache_parameters::Bool)
    Base.require_one_based_indexing(A.u)
    idxs = cache_parameters ? (1:(length(A.t) - 1)) : (1:0)
    return cumsum(_integral(A, idx, t1, t2)
    for (idx, t1, t2) in zip(idxs, @view(A.t[begin:(end - 1)]), @view(A.t[(begin + 1):end])))
end

function get_parameters(A::LinearInterpolation, idx)
    if A.cache_parameters
        A.p.slope[idx]
    else
        linear_interpolation_parameters(A.u, A.t, idx)
    end
end

function get_parameters(A::QuadraticInterpolation, idx)
    if A.cache_parameters
        A.p.α[idx], A.p.β[idx]
    else
        quadratic_interpolation_parameters(A.u, A.t, idx, A.mode)
    end
end

function get_parameters(A::QuadraticSpline, idx)
    if A.cache_parameters
        A.p.α[idx], A.p.β[idx]
    else
        quadratic_spline_parameters(A.u, A.t, A.k, A.c, A.sc, idx)
    end
end

function get_parameters(A::CubicSpline, idx)
    if A.cache_parameters
        A.p.c₁[idx], A.p.c₂[idx]
    else
        cubic_spline_parameters(A.u, A.h, A.z, idx)
    end
end

function get_parameters(A::CubicHermiteSpline, idx)
    if A.cache_parameters
        A.p.c₁[idx], A.p.c₂[idx]
    else
        cubic_hermite_spline_parameters(A.du, A.u, A.t, idx)
    end
end

function get_parameters(A::QuinticHermiteSpline, idx)
    if A.cache_parameters
        A.p.c₁[idx], A.p.c₂[idx], A.p.c₃[idx]
    else
        quintic_hermite_spline_parameters(A.ddu, A.du, A.u, A.t, idx)
    end
end

function du_PCHIP(u, t)
    h = diff(u)
    δ = h ./ diff(t)
    s = sign.(δ)

    function _du(k)
        sₖ₋₁, sₖ = if k == 1
            s[1], s[2]
        elseif k == lastindex(t)
            s[end - 1], s[end]
        else
            s[k - 1], s[k]
        end

        if sₖ₋₁ == 0 && sₖ == 0
            zero(eltype(δ))
        elseif sₖ₋₁ == sₖ
            if k == 1
                ((2 * h[1] + h[2]) * δ[1] - h[1] * δ[2]) / (h[1] + h[2])
            elseif k == lastindex(t)
                ((2 * h[end] + h[end - 1]) * δ[end] - h[end] * δ[end - 1]) /
                (h[end] + h[end - 1])
            else
                w₁ = 2h[k] + h[k - 1]
                w₂ = h[k] + 2h[k - 1]
                δ[k - 1] * δ[k] * (w₁ + w₂) / (w₁ * δ[k] + w₂ * δ[k - 1])
            end
        else
            zero(eltype(δ))
        end
    end

    return _du.(eachindex(t))
end

function integrate_cubic_polynomial(t1, t2, offset, a, b, c, d)
    t1_rel = t1 - offset
    t2_rel = t2 - offset
    t_sum = t1_rel + t2_rel
    t_sq_sum = t1_rel^2 + t2_rel^2
    Δt = t2 - t1
    Δt * (a + t_sum * (b / 2 + d * t_sq_sum / 4) + c * (t_sq_sum + t1_rel * t2_rel) / 3)
end

function integrate_quintic_polynomial(t1, t2, offset, a, b, c, d, e, f)
    t1_rel = t1 - offset
    t2_rel = t2 - offset
    t_sum = t1_rel + t2_rel
    t_sq_sum = t1_rel^2 + t2_rel^2
    t_cb_sum = t1_rel^3 + t2_rel^3
    Δt = t2 - t1
    cube_diff_factor = t_sq_sum + t1_rel * t2_rel
    Δt * (a + t_sum * (b / 2 + d * t_sq_sum / 4) +
     cube_diff_factor * (c / 3 + f * t_cb_sum / 6)) +
    e * (t2_rel^5 - t1_rel^5) / 5
end

function munge_extrapolation(extrapolation, extrapolation_left, extrapolation_right)
    if extrapolation == ExtrapolationType.None
        extrapolation_left, extrapolation_right
    else
        extrapolation, extrapolation
    end
end

function transformation_periodic(A, t)
    Δt = last(A.t) - first(A.t)
    n, t_ = fldmod(t - first(A.t), Δt)
    t_ += first(A.t)
    (n > 0) && (n -= 1)
    t_, n
end

function transformation_reflective(A, t)
    Δt = last(A.t) - first(A.t)
    n, t_ = fldmod(t - first(A.t), Δt)
    t_ = isodd(n) ? last(A.t) - t_ : first(A.t) + t_
    (n > 0) && (n -= 1)
    t_, n
end

typed_nan(::AbstractArray{T}) where {T <: AbstractFloat} = T(NaN)
typed_nan(::AbstractArray{T}) where {T <: Integer} = zero(T)
