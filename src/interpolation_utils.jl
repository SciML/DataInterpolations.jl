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
    # `N` is zeroed because BSpline derivative paths read the full vector
    # (see `_derivative(::BSplineInterpolation, …)` in `derivatives.jl`).
    # Positions outside the body's `(i-d):i` write window must be zero or
    # stale values from previous calls would leak in.
    N .= zero(u)
    if u == k[1]
        N[1] = one(u)
        return 1:1
    elseif u == k[end]
        N[end] = one(u)
        return length(N):length(N)
    else
        # `k` is sorted; the legacy `findfirst(x -> x > u, k) - 1` did an O(n)
        # linear scan. `searchsortedlast` returns the same index in O(log n).
        i = searchsortedlast(k, u)
        return _spline_coefficients_body!(N, d, k, u, i)
    end
end

# Body of `spline_coefficients!` after the locator index `i` has been
# determined. Used by `spline_coefficients!` (logarithmic locator) and by
# `quadratic_spline_params` (O(1) amortised running locator).
function _spline_coefficients_body!(N, d, k, u, i)
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

function spline_coefficients!(N, d, k, u::AbstractVector)
    for i in 1:size(N)[2]
        spline_coefficients!(view(N, i, :), d, k, u[i])
    end
    return nothing
end

function quadratic_spline_params(t::AbstractVector, sc::AbstractVector)
    # Duplicate time points make the collocation system singular
    if any(i -> t[i] == t[i + 1], 1:(length(t) - 1))
        throw(
            ArgumentError(
                "The time points `t` must be unique for `QuadraticSpline`, but duplicate values were found."
            )
        )
    end

    # Create knot vector
    # Don't use x[end-1] as knot to match number of degrees of freedom with data
    k = zeros(length(t) + 3)
    k[1:3] .= t[1]
    k[(end - 2):end] .= t[end]
    k[4:(end - 3)] .= t[2:(end - 2)]

    # Create linear system Ac = u, where:
    # - A consists of basis function evaluations in t
    # - c are 1D control points
    n = length(t)
    dtype_sc = typeof(one(eltype(t)) / one(eltype(t)))

    diag = Vector{dtype_sc}(undef, n)
    diag_hi = Vector{dtype_sc}(undef, n - 1)
    diag_lo = Vector{dtype_sc}(undef, n - 1)

    # `t` is sorted and `k` is built from `t`, so the locator
    # `searchsortedlast(k, tᵢ)` is non-decreasing in `i`. Maintain a running
    # pointer to advance amortised O(1) per knot — total O(n) instead of the
    # O(n²) `findfirst` scan or O(n log n) per-call `searchsortedlast`.
    nk = length(k)
    d = 2
    fill!(sc, zero(dtype_sc))
    locator = 1
    for (i, tᵢ) in enumerate(t)
        if tᵢ == k[1] || tᵢ == k[end]
            # `t[1] == k[1]` and `t[end] == k[end]` by construction, so this
            # branch only fires for `i == 1` (sc[1] = 1) and `i == n`
            # (sc[end] = 1). Read directly without touching `sc`.
            on_first = tᵢ == k[1]
            diag[i] = (on_first && i == 1) || (!on_first && i == length(sc)) ?
                one(dtype_sc) : zero(dtype_sc)
            (i > 1) && (diag_lo[i - 1] = zero(dtype_sc))
            (i < n) && (diag_hi[i] = zero(dtype_sc))
            continue
        end
        # Advance the running locator until `k[locator+1] > tᵢ` — equivalent
        # to `searchsortedlast(k, tᵢ)` on monotone-increasing `tᵢ` inputs.
        while locator < nk && k[locator + 1] <= tᵢ
            locator += 1
        end
        _spline_coefficients_body!(sc, d, k, tᵢ, locator)
        diag[i] = sc[i]
        (i > 1) && (diag_lo[i - 1] = sc[i - 1])
        (i < n) && (diag_hi[i] = sc[i + 1])
        # The body writes only `sc[locator-d:locator]`; zero those entries
        # so the next iteration starts with `sc .== 0` again (the body
        # assumes positions outside its write window are already zero).
        for j in (locator - d):locator
            sc[j] = zero(dtype_sc)
        end
    end

    A = Tridiagonal(diag_lo, diag, diag_hi)

    return k, A
end

# helper function for data manipulation
function munge_data(
        u::AbstractVector, t::AbstractVector;
        check_sorted = t, sorted_arg_name = ("second", "t")
    )
    length(t) == length(u) ||
        throw(ArgumentError("`u`, `t` length mismatch: length(t) ≠ length(u)"))

    Tu = nonmissingtype(eltype(u))
    Tt = nonmissingtype(eltype(t))

    if Tu === eltype(u) && Tt === eltype(t)
        if !issorted(check_sorted; by = ForwardDiff.value)
            # there is likely an user error
            msg = "The $(sorted_arg_name[1]) argument (`$(sorted_arg_name[2])`), which is used for the interpolation domain, is not sorted."
            if issorted(u)
                msg *= "\nIt looks like the arguments `u` and `$(sorted_arg_name[2])` were inversed, make sure you used the arguments in the correct order."
            end
            throw(ArgumentError(msg))
        end

        return u, t
    end

    non_missing_mask = map((ui, ti) -> !ismissing(ui) && !ismissing(ti), u, t)
    u = convert(AbstractVector{Tu}, u[non_missing_mask])
    t = convert(AbstractVector{Tt}, t[non_missing_mask])

    return u, t
end

function munge_data(U::AbstractMatrix, t::AbstractVector)
    length(t) == size(U, 2) ||
        throw(ArgumentError("`u`, `t` length mismatch: length(t) ≠ size(U, 2)"))

    TU = nonmissingtype(eltype(U))
    Tt = nonmissingtype(eltype(t))
    if TU === eltype(U) && Tt === eltype(t)
        return U, t
    end

    non_missing_mask = map(
        (uis, ti) -> !any(ismissing, uis) && !ismissing(ti), eachcol(U), t
    )
    U = convert(AbstractMatrix{TU}, U[:, non_missing_mask])
    t = convert(AbstractVector{Tt}, t[non_missing_mask])

    return U, t
end

function munge_data(U::AbstractArray{T, N}, t) where {T, N}
    length(t) == size(U, N) ||
        throw(ArgumentError("`u`, `t` length mismatch: length(t) ≠ size(U, N)"))

    TU = nonmissingtype(eltype(U))
    Tt = nonmissingtype(eltype(t))
    if TU === eltype(U) && Tt === eltype(t)
        return U, t
    end

    non_missing_mask = map(
        (uis, ti) -> !any(ismissing, uis) && !ismissing(ti), eachslice(U; dims = N), t
    )
    U = convert(AbstractArray{TU, N}, copy(selectdim(U, N, non_missing_mask)))
    t = convert(AbstractVector{Tt}, t[non_missing_mask])

    return U, t
end

# Resolve a concrete `FindFirstFunctions.SearchStrategy` for the given
# knot vector at construction time. Stored on every interpolation cache
# as `A.strategy` so that `get_idx`'s `searchsorted_last(A.strategy, …)` is
# fully static-dispatched — no per-query `_auto_pick` branch.
#
# We dispatch to `FindFirstFunctions.Auto(t)`, which:
#
#   - Resolves a concrete `StrategyKind` from `length(t)` + the
#     `SearchProperties{T}(t)` probe at construction.
#   - For uniformly-spaced data (any `AbstractRange` or a `Vector` whose
#     every element lies within ~1e-12 of the exactly-uniform line),
#     picks `KIND_UNIFORM_STEP` and bakes the precomputed `inv_step`
#     into `props`. The hot path is then one subtract, one multiply,
#     one truncate per query — no division, no logarithmic search.
#   - For non-uniform data with `length(t) ≤ 16`, picks `KIND_LINEAR_SCAN`.
#   - Otherwise picks `KIND_BRACKET_GALLOP` (the v2 default).
#
# `Auto{T}` is parametric on the data ratio type, so the cache's
# `strategyType` parameter resolves to a single concrete `Auto{T}` per
# `t` and dispatch stays type-stable. `Vector{Int}` and `Vector{Float64}`
# both ratio-promote to `Float64`, so `Auto{Float64}` covers the common
# Float-knot cases.
@inline _resolve_strategy(t::AbstractVector) = FindFirstFunctions.Auto(t)
# Props-aware form: reuses the already-computed (possibly caller-supplied)
# `SearchProperties` instead of re-probing `t` inside `Auto(t)`, and keeps
# `A.strategy.props` consistent with `A.t_props`.
@inline _resolve_strategy(t::AbstractVector, props::FindFirstFunctions.SearchProperties) =
    FindFirstFunctions.Auto(t, props)

# Static-uniformity tag for caches. `AbstractRange{<:Real}` is uniform at the
# type level — `Val(true)` is a compile-time constant. For `AbstractVector`
# we fall through to the runtime `t_props.is_uniform` flag, which makes the
# constructor's return type a `Union{LinearInterpolation{..., true},
# LinearInterpolation{..., false}}`. Each concrete instance is fully
# type-stable per query — only the construction boundary sees the union.
@inline _static_uniform_tag(
    ::AbstractRange{<:Real}, ::FindFirstFunctions.SearchProperties
) = Val(true)
@inline _static_uniform_tag(
    ::AbstractVector, props::FindFirstFunctions.SearchProperties
) = Val(props.is_uniform)

function get_idx(
        A::AbstractInterpolation, t, iguess::Integer; lb = 1,
        ub_shift = -1, idx_shift = 0, side = :last
    )
    tvec = A.t
    ub = length(tvec) + ub_shift
    # `A.strategy` is a concrete `Auto{T}` resolved at construction time;
    # its stored kind dispatches without any per-call re-probing.
    strat = A.strategy
    raw = if side == :last
        FindFirstFunctions.searchsorted_last(strat, tvec, t, iguess)
    elseif side == :first
        FindFirstFunctions.searchsorted_first(strat, tvec, t, iguess)
    else
        error("side must be :first or :last")
    end
    return clamp(raw + idx_shift, lb, ub)
end

function get_idx(
        A::AbstractInterpolation, t, iguess::Guesser; lb = 1,
        ub_shift = -1, idx_shift = 0, side = :last
    )
    tvec = A.t
    ub = length(tvec) + ub_shift
    strat = A.strategy
    # `iguess(t)` gives a linear-extrapolation hint when `t` looks linear and
    # falls back to the cached `idx_prev` otherwise.
    hint = iguess(t)
    raw = if side == :last
        FindFirstFunctions.searchsorted_last(strat, tvec, t, hint)
    elseif side == :first
        FindFirstFunctions.searchsorted_first(strat, tvec, t, hint)
    else
        error("side must be :first or :last")
    end
    idx = clamp(raw + idx_shift, lb, ub)
    iguess.idx_prev[] = idx
    return idx
end

cumulative_integral(::AbstractInterpolation, ::Bool) = nothing
function cumulative_integral(A::AbstractInterpolation{<:Number}, cache_parameters::Bool)
    Base.require_one_based_indexing(A.u)
    idxs = cache_parameters ? (1:(length(A.t) - 1)) : (1:0)
    return cumsum(
        _integral(A, idx, t1, t2)
            for (idx, t1, t2) in
            zip(idxs, @view(A.t[begin:(end - 1)]), @view(A.t[(begin + 1):end]))
    )
end

function get_parameters(A::LinearInterpolation, idx)
    return if A.cache_parameters
        A.p.slope[idx]
    else
        linear_interpolation_parameters(A.u, A.t, idx)
    end
end

function get_parameters(A::SmoothedConstantInterpolation, idx)
    return if A.cache_parameters
        d_lower = A.p.d[idx]
        d_upper = A.p.d[idx + 1]
        c_lower = A.p.c[idx]
        c_upper = A.p.c[idx + 1]
        d_lower, d_upper, c_lower, c_upper
    else
        d_lower,
            c_lower = smoothed_constant_interpolation_parameters(
            A.u, A.t, A.d_max, idx, A.extrapolation_left, A.extrapolation_right
        )
        d_upper,
            c_upper = smoothed_constant_interpolation_parameters(
            A.u, A.t, A.d_max, idx + 1, A.extrapolation_left, A.extrapolation_right
        )
        d_lower, d_upper, c_lower, c_upper
    end
end

function get_parameters(A::QuadraticInterpolation, idx)
    return if A.cache_parameters
        A.p.α[idx], A.p.β[idx]
    else
        quadratic_interpolation_parameters(A.u, A.t, idx, A.mode)
    end
end

function get_parameters(A::QuadraticSpline, idx)
    return if A.cache_parameters
        A.p.α[idx], A.p.β[idx]
    else
        quadratic_spline_parameters(A.u, A.t, A.k, A.c, A.sc, idx)
    end
end

function get_parameters(A::CubicSpline, idx)
    return if A.cache_parameters
        A.p.c₁[idx], A.p.c₂[idx]
    else
        cubic_spline_parameters(A.u, A.h, A.z, idx)
    end
end

function get_parameters(A::CubicHermiteSpline, idx)
    return if A.cache_parameters
        A.p.c₁[idx], A.p.c₂[idx]
    else
        cubic_hermite_spline_parameters(A.du, A.u, A.t, idx)
    end
end

function get_parameters(A::QuinticHermiteSpline, idx)
    return if A.cache_parameters
        A.p.c₁[idx], A.p.c₂[idx], A.p.c₃[idx]
    else
        quintic_hermite_spline_parameters(A.ddu, A.du, A.u, A.t, idx)
    end
end

function du_PCHIP(u, t)
    h = diff(t)
    δ = diff(u) ./ h
    s = sign.(δ)

    # Special handling of the slope at the endpoints, see
    # Cleve Moler, Numerical Computing with MATLAB, Chap 3.6 (file pchiptx.m, function pchipend())
    function _edge_case(h₁, h₂, δ₁, δ₂)
        d = ((2 * h₁ + h₂) * δ₁ - h₁ * δ₂) / (h₁ + h₂)
        return if sign(d) != sign(δ₁)
            zero(eltype(δ))
        elseif sign(δ₁) != sign(δ₂) && abs(d) > 3 * abs(δ₁)
            3 * δ₁
        else
            d
        end
    end

    function _du(k)
        sₖ₋₁, sₖ = if k == 1
            s[1], s[2]
        elseif k == lastindex(t)
            s[end - 1], s[end]
        else
            s[k - 1], s[k]
        end

        return if sₖ₋₁ == 0 && sₖ == 0
            zero(eltype(δ))
        elseif sₖ₋₁ == sₖ
            if k == 1
                _edge_case(h[1], h[2], δ[1], δ[2])
            elseif k == lastindex(t)
                _edge_case(h[end], h[end - 1], δ[end], δ[end - 1])
            else
                w₁ = 2h[k] + h[k - 1]
                w₂ = h[k] + 2h[k - 1]
                (w₁ + w₂) / (w₁ / δ[k - 1] + w₂ / δ[k])
            end
        else
            if k == 1
                _edge_case(h[1], h[2], δ[1], δ[2])
            elseif k == lastindex(t)
                _edge_case(h[end], h[end - 1], δ[end], δ[end - 1])
            else
                zero(eltype(δ))
            end
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
    return Δt * (a + t_sum * (b / 2 + d * t_sq_sum / 4) + c * (t_sq_sum + t1_rel * t2_rel) / 3)
end

function integrate_quintic_polynomial(t1, t2, offset, a, b, c, d, e, f)
    t1_rel = t1 - offset
    t2_rel = t2 - offset
    t_sum = t1_rel + t2_rel
    t_sq_sum = t1_rel^2 + t2_rel^2
    t_cb_sum = t1_rel^3 + t2_rel^3
    Δt = t2 - t1
    cube_diff_factor = t_sq_sum + t1_rel * t2_rel
    return Δt * (
        a + t_sum * (b / 2 + d * t_sq_sum / 4) +
            cube_diff_factor * (c / 3 + f * t_cb_sum / 6)
    ) +
        e * (t2_rel^5 - t1_rel^5) / 5
end

function munge_extrapolation(extrapolation, extrapolation_left, extrapolation_right)
    return if extrapolation == ExtrapolationType.None
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
    return t_, n
end

function transformation_reflective(A, t)
    Δt = last(A.t) - first(A.t)
    n, t_ = fldmod(t - first(A.t), Δt)
    t_ = isodd(n) ? last(A.t) - t_ : first(A.t) + t_
    (n > 0) && (n -= 1)
    return t_, n
end

typed_nan(::AbstractArray{T}) where {T <: AbstractFloat} = T(NaN)
typed_nan(::AbstractArray{T}) where {T <: Integer} = zero(T)

# Should be replaceable by LinearAlgebra function soon: https://github.com/JuliaLang/LinearAlgebra.jl/pull/1234
function euclidean(x::AbstractArray, y::AbstractArray)
    return sqrt(mapreduce((xi, yi) -> abs2(yi - xi), +, x, y))
end

function smooth_arc_length_params_1!(Δu, u, d, j)
    uⱼ = view(u, :, j)
    uⱼ₊₁ = view(u, :, j + 1)
    dⱼ = view(d, :, j)
    dⱼ₊₁ = view(d, :, j + 1)
    @. Δu = uⱼ₊₁ - uⱼ
    d_inner = dot(dⱼ, dⱼ₊₁)
    return uⱼ, uⱼ₊₁, dⱼ, dⱼ₊₁, d_inner
end

function smooth_arc_length_params_2(u_int, uⱼ, uⱼ₊₁)
    dist₁ = euclidean(u_int, uⱼ)
    dist₂ = euclidean(u_int, uⱼ₊₁)
    Δt_line_seg = abs(dist₂ - dist₁)
    short_side_left = false
    δⱼ = if dist₁ < dist₂
        short_side_left = true
        dist₁
    else
        dist₂
    end
    return δⱼ, short_side_left, Δt_line_seg
end

function get_transition_ts(A::SmoothedConstantInterpolation)
    out = similar(A.t, 3 * length(A.t) - 2)

    for idx in 1:(length(A.t) - 1)
        d_lower, d_upper, _, _ = get_parameters(A, idx)
        if idx == 1
            out[1] = A.t[1]
            out[2] = A.t[1] + d_lower
        else
            out[3idx - 3] = A.t[idx] - d_upper
            out[3idx - 2] = A.t[idx]
            out[3idx - 1] = A.t[idx] + d_lower
        end
    end

    d_lower, _, _, _ = get_parameters(A, 1)
    out[end - 1] = A.t[end] - d_lower
    out[end] = A.t[end]

    return out
end

get_transition_ts(A::AbstractInterpolation) = A.t
