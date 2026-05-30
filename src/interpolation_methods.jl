function _interpolate(A, t)
    return if t < first(A.t)
        _extrapolate_left(A, t)
    elseif t > last(A.t)
        _extrapolate_right(A, t)
    else
        _interpolate(A, t, A.iguesser)
    end
end

function _extrapolate_left(A, t)
    (; extrapolation_left) = A
    return if extrapolation_left == ExtrapolationType.None
        throw(LeftExtrapolationError())
    elseif extrapolation_left == ExtrapolationType.Constant
        slope = _derivative(A, first(A.t), 1)
        first(A.u) + zero(slope * t)
    elseif extrapolation_left == ExtrapolationType.Linear
        slope = _derivative(A, first(A.t), 1)
        first(A.u) + slope * (t - first(A.t))
    else
        _extrapolate_other(A, t, extrapolation_left)
    end
end

function _extrapolate_right(A, t)
    (; extrapolation_right) = A
    return if extrapolation_right == ExtrapolationType.None
        throw(RightExtrapolationError())
    elseif extrapolation_right == ExtrapolationType.Constant
        slope = _derivative(A, last(A.t), length(A.t))
        last(A.u) + zero(slope * t)
    elseif extrapolation_right == ExtrapolationType.Linear
        slope = _derivative(A, last(A.t), length(A.t))
        last(A.u) + slope * (t - last(A.t))
    else
        _extrapolate_other(A, t, extrapolation_right)
    end
end

function _extrapolate_other(A, t, extrapolation)
    return if extrapolation == ExtrapolationType.Extension
        _interpolate(A, t, A.iguesser)
    elseif extrapolation == ExtrapolationType.Periodic
        t_, _ = transformation_periodic(A, t)
        _interpolate(A, t_, A.iguesser)
    elseif extrapolation == ExtrapolationType.Reflective
        t_, _ = transformation_reflective(A, t)
        _interpolate(A, t_, A.iguesser)
    else
        throw(ExtrapolationNotImplementedError())
    end
end

@inline function _sorted_batch_extrapolation_ok(A)
    el = A.extrapolation_left
    er = A.extrapolation_right
    el_ok = el == ExtrapolationType.None ||
        el == ExtrapolationType.Constant ||
        el == ExtrapolationType.Linear ||
        el == ExtrapolationType.Extension
    er_ok = er == ExtrapolationType.None ||
        er == ExtrapolationType.Constant ||
        er == ExtrapolationType.Linear ||
        er == ExtrapolationType.Extension
    return el_ok && er_ok
end

@inline function _find_last_interior_index(
        tt::AbstractVector, i_first::Integer, m::Integer, tn
    )
    if i_first > m
        return i_first - 1
    end
    @inbounds if tt[m] <= tn
        return m
    end
    @inbounds if tt[i_first] > tn
        return i_first - 1
    end
    lo = i_first
    hi = m
    @inbounds while lo < hi
        mid = (lo + hi + 1) >> 1
        if tt[mid] <= tn
            lo = mid
        else
            hi = mid - 1
        end
    end
    return lo
end

@inline function _eval_interior_adaptive!(
        out::AbstractVector, t::AbstractVector, tt::AbstractVector,
        i_first::Integer, i_last::Integer, n::Integer,
        eval_at::F, strategy::FindFirstFunctions.SearchStrategy
    ) where {F}
    n_interior = i_last - i_first + 1
    n_interior <= 0 && return
    if 32 * n_interior < n
        idx_buf = Vector{Int}(undef, n_interior)
        FindFirstFunctions.searchsortedlast!(
            idx_buf, t, view(tt, i_first:i_last); strategy = strategy
        )
        @inbounds for k in 1:n_interior
            idx = idx_buf[k]
            if idx < 1
                idx = 1
            elseif idx >= n
                idx = n - 1
            end
            j = i_first + k - 1
            out[j] = eval_at(idx, tt[j])
        end
    else
        let idx = 1, i_local = i_first
            @inbounds while i_local <= i_last
                ttt = tt[i_local]
                while idx < n - 1 && ttt > t[idx + 1]
                    idx += 1
                end
                out[i_local] = eval_at(idx, ttt)
                i_local += 1
            end
        end
    end
    return
end

function _extrapolate_left(A::ConstantInterpolation, t)
    (; extrapolation_left) = A
    return if extrapolation_left == ExtrapolationType.None
        throw(LeftExtrapolationError())
    elseif extrapolation_left in (ExtrapolationType.Constant, ExtrapolationType.Linear)
        first(A.u)
    else
        _extrapolate_other(A, t, extrapolation_left)
    end
end

function _extrapolate_right(A::ConstantInterpolation, t)
    (; extrapolation_right) = A
    return if extrapolation_right == ExtrapolationType.None
        throw(RightExtrapolationError())
    elseif extrapolation_right in (ExtrapolationType.Constant, ExtrapolationType.Linear)
        last(A.u)
    else
        _extrapolate_other(A, t, extrapolation_right)
    end
end

function _extrapolate_right(A::SmoothedConstantInterpolation, t)
    return if A.extrapolation_right == ExtrapolationType.None
        throw(RightExtrapolationError())
    elseif A.extrapolation_right in (
            ExtrapolationType.Constant, ExtrapolationType.Extension,
        )
        d = min(A.t[end] - A.t[end - 1], 2A.d_max) / 2
        if A.t[end] + d < t
            A.u[end]
        else
            c = (A.u[end] - A.u[end - 1]) / 2
            A.u[end - 1] - c * (((t - A.t[end]) / d)^2 - 2 * ((t - A.t[end]) / d) - 1)
        end

    else
        _extrapolate_other(A, t, A.extrapolation_right)
    end
end

# Linear Interpolation
function _interpolate(A::LinearInterpolation{<:AbstractVector}, t::Number, iguess)
    if isnan(t)
        # For correct derivative with NaN
        idx = firstindex(A.u)
        t1 = t2 = oneunit(eltype(A.t))
        u1 = u2 = oneunit(eltype(A.u))
        slope = t / t * get_parameters(A, idx)
    else
        idx = get_idx(A, t, iguess)
        t1, t2 = A.t[idx], A.t[idx + 1]
        u1, u2 = A.u[idx], A.u[idx + 1]
        slope = get_parameters(A, idx)
    end

    Δt = t - t1
    Δu = slope * Δt
    val = u1
    Δu_nan = any(isnan.(Δu))
    if t == t2 && Δu_nan
        val = u2
    elseif !(iszero(Δt) && Δu_nan)
        val += Δu
    end
    val = oftype(Δu, val)

    return val
end

function _interpolate(A::LinearInterpolation{<:AbstractArray}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    Δt = t - A.t[idx]
    slope = get_parameters(A, idx)
    ax = axes(A.u)[1:(end - 1)]
    return A.u[ax..., idx] + slope * Δt
end

function (A::LinearInterpolation{<:AbstractVector{<:Number}})(
        out::AbstractVector, tt::AbstractVector
    )
    if length(out) != length(tt)
        throw(
            DimensionMismatch(
                "number of evaluation points and length of the result vector must be equal"
            )
        )
    end
    if _sorted_batch_extrapolation_ok(A) && issorted(tt)
        _linear_eval_sorted!(out, A, tt)
    else
        map!(A, out, tt)
    end
    return out
end

function _linear_eval_sorted!(
        out::AbstractVector, A::LinearInterpolation{<:AbstractVector{<:Number}}, tt::AbstractVector
    )
    u = A.u
    t = A.t
    el = A.extrapolation_left
    er = A.extrapolation_right
    n = length(t)
    m = length(tt)
    t1 = @inbounds t[1]
    tn = @inbounds t[n]

    i = 1

    # Left extrapolation
    if el == ExtrapolationType.None
        @inbounds if i <= m && tt[i] < t1
            throw(LeftExtrapolationError())
        end
    elseif el == ExtrapolationType.Constant
        u1 = @inbounds u[1]
        @inbounds while i <= m && tt[i] < t1
            out[i] = u1
            i += 1
        end
    else
        u1 = @inbounds u[1]
        slope1 = get_parameters(A, 1)
        @inbounds while i <= m && tt[i] < t1
            out[i] = u1 + slope1 * (tt[i] - t1)
            i += 1
        end
    end

    i_first_interior = i
    i_last_interior = _find_last_interior_index(tt, i_first_interior, m, tn)
    n_interior = i_last_interior - i_first_interior + 1
    if n_interior > 0
        if 32 * n_interior < n
            idx_buf = Vector{Int}(undef, n_interior)
            FindFirstFunctions.searchsortedlast!(
                idx_buf, t,
                view(tt, i_first_interior:i_last_interior);
                strategy = FindFirstFunctions.Auto(A.t_props)
            )
            @inbounds for k in 1:n_interior
                idx = idx_buf[k]
                if idx < 1
                    idx = 1
                elseif idx >= n
                    idx = n - 1
                end
                j = i_first_interior + k - 1
                slope = get_parameters(A, idx)
                out[j] = u[idx] + slope * (tt[j] - t[idx])
            end
        else
            let i_local = i_first_interior, idx = 1
                @inbounds while i_local <= i_last_interior
                    ttt = tt[i_local]
                    while idx < n - 1 && ttt > t[idx + 1]
                        idx += 1
                    end
                    slope = get_parameters(A, idx)
                    out[i_local] = u[idx] + slope * (ttt - t[idx])
                    i_local += 1
                end
            end
        end
    end
    i = i_last_interior + 1

    # Right extrapolation
    if er == ExtrapolationType.None
        @inbounds if i <= m
            throw(RightExtrapolationError())
        end
    elseif er == ExtrapolationType.Constant
        un = @inbounds u[n]
        @inbounds while i <= m
            out[i] = un
            i += 1
        end
    else
        un = @inbounds u[n]
        slope_n = get_parameters(A, n - 1)
        @inbounds while i <= m
            out[i] = un + slope_n * (tt[i] - tn)
            i += 1
        end
    end

    return nothing
end

# Quadratic Interpolation
function _interpolate(A::QuadraticInterpolation, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    Δt = t - A.t[idx]
    α, β = get_parameters(A, idx)
    out = A.u isa AbstractMatrix ? A.u[:, idx] : A.u[idx]
    out += @. Δt * (α * Δt + β)
    return out
end

function (A::QuadraticInterpolation{<:AbstractVector{<:Number}})(
        out::AbstractVector, tt::AbstractVector
    )
    if length(out) != length(tt)
        throw(
            DimensionMismatch(
                "number of evaluation points and length of the result vector must be equal"
            )
        )
    end
    if _sorted_batch_extrapolation_ok(A) && issorted(tt)
        _quadratic_eval_sorted!(out, A, tt)
    else
        map!(A, out, tt)
    end
    return out
end

function _quadratic_eval_sorted!(
        out::AbstractVector, A::QuadraticInterpolation{<:AbstractVector{<:Number}}, tt::AbstractVector
    )
    u = A.u
    t = A.t
    el = A.extrapolation_left
    er = A.extrapolation_right
    n = length(t)
    m = length(tt)
    t1 = @inbounds t[1]
    tn = @inbounds t[n]

    i = 1

    # Left extrapolation
    if el == ExtrapolationType.None
        @inbounds if i <= m && tt[i] < t1
            throw(LeftExtrapolationError())
        end
    elseif el == ExtrapolationType.Constant
        u1 = @inbounds u[1]
        @inbounds while i <= m && tt[i] < t1
            out[i] = u1
            i += 1
        end
    elseif el == ExtrapolationType.Linear
        u1 = @inbounds u[1]
        slope1 = _derivative(A, t1, 1)
        @inbounds while i <= m && tt[i] < t1
            out[i] = u1 + slope1 * (tt[i] - t1)
            i += 1
        end
    else
        u1 = @inbounds u[1]
        α1, β1 = get_parameters(A, 1)
        @inbounds while i <= m && tt[i] < t1
            Δt = tt[i] - t1
            out[i] = u1 + Δt * (α1 * Δt + β1)
            i += 1
        end
    end

    i_first_interior = i
    i_last_interior = _find_last_interior_index(tt, i_first_interior, m, tn)
    n_interior = i_last_interior - i_first_interior + 1
    if n_interior > 0
        if 32 * n_interior < n
            idx_buf = Vector{Int}(undef, n_interior)
            FindFirstFunctions.searchsortedlast!(
                idx_buf, t,
                view(tt, i_first_interior:i_last_interior);
                strategy = FindFirstFunctions.Auto(A.t_props)
            )
            @inbounds for k in 1:n_interior
                idx = idx_buf[k]
                if idx < 1
                    idx = 1
                elseif idx >= n
                    idx = n - 1
                end
                j = i_first_interior + k - 1
                ttt = tt[j]
                α, β = get_parameters(A, idx)
                Δt = ttt - t[idx]
                out[j] = u[idx] + Δt * (α * Δt + β)
            end
        else
            let i_local = i_first_interior, idx = 1
                @inbounds while i_local <= i_last_interior
                    ttt = tt[i_local]
                    while idx < n - 1 && ttt > t[idx + 1]
                        idx += 1
                    end
                    α, β = get_parameters(A, idx)
                    Δt = ttt - t[idx]
                    out[i_local] = u[idx] + Δt * (α * Δt + β)
                    i_local += 1
                end
            end
        end
    end
    i = i_last_interior + 1

    # Right extrapolation
    if er == ExtrapolationType.None
        @inbounds if i <= m
            throw(RightExtrapolationError())
        end
    elseif er == ExtrapolationType.Constant
        un = @inbounds u[n]
        @inbounds while i <= m
            out[i] = un
            i += 1
        end
    elseif er == ExtrapolationType.Linear
        un = @inbounds u[n]
        slope_n = _derivative(A, tn, n)
        @inbounds while i <= m
            out[i] = un + slope_n * (tt[i] - tn)
            i += 1
        end
    else
        tnm1 = @inbounds t[n - 1]
        unm1 = @inbounds u[n - 1]
        α, β = get_parameters(A, n - 1)
        @inbounds while i <= m
            Δt = tt[i] - tnm1
            out[i] = unm1 + Δt * (α * Δt + β)
            i += 1
        end
    end

    return nothing
end

# Lagrange Interpolation
function _interpolate(A::LagrangeInterpolation{<:AbstractVector}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    findRequiredIdxs!(A, t, idx)
    if A.t[A.idxs[1]] == t
        return A.u[A.idxs[1]]
    end
    N = zero(A.u[1])
    D = zero(A.t[1])
    tmp = N
    for i in 1:length(A.idxs)
        if isnan(A.bcache[A.idxs[i]])
            mult = one(A.t[1])
            for j in 1:(i - 1)
                mult *= (A.t[A.idxs[i]] - A.t[A.idxs[j]])
            end
            for j in (i + 1):length(A.idxs)
                mult *= (A.t[A.idxs[i]] - A.t[A.idxs[j]])
            end
            A.bcache[A.idxs[i]] = mult
        else
            mult = A.bcache[A.idxs[i]]
        end
        tmp = inv((t - A.t[A.idxs[i]]) * mult)
        D += tmp
        N += (tmp * A.u[A.idxs[i]])
    end
    return N / D
end

function _interpolate(A::LagrangeInterpolation{<:AbstractMatrix}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    findRequiredIdxs!(A, t, idx)
    if A.t[A.idxs[1]] == t
        return A.u[:, A.idxs[1]]
    end
    N = zero(A.u[:, 1])
    D = zero(A.t[1])
    tmp = D
    for i in 1:length(A.idxs)
        if isnan(A.bcache[A.idxs[i]])
            mult = one(A.t[1])
            for j in 1:(i - 1)
                mult *= (A.t[A.idxs[i]] - A.t[A.idxs[j]])
            end
            for j in (i + 1):length(A.idxs)
                mult *= (A.t[A.idxs[i]] - A.t[A.idxs[j]])
            end
            A.bcache[A.idxs[i]] = mult
        else
            mult = A.bcache[A.idxs[i]]
        end
        tmp = inv((t - A.t[A.idxs[i]]) * mult)
        D += tmp
        @. N += (tmp * A.u[:, A.idxs[i]])
    end
    return N / D
end

function _interpolate(A::AkimaInterpolation{<:AbstractVector}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    wj = t - A.t[idx]
    return @evalpoly wj A.u[idx] A.b[idx] A.c[idx] A.d[idx]
end

struct _AkimaEvaluator{Tu, Tt, Tb, Tc, Td}
    u::Tu
    bv::Tb
    cv::Tc
    dv::Td
    t::Tt
end
@inline function (e::_AkimaEvaluator)(idx::Integer, ttt)
    @inbounds wj = ttt - e.t[idx]
    return @inbounds @evalpoly wj e.u[idx] e.bv[idx] e.cv[idx] e.dv[idx]
end

function (A::AkimaInterpolation{<:AbstractVector})(
        out::AbstractVector, tt::AbstractVector
    )
    if length(out) != length(tt)
        throw(
            DimensionMismatch(
                "number of evaluation points and length of the result vector must be equal"
            )
        )
    end
    if _sorted_batch_extrapolation_ok(A) && issorted(tt)
        _akima_eval_sorted!(out, A, tt)
    else
        map!(A, out, tt)
    end
    return out
end

function _akima_eval_sorted!(
        out::AbstractVector, A::AkimaInterpolation{<:AbstractVector}, tt::AbstractVector
    )
    u = A.u
    t = A.t
    bv = A.b
    cv = A.c
    dv = A.d
    el = A.extrapolation_left
    er = A.extrapolation_right
    n = length(t)
    m = length(tt)
    t1 = @inbounds t[1]
    tn = @inbounds t[n]

    i = 1

    # Left extrapolation
    if el == ExtrapolationType.None
        @inbounds if i <= m && tt[i] < t1
            throw(LeftExtrapolationError())
        end
    elseif el == ExtrapolationType.Constant
        u1 = @inbounds u[1]
        @inbounds while i <= m && tt[i] < t1
            out[i] = u1
            i += 1
        end
    elseif el == ExtrapolationType.Linear
        u1 = @inbounds u[1]
        b1 = @inbounds bv[1]
        @inbounds while i <= m && tt[i] < t1
            out[i] = u1 + b1 * (tt[i] - t1)
            i += 1
        end
    else  # Extension
        u1 = @inbounds u[1]
        b1 = @inbounds bv[1]
        c1 = @inbounds cv[1]
        d1 = @inbounds dv[1]
        @inbounds while i <= m && tt[i] < t1
            wj = tt[i] - t1
            out[i] = @evalpoly wj u1 b1 c1 d1
            i += 1
        end
    end

    i_first_interior = i
    i_last_interior = _find_last_interior_index(tt, i_first_interior, m, tn)
    _eval_interior_adaptive!(
        out, t, tt, i_first_interior, i_last_interior, n,
        _AkimaEvaluator(u, bv, cv, dv, t),
        FindFirstFunctions.Auto(A.t_props)
    )
    i = i_last_interior + 1

    # Right extrapolation
    if er == ExtrapolationType.None
        @inbounds if i <= m
            throw(RightExtrapolationError())
        end
    elseif er == ExtrapolationType.Constant
        un = @inbounds u[n]
        @inbounds while i <= m
            out[i] = un
            i += 1
        end
    elseif er == ExtrapolationType.Linear
        un = @inbounds u[n]
        bn = @inbounds bv[n]
        @inbounds while i <= m
            out[i] = un + bn * (tt[i] - tn)
            i += 1
        end
    else  # Extension
        un1 = @inbounds u[n - 1]
        bn1 = @inbounds bv[n - 1]
        cn1 = @inbounds cv[n - 1]
        dn1 = @inbounds dv[n - 1]
        tn1 = @inbounds t[n - 1]
        @inbounds while i <= m
            wj = tt[i] - tn1
            out[i] = @evalpoly wj un1 bn1 cn1 dn1
            i += 1
        end
    end

    return nothing
end

# Constant Interpolation
function _interpolate(A::ConstantInterpolation{<:AbstractVector}, t::Number, iguess)
    if A.dir === :left
        # :left means that value to the left is used for interpolation
        idx = get_idx(A, t, iguess; lb = 1, ub_shift = 0)
    else
        # :right means that value to the right is used for interpolation
        idx = get_idx(A, t, iguess; side = :first, lb = 1, ub_shift = 0)
    end
    return A.u[idx]
end

function _interpolate(A::ConstantInterpolation{<:AbstractMatrix}, t::Number, iguess)
    if A.dir === :left
        # :left means that value to the left is used for interpolation
        idx = get_idx(A, t, iguess; lb = 1, ub_shift = 0)
    else
        # :right means that value to the right is used for interpolation
        idx = get_idx(A, t, iguess; side = :first, lb = 1, ub_shift = 0)
    end
    return A.u[:, idx]
end

function (A::ConstantInterpolation{<:AbstractVector})(
        out::AbstractVector, tt::AbstractVector
    )
    if length(out) != length(tt)
        throw(
            DimensionMismatch(
                "number of evaluation points and length of the result vector must be equal"
            )
        )
    end
    if _sorted_batch_extrapolation_ok(A) && issorted(tt)
        _constant_eval_sorted!(out, A, tt)
    else
        map!(A, out, tt)
    end
    return out
end

function _constant_eval_sorted!(
        out::AbstractVector, A::ConstantInterpolation{<:AbstractVector}, tt::AbstractVector
    )
    u = A.u
    t = A.t
    el = A.extrapolation_left
    er = A.extrapolation_right
    is_left = A.dir === :left
    n = length(t)
    m = length(tt)
    t1 = @inbounds t[1]
    tn = @inbounds t[n]
    u1 = @inbounds u[1]
    un = @inbounds u[n]

    i = 1

    # Left extrapolation
    if el == ExtrapolationType.None
        @inbounds if i <= m && tt[i] < t1
            throw(LeftExtrapolationError())
        end
    else
        @inbounds while i <= m && tt[i] < t1
            out[i] = u1
            i += 1
        end
    end

    i_first_interior = i
    i_last_interior = _find_last_interior_index(tt, i_first_interior, m, tn)
    n_interior = i_last_interior - i_first_interior + 1
    if n_interior > 0
        if 32 * n_interior < n
            idx_buf = Vector{Int}(undef, n_interior)
            FindFirstFunctions.searchsortedlast!(
                idx_buf, t,
                view(tt, i_first_interior:i_last_interior);
                strategy = FindFirstFunctions.Auto(A.t_props)
            )
            if is_left
                @inbounds for k in 1:n_interior
                    idx = idx_buf[k]
                    if idx < 1
                        idx = 1
                    elseif idx > n
                        idx = n
                    end
                    j = i_first_interior + k - 1
                    out[j] = u[idx]
                end
            else
                @inbounds for k in 1:n_interior
                    j = i_first_interior + k - 1
                    idx = idx_buf[k]
                    if idx < n && tt[j] != t[idx]
                        idx += 1
                    end
                    if idx < 1
                        idx = 1
                    elseif idx > n
                        idx = n
                    end
                    out[j] = u[idx]
                end
            end
        elseif is_left
            idx = 1
            @inbounds while i <= m && tt[i] <= tn
                ttt = tt[i]
                while idx < n && ttt >= t[idx + 1]
                    idx += 1
                end
                out[i] = u[idx]
                i += 1
            end
            i_last_interior = i - 1
        else
            idx = 1
            @inbounds while i <= m && tt[i] <= tn
                ttt = tt[i]
                while idx < n && ttt > t[idx]
                    idx += 1
                end
                out[i] = u[idx]
                i += 1
            end
            i_last_interior = i - 1
        end
    end
    i = i_last_interior + 1

    # Right extrapolation
    if er == ExtrapolationType.None
        @inbounds if i <= m
            throw(RightExtrapolationError())
        end
    else
        @inbounds while i <= m
            out[i] = un
            i += 1
        end
    end

    return nothing
end

# Smoothed constant Interpolation
function _interpolate(A::SmoothedConstantInterpolation{<:AbstractVector}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    d_lower, d_upper, c_lower, c_upper = get_parameters(A, idx)

    out = A.u[idx]

    if (t - A.t[idx]) < d_lower
        out -= c_lower * ((t - A.t[idx]) / d_lower - 1)^2
    elseif (A.t[idx + 1] - t) < d_upper
        out += c_upper * (1 - (A.t[idx + 1] - t) / d_upper)^2
    end

    return out
end

# QuadraticSpline Interpolation
function _interpolate(A::QuadraticSpline{<:AbstractVector}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    α, β = get_parameters(A, idx)
    uᵢ = A.u[idx]
    Δt_scaled = (t - A.t[idx]) / (A.t[idx + 1] - A.t[idx])
    return Δt_scaled * (α * Δt_scaled + β) + uᵢ
end

function (A::QuadraticSpline{<:AbstractVector{<:Number}})(
        out::AbstractVector, tt::AbstractVector
    )
    if length(out) != length(tt)
        throw(
            DimensionMismatch(
                "number of evaluation points and length of the result vector must be equal"
            )
        )
    end
    if _sorted_batch_extrapolation_ok(A) && issorted(tt)
        _quadraticspline_eval_sorted!(out, A, tt)
    else
        map!(A, out, tt)
    end
    return out
end

function _quadraticspline_eval_sorted!(
        out::AbstractVector, A::QuadraticSpline{<:AbstractVector{<:Number}}, tt::AbstractVector
    )
    u = A.u
    t = A.t
    el = A.extrapolation_left
    er = A.extrapolation_right
    n = length(t)
    m = length(tt)
    t1 = @inbounds t[1]
    tn = @inbounds t[n]

    i = 1

    # Left extrapolation
    if el == ExtrapolationType.None
        @inbounds if i <= m && tt[i] < t1
            throw(LeftExtrapolationError())
        end
    elseif el == ExtrapolationType.Constant
        u1 = @inbounds u[1]
        @inbounds while i <= m && tt[i] < t1
            out[i] = u1
            i += 1
        end
    elseif el == ExtrapolationType.Linear
        u1 = @inbounds u[1]
        slope1 = _derivative(A, t1, 1)
        @inbounds while i <= m && tt[i] < t1
            out[i] = u1 + slope1 * (tt[i] - t1)
            i += 1
        end
    else
        u1 = @inbounds u[1]
        tip1 = @inbounds t[2]
        seg_inv = inv(tip1 - t1)
        α, β = get_parameters(A, 1)
        @inbounds while i <= m && tt[i] < t1
            dts = (tt[i] - t1) * seg_inv
            out[i] = dts * (α * dts + β) + u1
            i += 1
        end
    end

    i_first_interior = i
    i_last_interior = _find_last_interior_index(tt, i_first_interior, m, tn)
    n_interior = i_last_interior - i_first_interior + 1
    if n_interior > 0
        if 32 * n_interior < n
            idx_buf = Vector{Int}(undef, n_interior)
            FindFirstFunctions.searchsortedlast!(
                idx_buf, t,
                view(tt, i_first_interior:i_last_interior);
                strategy = FindFirstFunctions.Auto(A.t_props)
            )
            @inbounds for k in 1:n_interior
                idx = idx_buf[k]
                if idx < 1
                    idx = 1
                elseif idx >= n
                    idx = n - 1
                end
                j = i_first_interior + k - 1
                ttt = tt[j]
                α, β = get_parameters(A, idx)
                seg_inv = inv(t[idx + 1] - t[idx])
                dts = (ttt - t[idx]) * seg_inv
                out[j] = dts * (α * dts + β) + u[idx]
            end
        else
            let i_local = i_first_interior, idx = 1
                @inbounds while i_local <= i_last_interior
                    ttt = tt[i_local]
                    while idx < n - 1 && ttt > t[idx + 1]
                        idx += 1
                    end
                    α, β = get_parameters(A, idx)
                    seg_inv = inv(t[idx + 1] - t[idx])
                    dts = (ttt - t[idx]) * seg_inv
                    out[i_local] = dts * (α * dts + β) + u[idx]
                    i_local += 1
                end
            end
        end
    end
    i = i_last_interior + 1

    # Right extrapolation
    if er == ExtrapolationType.None
        @inbounds if i <= m
            throw(RightExtrapolationError())
        end
    elseif er == ExtrapolationType.Constant
        un = @inbounds u[n]
        @inbounds while i <= m
            out[i] = un
            i += 1
        end
    elseif er == ExtrapolationType.Linear
        un = @inbounds u[n]
        slope_n = _derivative(A, tn, n)
        @inbounds while i <= m
            out[i] = un + slope_n * (tt[i] - tn)
            i += 1
        end
    else
        tnm1 = @inbounds t[n - 1]
        unm1 = @inbounds u[n - 1]
        seg_inv = inv(tn - tnm1)
        α, β = get_parameters(A, n - 1)
        @inbounds while i <= m
            dts = (tt[i] - tnm1) * seg_inv
            out[i] = dts * (α * dts + β) + unm1
            i += 1
        end
    end

    return nothing
end

# CubicSpline Interpolation
function _interpolate(A::CubicSpline{<:AbstractVector}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    Δt₁ = t - A.t[idx]
    Δt₂ = A.t[idx + 1] - t
    I = (A.z[idx] * Δt₂^3 + A.z[idx + 1] * Δt₁^3) / (6A.h[idx + 1])
    c₁, c₂ = get_parameters(A, idx)
    C = c₁ * Δt₁
    D = c₂ * Δt₂
    return I + C + D
end

function _interpolate(A::CubicSpline{<:AbstractArray}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    Δt₁ = t - A.t[idx]
    Δt₂ = A.t[idx + 1] - t
    ax = axes(A.z)[1:(end - 1)]
    I = (A.z[ax..., idx] * Δt₂^3 + A.z[ax..., idx + 1] * Δt₁^3) / (6A.h[idx + 1])
    c₁, c₂ = get_parameters(A, idx)
    C = c₁ * Δt₁
    D = c₂ * Δt₂
    return I + C + D
end

function (A::CubicSpline{<:AbstractVector{<:Number}})(
        out::AbstractVector, tt::AbstractVector
    )
    if length(out) != length(tt)
        throw(
            DimensionMismatch(
                "number of evaluation points and length of the result vector must be equal"
            )
        )
    end
    if _sorted_batch_extrapolation_ok(A) && issorted(tt)
        _cubicspline_eval_sorted!(out, A, tt)
    else
        map!(A, out, tt)
    end
    return out
end

@inline function _cubicspline_segment_eval(
        ttt, ti, tip1, zi, zip1, hinv6, c1, c2
    )
    dt1 = ttt - ti
    dt2 = tip1 - ttt
    return (zi * dt2 * dt2 * dt2 + zip1 * dt1 * dt1 * dt1) * hinv6 +
        c1 * dt1 + c2 * dt2
end

function _cubicspline_eval_sorted!(
        out::AbstractVector, A::CubicSpline{<:AbstractVector{<:Number}}, tt::AbstractVector
    )
    u = A.u
    t = A.t
    z = A.z
    h = A.h
    el = A.extrapolation_left
    er = A.extrapolation_right
    n = length(t)
    m = length(tt)
    t1 = @inbounds t[1]
    tn = @inbounds t[n]

    i = 1

    # Left extrapolation
    if el == ExtrapolationType.None
        @inbounds if i <= m && tt[i] < t1
            throw(LeftExtrapolationError())
        end
    elseif el == ExtrapolationType.Constant
        u1 = @inbounds u[1]
        @inbounds while i <= m && tt[i] < t1
            out[i] = u1
            i += 1
        end
    elseif el == ExtrapolationType.Linear
        u1 = @inbounds u[1]
        slope1 = _derivative(A, t1, 1)
        @inbounds while i <= m && tt[i] < t1
            out[i] = u1 + slope1 * (tt[i] - t1)
            i += 1
        end
    else  # Extension
        ti = t1
        tip1 = @inbounds t[2]
        zi = @inbounds z[1]
        zip1 = @inbounds z[2]
        hinv6 = inv(6 * @inbounds(h[2]))
        c1, c2 = get_parameters(A, 1)
        @inbounds while i <= m && tt[i] < t1
            out[i] = _cubicspline_segment_eval(
                tt[i], ti, tip1, zi, zip1, hinv6, c1, c2
            )
            i += 1
        end
    end

    i_first_interior = i
    i_last_interior = _find_last_interior_index(tt, i_first_interior, m, tn)
    n_interior = i_last_interior - i_first_interior + 1
    if n_interior > 0
        if 32 * n_interior < n
            idx_buf = Vector{Int}(undef, n_interior)
            FindFirstFunctions.searchsortedlast!(
                idx_buf, t,
                view(tt, i_first_interior:i_last_interior);
                strategy = FindFirstFunctions.Auto(A.t_props)
            )
            @inbounds for k in 1:n_interior
                idx = idx_buf[k]
                if idx < 1
                    idx = 1
                elseif idx >= n
                    idx = n - 1
                end
                j = i_first_interior + k - 1
                ttt = tt[j]
                c1, c2 = get_parameters(A, idx)
                hinv6 = inv(6 * h[idx + 1])
                out[j] = _cubicspline_segment_eval(
                    ttt, t[idx], t[idx + 1], z[idx], z[idx + 1], hinv6, c1, c2
                )
            end
        else
            let i_local = i_first_interior, idx = 1
                @inbounds while i_local <= i_last_interior
                    ttt = tt[i_local]
                    while idx < n - 1 && ttt > t[idx + 1]
                        idx += 1
                    end
                    c1, c2 = get_parameters(A, idx)
                    hinv6 = inv(6 * h[idx + 1])
                    out[i_local] = _cubicspline_segment_eval(
                        ttt, t[idx], t[idx + 1], z[idx], z[idx + 1], hinv6, c1, c2
                    )
                    i_local += 1
                end
            end
        end
    end
    i = i_last_interior + 1

    # Right extrapolation
    if er == ExtrapolationType.None
        @inbounds if i <= m
            throw(RightExtrapolationError())
        end
    elseif er == ExtrapolationType.Constant
        un = @inbounds u[n]
        @inbounds while i <= m
            out[i] = un
            i += 1
        end
    elseif er == ExtrapolationType.Linear
        un = @inbounds u[n]
        slope_n = _derivative(A, tn, n)
        @inbounds while i <= m
            out[i] = un + slope_n * (tt[i] - tn)
            i += 1
        end
    else  # Extension
        ti = @inbounds t[n - 1]
        tip1 = tn
        zi = @inbounds z[n - 1]
        zip1 = @inbounds z[n]
        hinv6 = inv(6 * @inbounds(h[n]))
        c1, c2 = get_parameters(A, n - 1)
        @inbounds while i <= m
            out[i] = _cubicspline_segment_eval(
                tt[i], ti, tip1, zi, zip1, hinv6, c1, c2
            )
            i += 1
        end
    end

    return nothing
end

# BSpline Curve Interpolation
function _interpolate(
        A::BSplineInterpolation{<:AbstractVector{<:Number}},
        t::Number,
        iguess
    )
    t < A.t[1] && return A.u[1]
    t > A.t[end] && return A.u[end]
    # change t into param [0 1]
    idx = get_idx(A, t, iguess)
    t = A.p[idx] + (t - A.t[idx]) / (A.t[idx + 1] - A.t[idx]) * (A.p[idx + 1] - A.p[idx])
    n = length(A.t)
    # Per-call scratch buffer: evaluation must be reentrant for thread safety (#532)
    sc = zeros(eltype(t), n)
    nonzero_coefficient_idxs = spline_coefficients!(sc, A.d, A.k, t)
    ucum = zero(eltype(A.u))
    for i in nonzero_coefficient_idxs
        ucum += sc[i] * A.c[i]
    end
    return ucum
end

function _interpolate(
        A::BSplineInterpolation{<:AbstractArray{<:Number}},
        t::Number,
        iguess
    )
    ax_u = axes(A.u)[1:(end - 1)]
    t < A.t[1] && return A.u[ax_u..., 1]
    t > A.t[end] && return A.u[ax_u..., end]
    # change t into param [0 1]
    idx = get_idx(A, t, iguess)
    t = A.p[idx] + (t - A.t[idx]) / (A.t[idx + 1] - A.t[idx]) * (A.p[idx + 1] - A.p[idx])
    n = length(A.t)
    # Per-call scratch buffer: evaluation must be reentrant for thread safety (#532)
    sc = zeros(eltype(t), n)
    nonzero_coefficient_idxs = spline_coefficients!(sc, A.d, A.k, t)
    ucum = zeros(eltype(A.u), size(A.u)[1:(end - 1)]...)
    for i in nonzero_coefficient_idxs
        ucum = ucum + (sc[i] * A.c[ax_u..., i])
    end
    return ucum
end

# BSpline Curve Approx
function _interpolate(A::BSplineApprox{<:AbstractVector{<:Number}}, t::Number, iguess)
    t < A.t[1] && return A.u[1]
    t > A.t[end] && return A.u[end]
    # change t into param [0 1]
    idx = get_idx(A, t, iguess)
    t = A.p[idx] + (t - A.t[idx]) / (A.t[idx + 1] - A.t[idx]) * (A.p[idx + 1] - A.p[idx])
    # Per-call scratch buffer: evaluation must be reentrant for thread safety (#532)
    sc = zeros(eltype(t), A.h)
    nonzero_coefficient_idxs = spline_coefficients!(sc, A.d, A.k, t)
    ucum = zero(eltype(A.u))
    for i in nonzero_coefficient_idxs
        ucum += sc[i] * A.c[i]
    end
    return ucum
end

function _interpolate(
        A::BSplineApprox{<:AbstractArray{<:Number}}, t::Number, iguess
    )
    ax_u = axes(A.u)[1:(end - 1)]
    t < A.t[1] && return A.u[ax_u..., 1]
    t > A.t[end] && return A.u[ax_u..., end]
    # change t into param [0 1]
    idx = get_idx(A, t, iguess)
    t = A.p[idx] + (t - A.t[idx]) / (A.t[idx + 1] - A.t[idx]) * (A.p[idx + 1] - A.p[idx])
    # Per-call scratch buffer: evaluation must be reentrant for thread safety (#532)
    sc = zeros(eltype(t), A.h)
    nonzero_coefficient_idxs = spline_coefficients!(sc, A.d, A.k, t)
    ucum = zeros(eltype(A.u), size(A.u)[1:(end - 1)]...)
    for i in nonzero_coefficient_idxs
        ucum = ucum + (sc[i] * A.c[ax_u..., i])
    end
    return ucum
end

# Cubic Hermite Spline
function _interpolate(
        A::CubicHermiteSpline{<:AbstractVector}, t::Number, iguess
    )
    idx = get_idx(A, t, iguess)
    Δt₀ = t - A.t[idx]
    Δt₁ = t - A.t[idx + 1]
    out = A.u[idx] + Δt₀ * A.du[idx]
    c₁, c₂ = get_parameters(A, idx)
    out += Δt₀^2 * (c₁ + Δt₁ * c₂)
    return out
end

@inline function _cubic_hermite_segment_eval(ttt, ti, tip1, ui, dui, c1, c2)
    dt0 = ttt - ti
    dt1 = ttt - tip1
    return ui + dt0 * dui + dt0 * dt0 * (c1 + dt1 * c2)
end

function (A::CubicHermiteSpline{<:AbstractVector{<:Number}})(
        out::AbstractVector, tt::AbstractVector
    )
    if length(out) != length(tt)
        throw(
            DimensionMismatch(
                "number of evaluation points and length of the result vector must be equal"
            )
        )
    end
    if _sorted_batch_extrapolation_ok(A) && issorted(tt)
        _cubic_hermite_eval_sorted!(out, A, tt)
    else
        map!(A, out, tt)
    end
    return out
end

function _cubic_hermite_eval_sorted!(
        out::AbstractVector, A::CubicHermiteSpline{<:AbstractVector{<:Number}}, tt::AbstractVector
    )
    u = A.u
    t = A.t
    du = A.du
    el = A.extrapolation_left
    er = A.extrapolation_right
    n = length(t)
    m = length(tt)
    t1 = @inbounds t[1]
    tn = @inbounds t[n]

    i = 1

    # Left extrapolation
    if el == ExtrapolationType.None
        @inbounds if i <= m && tt[i] < t1
            throw(LeftExtrapolationError())
        end
    elseif el == ExtrapolationType.Constant
        u1 = @inbounds u[1]
        @inbounds while i <= m && tt[i] < t1
            out[i] = u1
            i += 1
        end
    elseif el == ExtrapolationType.Linear
        u1 = @inbounds u[1]
        du1 = @inbounds du[1]
        @inbounds while i <= m && tt[i] < t1
            out[i] = u1 + du1 * (tt[i] - t1)
            i += 1
        end
    else  # Extension
        ui = @inbounds u[1]
        dui = @inbounds du[1]
        tip1 = @inbounds t[2]
        c1, c2 = get_parameters(A, 1)
        @inbounds while i <= m && tt[i] < t1
            out[i] = _cubic_hermite_segment_eval(tt[i], t1, tip1, ui, dui, c1, c2)
            i += 1
        end
    end

    i_first_interior = i
    i_last_interior = _find_last_interior_index(tt, i_first_interior, m, tn)
    n_interior = i_last_interior - i_first_interior + 1
    if n_interior > 0
        if 32 * n_interior < n
            idx_buf = Vector{Int}(undef, n_interior)
            FindFirstFunctions.searchsortedlast!(
                idx_buf, t,
                view(tt, i_first_interior:i_last_interior);
                strategy = FindFirstFunctions.Auto(A.t_props)
            )
            @inbounds for k in 1:n_interior
                idx = idx_buf[k]
                if idx < 1
                    idx = 1
                elseif idx >= n
                    idx = n - 1
                end
                j = i_first_interior + k - 1
                ttt = tt[j]
                c1, c2 = get_parameters(A, idx)
                out[j] = _cubic_hermite_segment_eval(
                    ttt, t[idx], t[idx + 1], u[idx], du[idx], c1, c2
                )
            end
        else
            let i_local = i_first_interior, idx = 1
                @inbounds while i_local <= i_last_interior
                    ttt = tt[i_local]
                    while idx < n - 1 && ttt > t[idx + 1]
                        idx += 1
                    end
                    c1, c2 = get_parameters(A, idx)
                    out[i_local] = _cubic_hermite_segment_eval(
                        ttt, t[idx], t[idx + 1], u[idx], du[idx], c1, c2
                    )
                    i_local += 1
                end
            end
        end
    end
    i = i_last_interior + 1

    # Right extrapolation
    if er == ExtrapolationType.None
        @inbounds if i <= m
            throw(RightExtrapolationError())
        end
    elseif er == ExtrapolationType.Constant
        un = @inbounds u[n]
        @inbounds while i <= m
            out[i] = un
            i += 1
        end
    elseif er == ExtrapolationType.Linear
        un = @inbounds u[n]
        dun = @inbounds du[n]
        @inbounds while i <= m
            out[i] = un + dun * (tt[i] - tn)
            i += 1
        end
    else  # Extension
        ti = @inbounds t[n - 1]
        ui = @inbounds u[n - 1]
        dui = @inbounds du[n - 1]
        c1, c2 = get_parameters(A, n - 1)
        @inbounds while i <= m
            out[i] = _cubic_hermite_segment_eval(tt[i], ti, tn, ui, dui, c1, c2)
            i += 1
        end
    end

    return nothing
end

# Quintic Hermite Spline
function _interpolate(
        A::QuinticHermiteSpline{<:AbstractVector}, t::Number, iguess
    )
    idx = get_idx(A, t, iguess)
    Δt₀ = t - A.t[idx]
    Δt₁ = t - A.t[idx + 1]
    out = A.u[idx] + Δt₀ * (A.du[idx] + A.ddu[idx] * Δt₀ / 2)
    c₁, c₂, c₃ = get_parameters(A, idx)
    out += Δt₀^3 * (c₁ + Δt₁ * (c₂ + c₃ * Δt₁))
    return out
end

@inline function _quintic_hermite_segment_eval(
        ttt, ti, tip1, ui, dui, ddui, c1, c2, c3
    )
    dt0 = ttt - ti
    dt1 = ttt - tip1
    quad = ui + dt0 * (dui + ddui * dt0 / 2)
    high = dt0 * dt0 * dt0 * (c1 + dt1 * (c2 + c3 * dt1))
    return quad + high
end

function (A::QuinticHermiteSpline{<:AbstractVector{<:Number}})(
        out::AbstractVector, tt::AbstractVector
    )
    if length(out) != length(tt)
        throw(
            DimensionMismatch(
                "number of evaluation points and length of the result vector must be equal"
            )
        )
    end
    if _sorted_batch_extrapolation_ok(A) && issorted(tt)
        _quintic_hermite_eval_sorted!(out, A, tt)
    else
        map!(A, out, tt)
    end
    return out
end

function _quintic_hermite_eval_sorted!(
        out::AbstractVector, A::QuinticHermiteSpline{<:AbstractVector{<:Number}}, tt::AbstractVector
    )
    u = A.u
    t = A.t
    du = A.du
    ddu = A.ddu
    el = A.extrapolation_left
    er = A.extrapolation_right
    n = length(t)
    m = length(tt)
    t1 = @inbounds t[1]
    tn = @inbounds t[n]

    i = 1

    if el == ExtrapolationType.None
        @inbounds if i <= m && tt[i] < t1
            throw(LeftExtrapolationError())
        end
    elseif el == ExtrapolationType.Constant
        u1 = @inbounds u[1]
        @inbounds while i <= m && tt[i] < t1
            out[i] = u1
            i += 1
        end
    elseif el == ExtrapolationType.Linear
        u1 = @inbounds u[1]
        du1 = @inbounds du[1]
        @inbounds while i <= m && tt[i] < t1
            out[i] = u1 + du1 * (tt[i] - t1)
            i += 1
        end
    else  # Extension
        ui = @inbounds u[1]
        dui = @inbounds du[1]
        ddui = @inbounds ddu[1]
        tip1 = @inbounds t[2]
        c1, c2, c3 = get_parameters(A, 1)
        @inbounds while i <= m && tt[i] < t1
            out[i] = _quintic_hermite_segment_eval(
                tt[i], t1, tip1, ui, dui, ddui, c1, c2, c3
            )
            i += 1
        end
    end

    i_first_interior = i
    i_last_interior = _find_last_interior_index(tt, i_first_interior, m, tn)
    n_interior = i_last_interior - i_first_interior + 1
    if n_interior > 0
        if 32 * n_interior < n
            idx_buf = Vector{Int}(undef, n_interior)
            FindFirstFunctions.searchsortedlast!(
                idx_buf, t,
                view(tt, i_first_interior:i_last_interior);
                strategy = FindFirstFunctions.Auto(A.t_props)
            )
            @inbounds for k in 1:n_interior
                idx = idx_buf[k]
                if idx < 1
                    idx = 1
                elseif idx >= n
                    idx = n - 1
                end
                j = i_first_interior + k - 1
                ttt = tt[j]
                c1, c2, c3 = get_parameters(A, idx)
                out[j] = _quintic_hermite_segment_eval(
                    ttt, t[idx], t[idx + 1], u[idx], du[idx], ddu[idx], c1, c2, c3
                )
            end
        else
            let i_local = i_first_interior, idx = 1
                @inbounds while i_local <= i_last_interior
                    ttt = tt[i_local]
                    while idx < n - 1 && ttt > t[idx + 1]
                        idx += 1
                    end
                    c1, c2, c3 = get_parameters(A, idx)
                    out[i_local] = _quintic_hermite_segment_eval(
                        ttt, t[idx], t[idx + 1], u[idx], du[idx], ddu[idx], c1, c2, c3
                    )
                    i_local += 1
                end
            end
        end
    end
    i = i_last_interior + 1

    if er == ExtrapolationType.None
        @inbounds if i <= m
            throw(RightExtrapolationError())
        end
    elseif er == ExtrapolationType.Constant
        un = @inbounds u[n]
        @inbounds while i <= m
            out[i] = un
            i += 1
        end
    elseif er == ExtrapolationType.Linear
        un = @inbounds u[n]
        dun = @inbounds du[n]
        @inbounds while i <= m
            out[i] = un + dun * (tt[i] - tn)
            i += 1
        end
    else  # Extension
        ti = @inbounds t[n - 1]
        ui = @inbounds u[n - 1]
        dui = @inbounds du[n - 1]
        ddui = @inbounds ddu[n - 1]
        c1, c2, c3 = get_parameters(A, n - 1)
        @inbounds while i <= m
            out[i] = _quintic_hermite_segment_eval(
                tt[i], ti, tn, ui, dui, ddui, c1, c2, c3
            )
            i += 1
        end
    end

    return nothing
end

function _interpolate(A::SmoothArcLengthInterpolation, t::Number, iguess)
    (; out, in_place) = A
    out = in_place ? out : similar(out, typeof(t))
    idx = get_idx(A, t, iguess)
    Δt_circ_seg = A.Δt_circle_segment[idx]
    Δt_line_seg = A.Δt_line_segment[idx]
    short_side_left = A.short_side_left[idx]
    Δt = t - A.t[idx]

    in_circle_arc = if short_side_left
        Δt < Δt_circ_seg
    else
        Δt > Δt_line_seg
    end

    if in_circle_arc
        t_circle_seg = short_side_left ? Δt : Δt - Δt_line_seg
        Rⱼ = A.radius[idx]
        S, C = sincos(t_circle_seg / Rⱼ)
        c = view(A.center, :, idx)
        v₁ = view(A.dir_1, :, idx)
        v₂ = view(A.dir_2, :, idx)
        @. out = c + C * v₁ + S * v₂
    else
        if short_side_left
            u₁ = view(A.u, :, idx + 1)
            d₁ = view(A.d, :, idx + 1)
            t_line_seg = A.t[idx + 1] - t
            @. out = u₁ - t_line_seg * d₁
        else
            u₀ = view(A.u, :, idx)
            d₀ = view(A.d, :, idx)
            @. out = u₀ + Δt * d₀
        end
    end

    return out
end
