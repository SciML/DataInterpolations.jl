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

# helper function for data manipulation
function munge_data(u::AbstractVector{<:Real}, t::AbstractVector{<:Real}, safetycopy::Bool)
    if safetycopy
        u = copy(u)
        t = copy(t)
    end
    return readonly_wrap(u), readonly_wrap(t)
end

function munge_data(u::AbstractVector, t::AbstractVector, safetycopy::Bool)
    Tu = Base.nonmissingtype(eltype(u))
    Tt = Base.nonmissingtype(eltype(t))
    @assert length(t) == length(u)
    non_missing_indices = collect(
        i for i in 1:length(t)
    if !ismissing(u[i]) && !ismissing(t[i])
    )

    if safetycopy
        u = Tu.([u[i] for i in non_missing_indices])
        t = Tt.([t[i] for i in non_missing_indices])
    else
        !isempty(non_missing_indices) && throw(MustCopyError())
    end

    return readonly_wrap(u), readonly_wrap(t)
end

function munge_data(U::StridedMatrix, t::AbstractVector, safetycopy::Bool)
    TU = Base.nonmissingtype(eltype(U))
    Tt = Base.nonmissingtype(eltype(t))
    @assert length(t) == size(U, 2)
    non_missing_indices = collect(
        i for i in 1:length(t)
    if !any(ismissing, U[:, i]) && !ismissing(t[i])
    )

    if safetycopy
        U = hcat([TU.(U[:, i]) for i in non_missing_indices]...)
        t = Tt.([t[i] for i in non_missing_indices])
    else
        !isempty(non_missing_indices) && throw(MustCopyError())
    end

    return readonly_wrap(U), readonly_wrap(t)
end

seems_linear(assume_linear_t::Bool, _) = assume_linear_t
seems_linear(assume_qr::Number, t) = looks_linear(t; threshold = assume_qr)

# """
# Determine if the abscissae are regularly distributed, taking the standard \
# deviation of the difference between the array of abscissae with respect to \
# the straight line linking its first and last elements
# """
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

# Don't nest ReadOnlyArrays
readonly_wrap(a::AbstractArray) = ReadOnlyArray(a)
readonly_wrap(a::ReadOnlyArray) = a

function get_idx(tvec, t, iguess; lb = 1, ub_shift = -1, idx_shift = 0, side = :last)
    ub = length(tvec) + ub_shift
    return if side == :last
        clamp(searchsortedlastcorrelated(tvec, t, iguess) + idx_shift, lb, ub)
    elseif side == :first
        clamp(searchsortedfirstcorrelated(tvec, t, iguess) + idx_shift, lb, ub)
    else
        error("side must be :first or :last")
    end
end

function cumulative_integral(A)
    if isempty(methods(_integral, (typeof(A), Any, Any)))
        return nothing
    end
    integral_values = [_integral(A, idx, A.t[idx + 1]) - _integral(A, idx, A.t[idx])
                       for idx in 1:(length(A.t) - 1)]
    pushfirst!(integral_values, zero(first(integral_values)))
    return cumsum(integral_values)
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
