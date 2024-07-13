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

function cumulative_integral(A, ::AbstractVector{<:Number})
    integral_prototype = _integral(A, 1, A.t[2])

    integral_values = [zero(integral_prototype),
        (_integral(A, idx, A.t[idx + 1]) - _integral(A, idx, A.t[idx])
        for idx in 1:(length(A.t) - 1))...]
    return cumsum(integral_values)
end

cumulative_integral(A, ::AbstractArray) = nothing
