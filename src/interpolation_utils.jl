function findRequiredIdxs(A::LagrangeInterpolation, t)
    idxs = sortperm(A.t, by = x -> abs(t - x))
    idxs[1:(A.n + 1)]
end

function spline_coefficients(n, d, k, u::Number)
    N = zeros(n)
    if u == k[1]
        N[1] = one(u)
    elseif u == k[end]
        N[end] = one(u)
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
    end
    N
end

function spline_coefficients(n, d, k, u::AbstractVector)
    N = zeros(eltype(u), n, n)
    for i in 1:n
        N[i, :] .= spline_coefficients(n, d, k, u[i])
    end
    N
end

# helper function for data manipulation
function munge_data(u::AbstractVector{<:Real}, t::AbstractVector{<:Real})
    return u, t
end

function munge_data(u::AbstractVector, t::AbstractVector)
    Tu = Base.nonmissingtype(eltype(u))
    Tt = Base.nonmissingtype(eltype(t))
    @assert length(t) == length(u)
    non_missing_indices = collect(i for i in 1:length(t)
                                      if !ismissing(u[i]) && !ismissing(t[i]))
    newu = Tu.([u[i] for i in non_missing_indices])
    newt = Tt.([t[i] for i in non_missing_indices])

    return newu, newt
end

function munge_data(U::StridedMatrix, t::AbstractVector)
    TU = Base.nonmissingtype(eltype(U))
    Tt = Base.nonmissingtype(eltype(t))
    @assert length(t) == size(U, 2)
    non_missing_indices = collect(i for i in 1:length(t)
                                      if !any(ismissing, U[:, i]) && !ismissing(t[i]))
    newUs = [TU.(U[:, i]) for i in non_missing_indices]
    newt = Tt.([t[i] for i in non_missing_indices])

    return hcat(newUs...), newt
end

"""
    bracketstrictlymontonic(v, x, guess; lt=<comparison>, by=<transform>, rev=false)

Starting from an initial `guess` index, find indices `(lo, hi)` such that `v[lo] ≤ x ≤
v[hi]` according to the specified order, assuming that `x` is actually within the range of
values found in `v`.  If `x` is outside that range, either `lo` will be `firstindex(v)` or
`hi` will be `lastindex(v)`.

Note that the results will not typically satisfy `lo ≤ guess ≤ hi`.  If `x` is precisely
equal to a value that is not unique in the input `v`, there is no guarantee that `(lo, hi)`
will encompass *all* indices corresponding to that value.

This algorithm is essentially an expanding binary search, which can be used as a precursor
to `searchsorted` and related functions, which can take `lo` and `hi` as arguments.  The
purpose of using this function first would be to accelerate convergence in those functions
by using correlated `guess`es for repeated calls.  The best `guess` for the next call of
this function would be the index returned by the previous call to `searchsorted`.

See [`sort!`](@ref) for an explanation of the keyword arguments `by`, `lt` and `rev`.
"""
function bracketstrictlymontonic(v::AbstractVector,
    x,
    guess::T,
    o::Base.Order.Ordering)::NTuple{2, keytype(v)} where {T <: Integer}
    bottom = firstindex(v)
    top = lastindex(v)
    if guess < bottom || guess > top
        return bottom, top
        # # NOTE: for cache efficiency in repeated calls, we avoid accessing the first and last elements of `v`
        # # on each call to this function.  This should only result in significant slow downs for calls with
        # # out-of-bounds values of `x` *and* bad `guess`es.
        # elseif lt(o, x, v[bottom])
        #     return bottom, bottom
        # elseif lt(o, v[top], x)
        #     return top, top
    else
        u = T(1)
        lo, hi = guess, min(guess + u, top)
        @inbounds if Base.Order.lt(o, x, v[lo])
            while lo > bottom && Base.Order.lt(o, x, v[lo])
                lo, hi = max(bottom, lo - u), lo
                u += u
            end
        else
            while hi < top && !Base.Order.lt(o, x, v[hi])
                lo, hi = hi, min(top, hi + u)
                u += u
            end
        end
    end
    return lo, hi
end

function searchsortedfirstcorrelated(v::AbstractVector, x, guess)
    lo, hi = bracketstrictlymontonic(v, x, guess, Base.Order.Forward)
    searchsortedfirst(v, x, lo, hi, Base.Order.Forward)
end

function searchsortedlastcorrelated(v::AbstractVector, x, guess)
    lo, hi = bracketstrictlymontonic(v, x, guess, Base.Order.Forward)
    searchsortedlast(v, x, lo, hi, Base.Order.Forward)
end

searchsortedfirstcorrelated(r::AbstractRange, x, _) = searchsortedfirst(r, x)
searchsortedlastcorrelated(r::AbstractRange, x, _) = searchsortedlast(r, x)
