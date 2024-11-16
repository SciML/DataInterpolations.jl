import Base: append!, push!

function add_integral_values!(A)
    pop!(A.I)
    integral_values = Vector{eltype(A.I)}(undef, (length(A.t) - 1) - length(A.I))
    prev_sum = last(A.I)
    for i in eachindex(integral_values)
        idx = length(A.I) + i
        new_sum = prev_sum + _integral(A, idx, A.t[idx], A.t[idx + 1])
        integral_values[i] = new_sum
        prev_sum = new_sum
    end

    append!(A.I, integral_values)
end

function push!(A::LinearInterpolation{U, T}, u::eltype(U), t::eltype(T)) where {U, T}
    push!(A.u, u)
    push!(A.t, t)
    if A.cache_parameters
        slope = linear_interpolation_parameters(A.u, A.t, length(A.t) - 1)
        push!(A.p.slope, slope)
        add_integral_values!(A)
    end
    A
end

function push!(A::QuadraticInterpolation{U, T}, u::eltype(U), t::eltype(T)) where {U, T}
    # Use functionality in push! method to update parameters of last 2 intervals
    u_new = [last(A.u), u]
    t_new = [last(A.t), t]
    pop!(A.u)
    pop!(A.t)
    append!(A, u_new, t_new)
end

function push!(A::ConstantInterpolation{U, T}, u::eltype(U), t::eltype(T)) where {U, T}
    push!(A.u, u)
    push!(A.t, t)
    if A.cache_parameters
        add_integral_values!(A)
    end
    A
end

function append!(
        A::LinearInterpolation{U, T}, u::U, t::T) where {
        U, T}
    length_old = length(A.t)
    u, t = munge_data(u, t)
    append!(A.u, u)
    append!(A.t, t)
    if A.cache_parameters
        slope = linear_interpolation_parameters.(
            Ref(A.u), Ref(A.t), length_old:(length(A.t) - 1))
        append!(A.p.slope, slope)
        add_integral_values!(A)
    end
    A
end

function append!(
        A::ConstantInterpolation{U, T}, u::U, t::T) where {
        U, T}
    u, t = munge_data(u, t)
    append!(A.u, u)
    append!(A.t, t)
    if A.cache_parameters
        add_integral_values!(A)
    end
    A
end

function append!(
        A::QuadraticInterpolation{U, T}, u::U, t::T) where {
        U, T}
    u, t = munge_data(u, t)
    append!(A.u, u)
    append!(A.t, t)
    if A.cache_parameters
        pop!(A.p.α)
        pop!(A.p.β)
        parameters = quadratic_interpolation_parameters.(
            Ref(A.u), Ref(A.t), (length(A.p.α) + 1):(length(A.t) - 1), A.mode)
        α, β = collect.(eachrow(hcat(collect.(parameters)...)))
        append!(A.p.α, α)
        append!(A.p.β, β)
        add_integral_values!(A)
    end
    A
end
