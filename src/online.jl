import Base: append!, push!

function push!(di::LinearInterpolation{U, T}, u::eltype(U), t::eltype(T)) where {U, T}
    push!(di.u, u)
    push!(di.t, t)
    di
end

function push!(di::QuadraticInterpolation{U, T}, u::eltype(U), t::eltype(T)) where {U, T}
    push!(di.u, u)
    push!(di.t, t)
    di
end

function push!(di::ConstantInterpolation{U, T}, u::eltype(U), t::eltype(T)) where {U, T}
    push!(di.u, u)
    push!(di.t, t)
    di
end

function append!(di::LinearInterpolation{U, T}, u::U, t::T) where {U, T}
    u, t = munge_data(u, t)
    append!(di.u, u)
    append!(di.t, t)
    di
end

function append!(di::QuadraticInterpolation{U, T}, u::U, t::T) where {U, T}
    u, t = munge_data(u, t)
    append!(di.u, u)
    append!(di.t, t)
    di
end

function append!(di::ConstantInterpolation{U, T}, u::U, t::T) where {U, T}
    u, t = munge_data(u, t)
    append!(di.u, u)
    append!(di.t, t)
    di
end
