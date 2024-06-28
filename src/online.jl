import Base: append!, push!

function push!(A::LinearInterpolation{U, T}, u::eltype(U), t::eltype(T)) where {U, T}
    push!(A.u.parent, u)
    push!(A.t.parent, t)
    slope = LinearInterpolationParameters(A.u, A.t, size(A.t)[1] - 1)
    push!(A.p.slope, slope)
    A
end

function push!(di::QuadraticInterpolation{U, T}, u::eltype(U), t::eltype(T)) where {U, T}
    push!(di.u, u)
    push!(di.t, t)
    di
end

function push!(A::ConstantInterpolation{U, T}, u::eltype(U), t::eltype(T)) where {U, T}
    push!(A.u.parent, u)
    push!(A.t.parent, t)
    A
end

function append!(
        A::LinearInterpolation{ReadOnlyVector{eltypeU, U}, ReadOnlyVector{eltypeT, T}}, u::U, t::T) where {
        eltypeU, U, eltypeT, T}
    length_old = size(A.t)[1]
    u, t = munge_data(u, t, true)
    append!(A.u.parent, u)
    append!(A.t.parent, t)
    slope = LinearInterpolationParameters.(
        Ref(A.u), Ref(A.t), length_old:(size(A.t)[1] - 1))
    append!(A.p.slope, slope)
    A
end

function append!(
    A::ConstantInterpolation{ReadOnlyVector{eltypeU, U}, ReadOnlyVector{eltypeT, T}}, u::U, t::T) where {
        eltypeU, U, eltypeT, T}
    u, t = munge_data(u, t, true)
    append!(A.parent, u)
    append!(A.parent.t, t)
    A
end

function append!(
        di::ConstantInterpolation{U, T}, u::U, t::T, safetycopy::Bool = false) where {U, T}
    u, t = munge_data(u, t, safetycopy)
    append!(di.u, u)
    append!(di.t, t)
    di
end
