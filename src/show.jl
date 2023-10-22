###################### Generic Dispatches ######################

function Base.show(io::IO, mime::MIME"text/plain", interp::AbstractInterpolation)
    print(io, get_show(interp))
    header = ["time", get_names(interp.u)...]
    data = hcat(interp.t, get_data(interp.u))
    pretty_table(io, data; header = header, vcrop_mode = :middle)
end

function get_show(interp::AbstractInterpolation)
    return string(nameof(typeof(interp))) * " with $(length(interp.t)) points\n"
end

function get_data(u::AbstractVector)
    return u
end

function get_data(u::AbstractVector{<:AbstractVector})
    return reduce(hcat, u)'
end

function get_data(u::AbstractMatrix)
    return u'
end

function get_names(u::AbstractVector)
    return ["u"]
end

function get_names(u::AbstractVector{<:AbstractVector})
    return ["u$i" for i in eachindex(first(u))]
end

function get_names(u::AbstractMatrix)
    return ["u$i" for i in axes(u, 1)]
end

###################### Specific Dispatches ######################

function get_show(interp::QuadraticInterpolation)
    return string(nameof(typeof(interp))) *
           " with $(length(interp.t)) points, $(interp.mode) mode\n"
end

function get_show(interp::LagrangeInterpolation)
    return string(nameof(typeof(interp))) *
           " with $(length(interp.t)) points, with order $(interp.n)\n"
end

function get_show(interp::ConstantInterpolation)
    return string(nameof(typeof(interp))) *
           " with $(length(interp.t)) points, in $(interp.dir) direction\n"
end

function get_show(interp::BSplineInterpolation)
    return string(nameof(typeof(interp))) *
           " with $(length(interp.t)) points, with degree $(interp.d)\n"
end

function get_show(interp::BSplineApprox)
    return string(nameof(typeof(interp))) *
           " with $(length(interp.t)) points, with degree $(interp.d), number of control points $(interp.h)\n"
end
