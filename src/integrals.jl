function integral(A::AbstractInterpolation, t::Number)
    bw, fw = samples(A)
    idx = max(1 + bw, min(searchsortedlast(A.t, t), length(A.t) - fw))
    _integral(A, idx, t)
end

function integral(A::AbstractInterpolation, t1::Number, t2::Number)
    bw, fw = samples(A)
    # the index less than or equal to t1
    idx1 = max(1 + bw, min(searchsortedlast(A.t, t1), length(A.t) - fw))
    # the index less than t2
    idx2 = max(2 + bw, min(searchsortedlast(A.t, t2), length(A.t) - fw))
    if A.t[idx2] == t2
        idx2 -= 1
    end
    total = zero(eltype(A))
    for idx in idx1:idx2
        lt1 = idx == idx1 ? t1 : A.t[idx]
        lt2 = idx == idx2 ? t2 : A.t[idx + 1]
        total += _integral(A, idx, lt2) - _integral(A, idx, lt1)
    end
    total
end

samples(A::LinearInterpolation{<:AbstractVector}) = (0, 1)
function _integral(A::LinearInterpolation{<:AbstractVector{<:Number}},
    idx::Number,
    t::Number)
    t1 = A.t[idx]
    t2 = A.t[idx + 1]
    u1 = A.u[idx]
    u2 = A.u[idx + 1]
    t^2 * (u1 - u2) / (2 * t1 - 2 * t2) + t * (t1 * u2 - t2 * u1) / (t1 - t2)
end

samples(A::ConstantInterpolation{<:AbstractVector}) = (0, 1)
function _integral(A::ConstantInterpolation{<:AbstractVector}, idx::Number, t::Number)
    if A.dir === :left
        # :left means that value to the left is used for interpolation
        return A.u[idx] * t
    else
        # :right means that value to the right is used for interpolation
        return A.u[idx + 1] * t
    end
end

samples(A::QuadraticInterpolation{<:AbstractVector}) = (0, 1)
function _integral(A::QuadraticInterpolation{<:AbstractVector{<:Number}},
    idx::Number,
    t::Number)
    A.mode == :Backward && idx > 1 && (idx -= 1)
    idx = min(length(A.t) - 2, idx)
    t1 = A.t[idx]
    t2 = A.t[idx + 1]
    t3 = A.t[idx + 2]
    u1 = A.u[idx]
    u2 = A.u[idx + 1]
    u3 = A.u[idx + 2]
    (t^3 * (-t1 * u2 + t1 * u3 + t2 * u1 - t2 * u3 - t3 * u1 + t3 * u2) /
     (3 * t1^2 * t2 - 3 * t1^2 * t3 - 3 * t1 * t2^2 + 3 * t1 * t3^2 + 3 * t2^2 * t3 -
      3 * t2 * t3^2) +
     t^2 * (t1^2 * u2 - t1^2 * u3 - t2^2 * u1 + t2^2 * u3 + t3^2 * u1 - t3^2 * u2) /
     (2 * t1^2 * t2 - 2 * t1^2 * t3 - 2 * t1 * t2^2 + 2 * t1 * t3^2 + 2 * t2^2 * t3 -
      2 * t2 * t3^2) +
     t *
     (t1^2 * t2 * u3 - t1^2 * t3 * u2 - t1 * t2^2 * u3 + t1 * t3^2 * u2 + t2^2 * t3 * u1 -
      t2 * t3^2 * u1) /
     (t1^2 * t2 - t1^2 * t3 - t1 * t2^2 + t1 * t3^2 + t2^2 * t3 - t2 * t3^2))
end

samples(A::QuadraticSpline{<:AbstractVector{<:Number}}) = (0, 1)
function _integral(A::QuadraticSpline{<:AbstractVector{<:Number}}, idx::Number, t::Number)
    t1 = A.t[idx]
    t2 = A.t[idx + 1]
    u1 = A.u[idx]
    z1 = A.z[idx]
    z2 = A.z[idx + 1]
    t^3 * (z1 - z2) / (6 * t1 - 6 * t2) + t^2 * (t1 * z2 - t2 * z1) / (2 * t1 - 2 * t2) +
    t * (-t1^2 * z1 - t1^2 * z2 + 2 * t1 * t2 * z1 + 2 * t1 * u1 - 2 * t2 * u1) /
    (2 * t1 - 2 * t2)
end

samples(A::CubicSpline{<:AbstractVector{<:Number}}) = (0, 1)
function _integral(A::CubicSpline{<:AbstractVector{<:Number}}, idx::Number, t::Number)
    t1 = A.t[idx]
    t2 = A.t[idx + 1]
    u1 = A.u[idx]
    u2 = A.u[idx + 1]
    z1 = A.z[idx]
    z2 = A.z[idx + 1]
    h2 = A.h[idx + 1]
    (t^4 * (-z1 + z2) / (24 * h2) + t^3 * (-t1 * z2 + t2 * z1) / (6 * h2) +
     t^2 * (h2^2 * z1 - h2^2 * z2 + 3 * t1^2 * z2 - 3 * t2^2 * z1 - 6 * u1 + 6 * u2) /
     (12 * h2) +
     t *
     (h2^2 * t1 * z2 - h2^2 * t2 * z1 - t1^3 * z2 - 6 * t1 * u2 + t2^3 * z1 + 6 * t2 * u1) /
     (6 * h2))
end
