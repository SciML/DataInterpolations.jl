import Base: append!

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