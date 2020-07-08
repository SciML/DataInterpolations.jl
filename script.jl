function munge_data_mine(u::AbstractVector{<:Real}, t::AbstractVector{<:Real})
  return u, t
end

function munge_data_mine(u::AbstractVector{<:Union{Real,Missing}}, t::AbstractVector{<:Union{Real,Missing}})
  Tu = Base.nonmissingtype(eltype(u))
  Tt = Base.nonmissingtype(eltype(t))
  @assert length(t) == length(u)
  non_missing_indices = collect(i for i in 1:length(t) if !ismissing(u[i]) && !ismissing(t[i]))
  newu = Tu.([u[i] for i in non_missing_indices])
  newt = Tt.([t[i] for i in non_missing_indices])

  return newu, newt
end

function munge_data_mine(U::StridedMatrix, t::AbstractVector)
  TU = Base.nonmissingtype(eltype(U))
  Tt = Base.nonmissingtype(eltype(t))
  @assert length(t) == size(U,2)
  non_missing_indices = collect(i for i in 1:length(t) if !any(ismissing,U[:,i]) && !ismissing(t[i]))
  newUs = [TU.(U[:,i]) for i in non_missing_indices]
  newt= Tt.([t[i] for i in non_missing_indices])

  return hcat(newUs...), newt
end

A = [1;1]
B = [2;2]
f = x -> munge_data_mine(x, A)[1][1]
Zygote.gradient(f, [1;1])
