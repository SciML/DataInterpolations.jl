module DataInterpolationsSymbolicsExt

if isdefined(Base, :get_extension)
    using DataInterpolations: AbstractInterpolation
    import DataInterpolations: derivative
    using Symbolics
    using Symbolics: Num, unwrap, SymbolicUtils
else
    using ..DataInterpolations: AbstractInterpolation
    import ..DataInterpolations: derivative
    using ..Symbolics
    using ..Symbolics: Num, unwrap, SymbolicUtils
end

@register_symbolic (interp::AbstractInterpolation)(t)
Base.nameof(interp::AbstractInterpolation) = :Interpolation

function derivative(interp::AbstractInterpolation, t::Num, order = 1)
    Symbolics.wrap(SymbolicUtils.term(derivative, interp, unwrap(t), order))
end
SymbolicUtils.promote_symtype(::typeof(derivative), _...) = Real

function Symbolics.derivative(::typeof(derivative), args::NTuple{3, Any}, ::Val{2})
    Symbolics.unwrap(derivative(args[1], Symbolics.wrap(args[2]), args[3] + 1))
end

function Symbolics.derivative(interp::AbstractInterpolation, args::NTuple{1, Any}, ::Val{1})
    Symbolics.unwrap(derivative(interp, Symbolics.wrap(args[1])))
end

end # module
