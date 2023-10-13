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

(interp::AbstractInterpolation)(t::Num) = SymbolicUtils.term(interp, unwrap(t))
SymbolicUtils.promote_symtype(t::AbstractInterpolation, _...) = Real
Base.nameof(interp::AbstractInterpolation) = :Interpolation

function derivative(interp::AbstractInterpolation, t::Num)
    SymbolicUtils.term(derivative, interp, unwrap(t))
end
SymbolicUtils.promote_symtype(::typeof(derivative), _...) = Real

function Symbolics.derivative(interp::AbstractInterpolation, args::NTuple{1, Any}, ::Val{1})
    Symbolics.unwrap(derivative(interp, Symbolics.wrap(args[1])))
end

end # module
