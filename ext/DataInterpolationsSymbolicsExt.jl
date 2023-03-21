module DataInterpolationsSymbolicsExt

if isdefined(Base, :get_extension)
    using DataInterpolations: AbstractInterpolation
    using Symbolics: Num, unwrap, SymbolicUtils
else
    using ..DataInterpolations: AbstractInterpolation
    using ..Symbolics: Num, unwrap, SymbolicUtils
end

(interp::AbstractInterpolation)(t::Num) = SymbolicUtils.term(interp, unwrap(t))
SymbolicUtils.promote_symtype(t::AbstractInterpolation, _...) = Real
Base.nameof(interp::AbstractInterpolation) = :Interpolation

end # module
