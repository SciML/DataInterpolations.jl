module DataInterpolationsSymbolicsExt

using DataInterpolations: AbstractInterpolation
using Symbolics: Num, unwrap, SymbolicUtils
(interp::AbstractInterpolation)(t::Num) = SymbolicUtils.term(interp, unwrap(t))
SymbolicUtils.promote_symtype(t::AbstractInterpolation, _...) = Real
Base.nameof(interp::AbstractInterpolation) = :Interpolation

end # module
