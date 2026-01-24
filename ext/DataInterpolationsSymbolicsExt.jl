module DataInterpolationsSymbolicsExt

using DataInterpolations: AbstractInterpolation
import DataInterpolations: derivative
using Symbolics
using Symbolics: Num, unwrap, SymbolicUtils

@register_symbolic (interp::AbstractInterpolation)(t)
Base.nameof(interp::AbstractInterpolation) = :Interpolation

@static if pkgversion(Symbolics) >= v"7"
    @register_symbolic derivative(interp::AbstractInterpolation, t, order::Integer) false
    function SymbolicUtils.promote_symtype(
            ::typeof(derivative), Ti::SymbolicUtils.TypeT,
            Tt::SymbolicUtils.TypeT,
            To::SymbolicUtils.TypeT
        )
        @assert Ti <: AbstractInterpolation
        @assert Tt <: Real
        @assert To <: Integer
        Real
    end
    function SymbolicUtils.promote_shape(
            ::typeof(derivative),
            @nospecialize(shi::SymbolicUtils.ShapeT),
            @nospecialize(sht::SymbolicUtils.ShapeT),
            @nospecialize(sho::SymbolicUtils.ShapeT)
        )
        @assert !SymbolicUtils.is_array_shape(shi)
        @assert !SymbolicUtils.is_array_shape(sht)
        @assert !SymbolicUtils.is_array_shape(sho)
        return SymbolicUtils.ShapeVecT()
    end

    @register_derivative derivative(interp, t, ord) 2 derivative(interp, t, ord + 1)
    @register_derivative (interp::AbstractInterpolation)(t) 1 derivative(interp, t, 1)
else
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
end

end # module
