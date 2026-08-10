using SciMLTesting, DataInterpolations, Test

# ExplicitImports only sees an extension once the extension module exists, and an
# extension module exists only once its trigger packages are loaded. Loading every
# weakdep here is what makes the checks below cover `ext/` and not just `src/`.
# Zygote is a weakdep of DataInterpolations but triggers no extension, so loading it
# would add nothing.
using ChainRulesCore, Makie, Mooncake, Optim, SparseConnectivityTracer, Symbolics

# Names DataInterpolations owns and its own extensions use. ExplicitImports already
# permits a package to use its own internals, but it decides "internal" via
# `Base.moduleroot`, and an extension module's root is itself rather than the package
# it extends, so every such use is reported. These are all genuinely private
# implementation details -- promoting them to public API would be wrong -- so they are
# allowed here instead.
const DATAINTERPOLATIONS_INTERNALS = (
    :AbstractInterpolation, :CurvefitCache, :DerivativeNotFoundError,
    :ExtrapolationError, :IntegralNotFoundError, :_interpolate, :derivative, :get_idx,
    :get_show, :integral, :munge_data, :to_plottable, :warn_if_ill_conditioned,
)

run_qa(
    DataInterpolations;
    ei_kwargs = (;
        # `@enumx ExtrapolationType ...` generates a submodule with dynamic includes
        # that ExplicitImports cannot statically analyze; allow it to be unanalyzable
        # rather than failing the import-graph checks.
        no_implicit_imports = (; allow_unanalyzable = (DataInterpolations.ExtrapolationType,)),
        no_stale_explicit_imports = (;
            allow_unanalyzable = (DataInterpolations.ExtrapolationType,),
            # DataInterpolationsSparseConnectivityTracerExt defines its tracer overloads
            # with `@eval function (interp::$(I){uType})(...)` over a list of
            # interpolation type names, so these imports are only used inside
            # interpolated macro output where static analysis cannot see them.
            ignore = (
                :AkimaInterpolation, :BSplineApprox, :BSplineInterpolation,
                :ConstantInterpolation, :CubicHermiteSpline, :CubicSpline,
                :LagrangeInterpolation, :LinearInterpolation, :QuadraticInterpolation,
                :QuadraticSpline, :QuinticHermiteSpline,
            ),
        ),
        all_qualified_accesses_are_public = (;
            ignore = (
                DATAINTERPOLATIONS_INTERNALS...,
                # ForwardDiff.{Dual,derivative,value} are not marked public in
                # ForwardDiff, so the public-API check flags them; keep them allowed
                # until ForwardDiff declares them public.
                :Dual, :value,
                # Mooncake's tangent representation, which any package defining
                # Mooncake rules has to reach for; not public in Mooncake.
                :FData, :NoRData, :RData, :increment_and_get_rdata!,
                # Makie's plot-recipe hooks: the documented way to teach Makie about a
                # new type, but not declared public in Makie.
                :SpecApi, :plottype,
                # Optim result accessor; not public in Optim.
                :minimizer,
                # SymbolicUtils' symtype/shape promotion interface, which the methods
                # registered by `@register_symbolic` must extend; not public in
                # SymbolicUtils.
                :ShapeT, :ShapeVecT, :TypeT, :is_array_shape, :promote_shape,
                :promote_symtype,
            ),
        ),
        all_explicit_imports_are_public = (;
            ignore = (
                DATAINTERPOLATIONS_INTERNALS...,
                # Mooncake's rule-authoring interface (`@from_chainrules`,
                # `@zero_adjoint` and the contexts they dispatch on); not public in
                # Mooncake.
                :DefaultCtx, :MinimalCtx, Symbol("@from_chainrules"),
                Symbol("@zero_adjoint"),
                # SparseConnectivityTracer's tracer types plus the 1-to-1 overload
                # helpers a package must call to declare its own operators; not public
                # in SparseConnectivityTracer.
                :AbstractTracer, :Dual, :GradientTracer, :HessianTracer,
                :gradient_tracer_1_to_1, :hessian_tracer_1_to_1, :primal, :tracer,
            ),
        ),
    ),
)
