using SciMLTesting, DataInterpolations, Test

run_qa(
    DataInterpolations;
    explicit_imports = true,
    ei_kwargs = (;
        # `@enumx ExtrapolationType ...` generates a submodule with dynamic includes
        # that ExplicitImports cannot statically analyze; allow it to be unanalyzable
        # rather than failing the import-graph checks.
        no_implicit_imports = (; allow_unanalyzable = (DataInterpolations.ExtrapolationType,)),
        no_stale_explicit_imports = (; allow_unanalyzable = (DataInterpolations.ExtrapolationType,)),
        # Non-public qualified accesses we cannot annotate ourselves:
        #   * ForwardDiff.{Dual,derivative,value} are not marked public in ForwardDiff.
        #   * Base.require_one_based_indexing is a Base internal still flagged as
        #     non-public on Julia 1.11; newer Julia (>= 1.12) marks it public.
        all_qualified_accesses_are_public = (;
            ignore = (:Dual, :derivative, :value, :require_one_based_indexing),
        ),
    ),
)
