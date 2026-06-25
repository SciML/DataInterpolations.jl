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
        # Non-public names accessed by qualified-access that vary by package version:
        #   * ForwardDiff.{Dual,derivative,value} are not marked public in ForwardDiff.
        #   * Base.{front,require_one_based_indexing} are Base internals flagged as
        #     non-public on Julia 1.10 (the LTS lane); newer Julia marks them public.
        all_qualified_accesses_are_public = (;
            ignore = (:Dual, :derivative, :value, :front, :require_one_based_indexing),
        ),
    ),
)
