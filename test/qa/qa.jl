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
        # ForwardDiff.{Dual,derivative,value} are not marked public in ForwardDiff,
        # so the public-API check flags them; keep them allowed until ForwardDiff
        # declares them public.
        all_qualified_accesses_are_public = (; ignore = (:Dual, :derivative, :value)),
    ),
)
