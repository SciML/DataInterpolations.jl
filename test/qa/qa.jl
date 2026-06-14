using DataInterpolations, Aqua
@testset "Aqua" begin
    Aqua.find_persistent_tasks_deps(DataInterpolations)
    Aqua.test_ambiguities(DataInterpolations, recursive = false)
    Aqua.test_deps_compat(DataInterpolations)
    Aqua.test_piracies(DataInterpolations)
    Aqua.test_project_extras(DataInterpolations)
    Aqua.test_stale_deps(DataInterpolations)
    Aqua.test_unbound_args(DataInterpolations)
    Aqua.test_undefined_exports(DataInterpolations)
end
