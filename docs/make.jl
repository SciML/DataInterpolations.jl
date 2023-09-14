using Documenter, DataInterpolations

ENV["GKSwstype"] = "100"

makedocs(
    modules = [DataInterpolations],
    sitename = "DataInterpolations.jl",
    strict = [:doctest, :linkcheck, :parse_error, :example_block],
    format = Documenter.HTML(),
    pages = [
        "index.md"
        "Tutorial" => ["tutorial.md"]
    ],
)

deploydocs(repo = "github.com/PumasAI/DataInterpolations.jl"; push_preview = true)
