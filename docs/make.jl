using Documenter, DataInterpolations

ENV["GKSwstype"] = "100"

makedocs(
    modules = [DataInterpolations],
    sitename = "DataInterpolations.jl",
    warnonly = [:missing_docs],
    format = Documenter.HTML(),
    pages = [
        "index.md",
        "Methods" => "methods.md",
        "Interface" => "interface.md"
    ],
)

deploydocs(repo = "github.com/PumasAI/DataInterpolations.jl"; push_preview = true)
