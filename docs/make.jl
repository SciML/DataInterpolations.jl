using Documenter, DataInterpolations

ENV["GKSwstype"] = "100"

makedocs(
    modules = [DataInterpolations],
    sitename = "DataInterpolations.jl",
    warnonly = [:missing_docs],
    format = Documenter.HTML(
        analytics = "UA-90474609-3",
        assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/DataInterpolations/stable/",
    ),
    pages = [
        "index.md",
        "Methods" => "methods.md",
        "Interface" => "interface.md"
    ],
)

deploydocs(repo = "github.com/SciML/DataInterpolations.jl"; push_preview = true)
