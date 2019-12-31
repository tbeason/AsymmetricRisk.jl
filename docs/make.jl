using Documenter
using AsymmetricRisk

makedocs(
    sitename = "AsymmetricRisk.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    modules = [AsymmetricRisk],
    pages = ["Introduction" => "index.md","univariate.md","bivariate.md"]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/tbeason/AsymmetricRisk.jl.git",
)
