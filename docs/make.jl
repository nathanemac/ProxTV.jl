using ProxTV
using Documenter

DocMeta.setdocmeta!(ProxTV, :DocTestSetup, :(using ProxTV); recursive = true)

makedocs(
  modules = [ProxTV],
  sitename = "ProxTV.jl",
  format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    canonical = "https://nathanemac.github.io/ProxTV.jl",
    assets = ["assets/favicon.ico"],
  ),
  repo = "https://github.com/nathanemac/ProxTV.jl",
  authors = "Nathan Allaire <nathan.allaire@polymtl.ca> and contributors",
  pages = [
    "Home" => "index.md",
    "Getting Started" => "00-getting-started.md",
    "Examples" => "10-examples.md",
    "Build your own prox API" => "20-build-your-own-prox-api.md",
    "Contributing" =>
      ["Guidelines" => "90-contributing.md", "Developer Guide" => "91-developer.md"],
    "References" => "95-reference.md",
  ],
)

deploydocs(
  repo = "github.com/nathanemac/ProxTV.jl",
  devbranch = "main",
  push_preview = false,  # désactive les previews
)
