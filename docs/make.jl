using ProxTV
using Documenter

DocMeta.setdocmeta!(ProxTV, :DocTestSetup, :(using ProxTV); recursive = true)

makedocs(
  modules = [ProxTV],
  sitename = "ProxTV.jl",
  format = Documenter.HTML(),
  repo = "https://github.com/nathanemac/ProxTV.jl",
  authors = "Nathan Allaire <nathan.allaire@polymtl.ca> and contributors",
  pages = [
    "index.md"
    [
      file for file in readdir(joinpath(@__DIR__, "src")) if
      file != "index.md" && endswith(file, ".md")
    ]
  ],
)

deploydocs(
  repo = "github.com/nathanemac/ProxTV.jl",
  devbranch = "main",
  push_preview = false,  # désactive les previews
)
