using ProxTV
using Documenter

DocMeta.setdocmeta!(ProxTV, :DocTestSetup, :(using ProxTV); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers

makedocs(;
  modules = [ProxTV],
  authors = "Nathan Allaire <nathan.allaire@polymtl.ca> and contributors",
  repo = "https://github.com/nathanemac/ProxTV.jl/blob/{commit}{path}#{line}",
  sitename = "ProxTV.jl",
  format = Documenter.HTML(; canonical = "https://nathanemac.github.io/ProxTV.jl"),
  pages = [
    "index.md"
    [
      file for file in readdir(joinpath(@__DIR__, "src")) if
      file != "index.md" && splitext(file)[2] == ".md"
    ]
  ],
)

deploydocs(; repo = "github.com/nathanemac/ProxTV.jl")
