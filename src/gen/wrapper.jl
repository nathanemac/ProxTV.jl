# Script to parse proxTV headers and generate Julia wrappers.
using proxTV_jll
using Clang
using Clang.Generators
using JuliaFormatter

function main()
  cd(@__DIR__)
  include_dir = joinpath(proxTV_jll.artifact_dir, "include")
  headers = [
    joinpath(include_dir, header) for
    header in readdir(include_dir) if endswith(header, ".h")
  ]

  options = load_options(joinpath(@__DIR__, "proxtv.toml"))
  options["general"]["output_file_path"] = joinpath("..", "libproxtv.jl")
  options["general"]["output_ignorelist"] =
    ["mxGetInf", "sign", "min", "max", "dpttrs_", "dpttrf_"]  # example: "RC_OK"

  args = get_default_args()
  push!(args, "-I$include_dir")
  push!(args, "-DNOMATLAB")

  ctx = create_context(headers, args, options)
  build!(ctx)

  path = options["general"]["output_file_path"]
  code = read(path, String)
  code = "using proxTV_jll\n\n" * code
  code = replace(code, "Cdouble" => "Float64")
  code = replace(code, "Cint" => "Int32")
  write(path, code)

  format_file(path, YASStyle())
  return nothing
end

# If we want to use the file as a script with `julia wrapper.jl`
if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
