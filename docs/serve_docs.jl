#!/usr/bin/env julia

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

Pkg.develop(PackageSpec(path = dirname(@__DIR__)))

using LiveServer

# Build docs
include("make.jl")

# Serve docs
serve(dir = joinpath(@__DIR__, "build"), host = "0.0.0.0", port = 8000)
