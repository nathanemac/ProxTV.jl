#!/bin/bash

# Build and serve ProxTV.jl documentation

# Exit on error
set -e

echo "Setting up documentation environment..."
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'

echo "Building and serving documentation..."
julia --project=docs docs/serve_docs.jl

# This script will keep running until interrupted with Ctrl+C
