#!/bin/bash

# Check that ProxTV.jl documentation builds correctly

# Exit on error
set -e

echo "Setting up documentation environment..."
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'

echo "Building documentation..."
julia --project=docs docs/make.jl

# Check if docs built successfully
if [ -f "docs/build/index.html" ]; then
  echo "✅ Documentation built successfully!"
  echo "You can view it by opening docs/build/index.html in your browser."
else
  echo "❌ Documentation build failed!"
  exit 1
fi
