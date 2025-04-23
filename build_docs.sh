#!/bin/bash

# Build ProxTV.jl documentation

# Exit on error
set -e

echo "Setting up documentation environment..."
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'

echo "Building documentation..."
julia --project=docs docs/make.jl

echo "Documentation built successfully."
echo "You can view it by opening docs/build/index.html in your browser."

# Optional: open the documentation in the default browser
if [[ "$OSTYPE" == "darwin"* ]]; then
  echo "Opening documentation in browser..."
  open docs/build/index.html
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
  echo "To view documentation, run: xdg-open docs/build/index.html"
elif [[ "$OSTYPE" == "msys" ]]; then
  echo "To view documentation, run: start docs/build/index.html"
fi
