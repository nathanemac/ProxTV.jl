```@meta
CurrentModule = ProxTV
```

# Getting Started

`ProxTV.jl` is a Julia package that provides efficient implementations of proximal operators for total variation (TV)_p regularization and p-norms.
It is designed to be used in conjunction with [ShiftedProximalOperators.jl](https://github.com/JuliaSmoothOptimizers/ShiftedProximalOperators.jl) for use in [RegularizedOptimization.jl](https://github.com/JuliaSmoothOptimizers/RegularizedOptimization.jl).

## Installation

You can install `ProxTV.jl` using the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia
pkg> add ProxTV
```

Alternatively, you can use the Pkg API:

```julia
using Pkg
Pkg.add("ProxTV")
```

## Basic Usage

The basic functionality of `ProxTV.jl` is to compute the proximal operator of the total variation regularization. The simplest usage is with 1D signals:

```julia
using ProxTV

n = 1000
y = rand(n)
x = similar(y)
λ = 0.1 # regularization parameter
p = 1.37 # p-norm parameter
ν = 1e-2 # proximal parameter

h = NormTVp(λ, p, n) # builds the norm object
prox!(x, h, y, ν) # computes the proximal operator in place

```

This computes the solution to:

$$\min_x \frac{1}{2\nu}\|x-y\|_2^2 + \lambda\sum_{i=1}^{n-1} |x_{i+1} - x_i|^p$$

which is the proximal operator of the TVp-norm with parameter λ.

## Understanding Total Variation

Total Variation (TV) regularization helps to denoise signals and images while preserving sharp edges or transitions. It penalizes the sum of absolute differences between adjacent values in the signal.

## More Examples

See the [Examples](10-examples.md) page for more detailed usage examples, including:

- Lp-norm regularization with different p-norms
- Shifted versions of both the TVp-norm and the Lp-norm
- Integration with `ShiftedProximalOperators.jl` for use in `RegularizedOptimization.jl`

## Build your prox API

If you want to integrate your proximal operator in `RegularizedOptimization.jl`, see the [Build your own prox API](20-build-your-own-prox-api.md) page.
