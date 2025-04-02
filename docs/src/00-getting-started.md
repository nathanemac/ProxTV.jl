```@meta
CurrentModule = ProxTV
```

# Getting Started

ProxTV.jl is a Julia package that provides efficient implementations of proximal operators for total variation (TV)_p regularization and p-norms.

## Installation

You can install ProxTV.jl using the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia
pkg> add ProxTV
```

Alternatively, you can use the Pkg API:

```julia
using Pkg
Pkg.add("ProxTV")
```

## Basic Usage

The basic functionality of ProxTV.jl is to compute the proximal operator of the total variation regularization. The simplest usage is with 1D signals:

```julia
using ProxTV

# Generate a random signal
n = 100
y = cumsum(randn(n)) # Random walk signal (noisy)

# Set regularization parameter
lambda = 2.0

# Allocate output buffer
x = zeros(n)

# Compute proximal operator of TV-L1 norm (default)
ProxTV.TV(y, lambda, x)

# x now contains the denoised signal
```

This computes the solution to:

$$\min_x \frac{1}{2}\|x-y\|_2^2 + \lambda\sum_{i=1}^{n-1} |x_{i+1} - x_i|$$

which is the proximal operator of the L1 total variation norm with parameter λ.

## More Examples

See the [Examples](10-examples.md) page for more detailed usage examples, including:

- TV regularization with different p-norms
- 2D total variation for image processing
- Weighted total variation
- Integration with ShiftedProximalOperators.jl
