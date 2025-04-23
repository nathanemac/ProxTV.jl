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

## Understanding Total Variation

Total Variation (TV) regularization helps to denoise signals and images while preserving sharp edges or transitions. It penalizes the sum of absolute differences between adjacent values in the signal.

The p-norm version generalizes this to:

$$\min_x \frac{1}{2}\|x-y\|_2^2 + \lambda\sum_{i=1}^{n-1} |x_{i+1} - x_i|^p$$

- For p=1: L1-norm (promotes piecewise constant solutions)
- For p=2: L2-norm (promotes piecewise linear solutions)
- For other p values: Custom behaviors between L1 and L2

## Different Norms

ProxTV.jl supports different p-norms for TV regularization:

```julia
# L1 norm (default)
ProxTV.TV(y, lambda, x, 1.0)

# L2 norm (quadratic variation)
ProxTV.TV(y, lambda, x, 2.0)

# Custom p-norm (between L1 and L2)
ProxTV.TV(y, lambda, x, 1.5)
```

## Weighted Regularization

You can apply weights to different parts of the signal:

```julia
# Define weights for each pair of adjacent points
weights = ones(n-1)
weights[40:60] .= 3.0  # Stronger regularization in the middle segment

# Apply weighted TV regularization
ProxTV.TV(y, lambda, weights, x)
```

## More Examples

See the [Examples](10-examples.md) page for more detailed usage examples, including:

- TV regularization with different p-norms
- 2D total variation for image processing
- Weighted total variation
- Integration with ShiftedProximalOperators.jl
