```@meta
CurrentModule = ProxTV
```

# Examples

This page contains examples of how to use the ProxTV.jl package for various total variation proximal operator computations.

## 1D Total Variation

### Basic 1D TV with L1 norm

```julia
using ProxTV

# Generate random signal
n = 100
y = cumsum(randn(n)) # Random walk signal
lambda = 2.0 # Regularization parameter
x = zeros(n) # Output buffer

# Compute proximal operator of TV-L1 norm
ProxTV.TV(y, lambda, x, 1.0)

# x now contains the solution
```

### L2 norm (Quadratic Variation)

```julia
using ProxTV

n = 100
y = cumsum(randn(n))
lambda = 2.0
x = zeros(n)

# Compute proximal operator with L2 norm
ProxTV.TV(y, lambda, x, 2.0)
```

### Custom p-norm

```julia
using ProxTV

n = 100
y = cumsum(randn(n))
lambda = 2.0
x = zeros(n)
p = 1.5 # Custom p-value between 1 and 2

# Compute proximal operator with custom p-norm
ProxTV.TV(y, lambda, x, p)
```

## Weighted Total Variation

```julia
using ProxTV

n = 100
y = cumsum(randn(n))
w = ones(n-1) # Weights
w[40:60] .= 3.0 # Increase regularization in the middle
lambda = 1.0
x = zeros(n)

# Compute proximal operator with weighted TV
ProxTV.TV(y, lambda, w, x)
```

## Integration with ShiftedProximalOperators

ProxTV.jl integrates with the ShiftedProximalOperators.jl package:

```julia
using ProxTV
using ShiftedProximalOperators

# Create a TV norm function
n = 100
tv_func = NormTVp(n, 1.0, 1.0) # n, lambda, p

# Create input signal
y = cumsum(randn(n))

# Compute prox directly
x = similar(y)
prox!(x, tv_func, y, 1.0)
```
