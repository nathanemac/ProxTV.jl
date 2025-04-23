```@meta
CurrentModule = ProxTV
```

# API Reference

This page contains reference documentation for the ProxTV.jl package.

## Core Functions

### TV

```julia
TV(y, lambda, x, [p=1.0])
```

Computes the proximal operator of the p-norm total variation. This solves:

$$\min_x \frac{1}{2}\|x-y\|_2^2 + \lambda\sum_{i=1}^{n-1} |x_{i+1} - x_i|^p$$

Parameters:

- `y`: Input signal
- `lambda`: Regularization parameter
- `x`: Output buffer (will be overwritten)
- `p`: Norm parameter (default: 1.0)

### TV with Weights

```julia
TV(y, lambda, w, x, [p=1.0])
```

Computes the proximal operator of the weighted p-norm total variation. This solves:

$$\min_x \frac{1}{2}\|x-y\|_2^2 + \lambda\sum_{i=1}^{n-1} w_i|x_{i+1} - x_i|^p$$

Parameters:

- `y`: Input signal
- `lambda`: Regularization parameter
- `w`: Weights vector (length n-1)
- `x`: Output buffer (will be overwritten)
- `p`: Norm parameter (default: 1.0)

### TVp_norm

```julia
TVp_norm(x, p)
```

Computes the p-norm total variation of a signal:

$$\sum_{i=1}^{n-1} |x_{i+1} - x_i|^p$$

Parameters:

- `x`: Input signal
- `p`: Norm parameter

## Integration with ShiftedProximalOperators

ProxTV.jl integrates with the ShiftedProximalOperators.jl package to provide efficient proximal operators for optimization algorithms.

### Types

- `ProxTVContext`: Context for TV proximal operators
- `InexactShiftedProximableFunction`: Base type for proximal functions without closed-form solutions
- `NormTVp`: Total variation norm with p-norm
- `ShiftedNormTVp`: Shifted version of the TV p-norm

### Methods

- `prox!`: Computes the proximal operator
- `shifted`: Creates a shifted version of a proximal function
- `shift!`: Shifts a proximal function in-place

## Function Index

```@index
Modules = [ProxTV]
```

## Function Details

Here are the main functions and their expected parameters:

### TV Function Signature

```julia
TV(y, lambda, x, [p=1.0])
```

Computes the proximal operator of the p-norm total variation. This solves:

$$\min_x \frac{1}{2}\|x-y\|_2^2 + \lambda\sum_{i=1}^{n-1} |x_{i+1} - x_i|^p$$

Parameters:

- `y`: Input signal
- `lambda`: Regularization parameter
- `x`: Output buffer (will be overwritten)
- `p`: Norm parameter (default: 1.0)

### Weighted TV Function Signature

```julia
TV(y, lambda, w, x, [p=1.0])
```

Computes the proximal operator of the weighted p-norm total variation. This solves:

$$\min_x \frac{1}{2}\|x-y\|_2^2 + \lambda\sum_{i=1}^{n-1} w_i|x_{i+1} - x_i|^p$$

Parameters:

- `y`: Input signal
- `lambda`: Regularization parameter
- `w`: Weights vector (length n-1)
- `x`: Output buffer (will be overwritten)
- `p`: Norm parameter (default: 1.0)

### TVp_norm Function Signature

```julia
TVp_norm(x, p)
```

Computes the p-norm total variation of a signal:

$$\sum_{i=1}^{n-1} |x_{i+1} - x_i|^p$$

Parameters:

- `x`: Input signal
- `p`: Norm parameter
