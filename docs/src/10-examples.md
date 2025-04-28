```@meta
CurrentModule = ProxTV
```

# Examples

This page contains examples of how to use the `ProxTV.jl` package for various total variation proximal operator computations.

## Exact proximal operators

```julia

n = 100
y = cumsum(randn(n))
λ = 0.1
info = zeros(3) # to store algorithm information (iterations, residual, time)
ws = ProxTV.newWorkspace(n) # allocates workspace for C++ calls

x1 = similar(y)
x2 = similar(y)
xinf = similar(y)

ProxTV.PN_LP1(y, λ, x1, info, n) # L1 norm
ProxTV.PN_LP2(y, λ, x2, info, n) # L2 norm
ProxTV.PN_LPinf(y, λ, xinf, info, n, ws) # L-inf norm
```

## Inexact proximal operators

```julia
x = zeros(n)
p = 1.6 # Custom p-value (can be anything between 1 and Inf)
objGap = 1e-8 # Primal-dual gap (stopping criterion)

# Compute proximal operator with custom p-norm
ProxTV.PN_LPp(y, λ, x, p, objGap)
```

## Integration with `ShiftedProximalOperators.jl`

### Norm objects

`ProxTV.jl` follows the design of `ShiftedProximalOperators.jl`.
You can notably define a `Norm` object, `shift` it, call it as a function, evaluate its proximal operator, etc.

```julia
n = 100
λ = 0.1
p = 1.6
x = rand(n)

# Define the norm objects.
# The first argument is the regularization parameter, the second is the p-norm parameter, the third is the problem size.
h_lp = NormLp(λ, p, n)
h_tvp = NormTVp(λ, p, n)

# Norm objects are callable.
h_lp(x) # evaluates the Lp-norm of x times λ
h_tvp(x) # evaluates the TVp-norm of x times λ
```

Inside `h`, the `Context` object is stored.
It contains all buffers and parameters for the proximal operator computation and is specific to the chosen norm.
It is automatically managed and updated by `ProxTV.jl` and does not need to be modified.

### Shifted norm objects

```julia
shift = ones(n)

ψ_lp = shifted(h_lp, shift)
ψ_tvp = shifted(h_tvp, shift)

# Shifted norm objects are callable.
ψ_lp(x) # evaluates the Lp-norm of x + shift, times λ
ψ_tvp(x) # evaluates the TVp-norm of x + shift, times λ

# The shift can be modified in place.
second_shift = 5*ones(n)
shift!(ψ_lp, second_shift)
shift!(ψ_tvp, second_shift)
# Now, calling ψ_lp and ψ_tvp evaluates the Lp-norm and TVp-norm of x + second_shift, times λ.
```

The `shifted` function returns a `ShiftedNorm` object, a subtype of `InexactShiftedProximableFunction`, that extends `ShiftedProximableFunction` from `ShiftedProximalOperators.jl` .

### Proximal operator

The proximal operator of a (shifted or not) norm object can be computed in place by calling the `prox!` function.

```julia
res = similar(x) # allocate a vector to store the result
prox!(res, h_lp, x, 1.0)
prox!(res, h_tvp, x, 1.0)

prox!(res, ψ_lp, x, 1.0) # also works with shifted norm objects
prox!(res, ψ_tvp, x, 1.0)
```

The `prox!` function is a wrapper around the C++ proximal operator functions of `proxTV` and follows the signature of the `ShiftedProximalOperators.jl` `prox!` function.
It takes as arguments the vector to store the result, the (shifted) norm object, the vector at which to evaluate the proximal operator, and the regularization parameter (ν).

## Usage in `RegularizedOptimization.jl`

The package is designed to be used in conjunction with `RegularizedOptimization.jl` and `ShiftedProximalOperators.jl`.

Here is a simple example of how to use the package to solve a complex regularized optimization problem.

The following example solves the problem:

```math
\min \|Ax - b\|_2^2 + \lambda \|x\|_p,
```

where `A` is a matrix, `b` is a vector, and `x` is the variable to optimize over.
The problem is solved using the `iR2N` algorithm from `RegularizedOptimization.jl`.

```julia
using RegularizedOptimization, ProxTV, NLPModels, ADNLPModels, RegularizedProblems, NLPModelsModifiers

model, nls_model, sol = bpdn_model(1)
nlp = LBFGSModel(model)
h = NormLp(0.1, 1.6, nlp.meta.nvar)
reg_nlp = RegularizedNLPModel(nlp, h)
solver = iR2NSolver(reg_nlp)
stats = RegularizedExecutionStats(reg_nlp)
RegularizedOptimization.solve!(
  solver,
  reg_nlp,
  stats,
  ν=1.0,
  atol=1e-6,
  rtol=1e-6,
  max_iter=1000,
  verbose=1
)
```

For further concerns, please open an issue on the Github repository or contact [Nathan Allaire](mailto:nathan.allaire@polymtl.ca).
