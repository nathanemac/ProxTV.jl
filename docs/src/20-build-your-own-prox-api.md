```@meta
CurrentModule = ProxTV
```

# Build your own prox API

!!! note "This page is under construction, and should be used when `iR2N` and `iR2` will be merged into `RegularizedOptimization.jl`".

## Motivation

Assume you want to solve a regularized optimization problem of the form:

```math
\min_x f(x) + h(x),
```

where `f` is smooth and `h` is possibly non-smooth.

Methods such as [R2N](https://github.com/JuliaSmoothOptimizers/RegularizedOptimization.jl/blob/master/src/R2N.jl) or [R2](https://github.com/JuliaSmoothOptimizers/RegularizedOptimization.jl/blob/master/src/R2_alg.jl) and their inexact variants `iR2N` and `iR2` solve this problem by iteratively minimizing a model of the form:

```math
\min_s \frac{1}{2}\nu^{-1}\|s + \nu \nabla f(x_k) - \|_2^2 + \psi(s),
```

which is $\textrm{prox}_{\nu \psi}(-\nu \nabla f(x_k))$.

Although `ProxTV.jl`, `ShiftedProximalOperators.jl` and `IRBP.jl` provide a various set of proximal operators, you might need to use one that is not implemented yet.

In this case, you can build your own proximal operator by following the steps below.

## Step by step guide

### Step 1: Define the `Context`

The `Context` is a structure that contains all buffers, parameters and structures that are used to compute the proximal operator.
For instance, the `Context` for the proximal operators in `ProxTV.jl` is defined as follows:

```julia
mutable struct ProxTVContext{F}
  hk::Float64
  h_symb::Symbol
  ∇fk::Vector{Float64}
  h_fun::F
  p::Real
  κξ::Float64
  shift::Vector{Float64}
  s_k_unshifted::Vector{Float64}
  dualGap::Float64
  prox_stats::Vector{Int64}
  callback_pointer::Ptr{Cvoid}
  info::Vector{Float64}
  temp_x::Vector{Float64}
  y_shifted::Vector{Float64}
  s::Vector{Float64}
  s_k::Vector{Float64}
end
```

You might need a totally different structure, but the idea is that it contains all the information you need along the iterations to compute the proximal operator.
You might even not need a context at all, and compute the proximal operator directly in the `prox!` function (which is the case for most operators in [ShiftedProximalOperators.jl](https://github.com/JuliaSmoothOptimizers/ShiftedProximalOperators.jl/tree/master/src)).

!!! note "Allocation-free"
    Note that `ProxTV.jl` is memory allocation free, essentially thanks to the `Context` structure.

### Step 2: Define the regularizer `h`

The regularizer is a callable object that contains at least the context (if needed) and the regularization parameter.

In the case of the `TVp-Norm` in `ProxTV.jl`, the regularizer is defined as follows:

```julia
h = NormTVp(λ, p, n)
```

The context is initialized when building the regularizer and accessible through `h.context`.

### Step 3: Define the shifted regularizer

In `RegularizedOptimization.jl` methods, a shifted version of the regularizer is used.
The shifted regularizer is a callable object `ψ` that contains **at least** the regularizer and the shift.

In the case of the `TVp-Norm` in `ProxTV.jl`, the shifted regularizer is initialized as follows:

```julia
ψ = shifted(h, shift)
```

!!! tip "`h` and `ψ` are callable"
    Make sure the regularizer `h` and the shifted regularizer `ψ` are callable as `h(x)` and `ψ(x)`.

Eventually, you need to define a method `shift!` that updates the shift in place: `shift!(ψ, xk)`.

### Step 4: Define the proximal operator

Now that structures are defined, you need to implement the proximal operator functions for both the regularizer and the shifted regularizer.
They must follow the signature:

```julia
function prox!(y, h::YourRegularizerType, q, ν)
    # implement the proximal operator for h
end

function prox!(y, ψ::YourShiftedRegularizerType, q, ν)
    # implement the proximal operator for ψ (which is the proximal operator of h shifted by ψ.shift)
end
```

- `y` is the buffer in which the result is stored
- `h` is the regularizer
- `q` is the point at which the proximal operator is computed
- `ν` is the regularization parameter

### Step 5: Update the context

Your `Context` object might need to be updated at each iteration.
For instance, the `ProxTVContext` object needs to update the value `h(xk)`, the gradient `∇f(xk)` and the shift to compute the proximal operator.

At each step in `iR2N` and `iR2`, the method `update_prox_context!` is called to update the context.

You need to define a method `update_prox_context!` that suits your needs and should follow the signature:

```julia
function update_prox_context!(solver, stats, ψ::YourShiftedRegularizerType)
    # update the context for ψ
end
```

It can be a no-op if you don't need to update the context or don't have a context.

### Conclusion

🎉 Tada!  If you followed all the steps above, you can now use your proximal operator in `RegularizedOptimization.jl`!

Here is a simple example of how to use your proximal operator with `iR2N`:

```julia
model, nls_model, sol = bpdn_model(1)
nlp = LBFGSModel(model)
n = nlp.meta.nvar
h = YourRegularizer(0.1, 1.6, n)
reg_nlp = RegularizedNLPModel(nlp, h)
solver = iR2NSolver(reg_nlp)
stats = RegularizedExecutionStats(reg_nlp)
RegularizedOptimization.solve!(
  solver,
  reg_nlp,
  stats,
  kwargs...
)
```

If you want to contribute to `ProxTV.jl` or `ShiftedProximalOperators.jl` with your fresh new proximal operator, feel free to open a PR!

If you have any question, feel free to open an issue on the [Github repository](https://github.com/nathanemac/ProxTV.jl) or to contact [Nathan Allaire](mailto:nathan.allaire@polymtl.ca).
