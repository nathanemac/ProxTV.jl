```@meta
CurrentModule = ProxTV
```

# ProxTV.jl

ProxTV.jl is a Julia package that provides efficient implementations of proximal operators for total variation (TV) regularization with various p-norms.

## Overview

Total Variation regularization is widely used in signal and image processing for denoising, smoothing, and feature detection while preserving sharp transitions. This package offers:

- Efficient computation of TV proximal operators for 1D and 2D signals
- Support for different p-norms (L1, L2, and custom p-norms)
- Weighted regularization
- Integration with the ShiftedProximalOperators.jl package

## Installation

```julia
using Pkg
Pkg.add("ProxTV")
```

## Quick Example

```julia
using ProxTV

# Generate a noisy signal
n = 100
true_signal = vcat(zeros(30), ones(40), zeros(30))
noisy_signal = true_signal + 0.2 * randn(n)

# Denoise using TV-L1
lambda = 1.0
denoised = zeros(n)
ProxTV.TV(noisy_signal, lambda, denoised, 1.0)

# denoised now contains the TV-L1 regularized signal
```

## Documentation

- [Getting Started](00-getting-started.md) - Installation and basic usage
- [Examples](10-examples.md) - Practical examples demonstrating different use cases
- [API Reference](20-api-reference.md) - Complete reference of functions and types
- [Contributing](90-contributing.md) - Guidelines for contributing to the project
- [Developer Guide](91-developer.md) - Technical information for developers

## Package Features

- **Efficient Implementation**: Based on algorithms from the original ProxTV package
- **Versatile Interface**: Support for different norms and problem dimensions
- **Julia Integration**: Designed to work well with Julia's type system and other packages

## Citation

If you use ProxTV.jl in your research, please cite:

```bibtex
@software{allaire2025proxtv,
  author = {Allaire, Nathan},
  title = {ProxTV.jl: Efficient Proximal Operators for Total Variation Regularization},
  url = {https://github.com/nathanemac/ProxTV.jl},
  year = {2025}
}
```

## Contributors

```@raw html
<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
```
