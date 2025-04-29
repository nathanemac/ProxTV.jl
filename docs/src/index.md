```@meta
CurrentModule = ProxTV
```

# ProxTV.jl

`ProxTV.jl` is a Julia package that provides efficient implementations of proximal operators for total variation (TV)_p regularization and p-norms.
It is designed to be used in conjunction with [ShiftedProximalOperators.jl](https://github.com/JuliaSmoothOptimizers/ShiftedProximalOperators.jl) for use in [RegularizedOptimization.jl](https://github.com/JuliaSmoothOptimizers/RegularizedOptimization.jl).

## Overview

Total Variation regularization is widely used in signal and image processing for denoising, smoothing, and feature detection while preserving sharp transitions. This package offers:

- Efficient computation of TV proximal operators for 1D, 2D, and nD signals
- Support for different p-norms (L1, L2, and custom p-norms)
- Weighted regularization
- Integration with the [ShiftedProximalOperators.jl](https://github.com/JuliaSmoothOptimizers/ShiftedProximalOperators.jl) package

## Installation

```julia
using Pkg
Pkg.add("ProxTV")
```

## Documentation

- [Getting Started](00-getting-started.md) - Installation and basic usage
- [Examples](10-examples.md) - Practical examples demonstrating different use cases
- [Build your own prox API](20-build-your-own-prox-api.md) - Use your own proximal operator in `RegularizedOptimization.jl`
- [Contributing](90-contributing.md) - Guidelines for contributing to the project
- [Developer Guide](91-developer.md) - Technical information for developers

## Package Features

- **Efficient Implementation**: Based on algorithms from the original [ProxTV](https://github.com/albarji/proxTV) package
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
