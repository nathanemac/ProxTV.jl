# ProxTV

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://nathanemac.github.io/ProxTV.jl/stable)
[![Build Status](https://github.com/nathanemac/ProxTV.jl/workflows/Test/badge.svg)](https://github.com/nathanemac/ProxTV.jl/actions)
[![Test workflow status](https://github.com/nathanemac/ProxTV.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/nathanemac/ProxTV.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Lint workflow Status](https://github.com/nathanemac/ProxTV.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/nathanemac/ProxTV.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/nathanemac/ProxTV.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/nathanemac/ProxTV.jl/actions/workflows/Docs.yml?query=branch%3Amain)

[![Coverage](https://codecov.io/gh/nathanemac/ProxTV.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/nathanemac/ProxTV.jl)
[![DOI](https://zenodo.org/badge/DOI/FIXME)](https://doi.org/FIXME)

[![All Contributors](https://img.shields.io/github/all-contributors/nathanemac/ProxTV.jl?labelColor=5e1ec7&color=c0ffee&style=flat-square)](#contributors)

ProxTV.jl is a Julia package that provides a collection of exact and inexact proximal operators. This includes the Total Variation (TV) regularization with any p-norm.

This package is a Julia implementation of the ProxTV package for MATLAB and Python which is available [here](https://github.com/albarji/proxTV).

## How to Use

The package is designed to be easy to use and to provide a consistent interface for all the implemented proximal operators.

### Installation

You can install ProxTV.jl using the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia
pkg> add ProxTV
```

### Example

Here is an example of how to use ProxTV.jl to compute the proximal operator of the Total Variation (TV) regularization with a p-norm on a 1D signal.

```julia

using ProxTV

n = rand(10:100)
lambda = 0.15
y = rand(n)
x = zeros(n)
p = 1.32 # inexact prox computation : no closed-form for p = 1.32
ProxTV.TV(y, lambda, x, p)
```

Other examples can be found in the [documentation](https://nathanemac.github.io/ProxTV.jl/stable).

## How to Cite

If you use ProxTV.jl in your work, please cite using the reference given in [CITATION.cff](https://github.com/nathanemac/ProxTV.jl/blob/main/CITATION.cff).

## Contributing

If you want to make contributions of any kind, please first that a look into our [contributing guide directly on GitHub](docs/src/90-contributing.md).

---

### Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
