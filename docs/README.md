# ProxTV.jl Documentation

This directory contains the documentation for the ProxTV.jl package.

## Documentation Structure

- `src/` - Contains the Markdown files that make up the documentation
- `make.jl` - The script that builds the documentation
- `Project.toml` - The project file for the documentation
- `serve_docs.jl` - Script to build and serve documentation locally

## Documentation Pages

- `index.md` - The home page with package overview
- `00-getting-started.md` - Installation and basic usage
- `10-examples.md` - Examples of using the package
- `20-api-reference.md` - API reference (manual documentation of functions)
- `90-contributing.md` - Contributing guidelines
- `91-developer.md` - Developer guide
- `95-reference.md` - API reference (automatically generated from docstrings)

## Building the Documentation Locally

To build the documentation locally, you have several options:

### Option 1: Using the build_docs.sh script

```bash
./build_docs.sh
```

This will install dependencies, build the documentation, and open it in your default browser.

### Option 2: Using the serve_docs.sh script

```bash
./serve_docs.sh
```

This will build the documentation and start a local web server at localhost:8000 where you can view it. The documentation will update automatically when you make changes to the source files.

### Option 3: Manual build

```julia
# First, install dependencies
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'

# Then build the documentation
julia --project=docs docs/make.jl
```

The built documentation will be in the `docs/build` directory.

## Deployment

The documentation is automatically deployed when changes are pushed to the main branch. This is handled by the GitHub workflow defined in `.github/workflows/Docs.yml`.

A GitHub action builds the documentation and deploys it to:
<https://nathanemac.github.io/ProxTV.jl/>

## Customization

Feel free to improve this documentation by:

1. Adding more examples
2. Improving explanations
3. Adding tutorials
4. Adding docstrings to functions in the source code
5. Fixing any errors or typos
