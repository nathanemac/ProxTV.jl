# [Reference](@id reference)

## Contents

```@contents
Pages = ["95-reference.md"]
```

## Index

```@index
Pages = ["95-reference.md"]
```

## All Package Exports

```@autodocs
Modules = [ProxTV]
Order = [:constant, :function, :macro]
Filter = t -> !(t in [
    ProxTV.ModelFunction,
    ProxTV.ProxTVContext,
    ProxTV.InexactShiftedProximableFunction,
    ProxTV.NormLp,
    ProxTV.ShiftedNormLp,
    ProxTV.NormTVp,
    ProxTV.ShiftedNormTVp
]) && !startswith(string(t), "prox!") &&
       !startswith(string(t), "shifted") &&
       !startswith(string(t), "shift!")
```

## Internal Types and Functions

This section contains documentation for internal functions not typically used directly by users.

```@autodocs
Modules = [ProxTV]
Public = false
```
