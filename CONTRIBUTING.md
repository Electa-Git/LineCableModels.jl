# Contributing

Use Julia 1.12 and instantiate the package environment before making changes:

```julia
using Pkg
Pkg.instantiate()
Pkg.test()
```

Format maintained Julia files with `JuliaFormatter.format(".")`. Commit
subjects must satisfy `.gitlint`: a scoped Conventional Commit, a lowercase
description, and no more than 72 characters.

Keep pull requests focused. Add tests for changed behaviour and update public
documentation when an API changes. Optional plotting integrations must remain
outside core loading and must be checked with their dedicated workflows.
