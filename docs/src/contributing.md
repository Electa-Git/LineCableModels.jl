# Contributing

Use Julia 1.12 and instantiate the package environment before making changes:

```julia
using Pkg
Pkg.instantiate()
Pkg.test()
```

Format maintained Julia files with `JuliaFormatter.format(".")`.
Use `type(scope): description` for commit subjects, with a lowercase imperative
description. Keep the complete subject within 72 characters. A breaking change
may use `type(scope)!: description`. Commit bodies are optional.

Install Gitlint and enable its native commit hook once per clone:

```bash
uv tool install gitlint-core
gitlint install-hook
```

The hook checks messages against `.gitlint` before creating a commit.
Gitlint keeps its default exemptions for merge, revert and fixup commits.

Keep pull requests focused. Add tests for changed behavior and update public
documentation when an API changes. Optional plotting integrations must remain
outside core loading and must be checked with their dedicated workflows.
