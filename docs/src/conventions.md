# Development conventions

LineCableModels follows the SciML formatter style. Run

```julia
using JuliaFormatter
format(".")
```

before committing.

Versions follow [Semantic Versioning](https://semver.org/). Public behavior is
kept compatible within a minor release; deprecations must provide a migration
path before removal.

Commit subjects use scoped Conventional Commits, start with a lowercase
description, and remain within 72 characters. For example:

```text
fix(engine): reject unsupported formulation options
```

Changes must include tests at the lowest useful level. Core tests must not load
optional integrations. CairoMakie is verified in a separate gate.
Examples in docstrings should be executable and self-contained; examples that
require fixtures, user interfaces, or external executables belong in integration
documentation instead.

Physical quantities state their SI units. A docstring for implemented
physically meaningful mathematics includes a `# Notes` section with the
equation and definitions used by the code. This requirement follows the
implementation, not a naming prefix.

Docstrings use `DocStringExtensions` abbreviations as the codebase-wide default
for generated signatures, type declarations, fields, and module inventories.
The complete templates, unit notation, mathematical requirements, and example
rules are documented in [Docstrings](@ref).
