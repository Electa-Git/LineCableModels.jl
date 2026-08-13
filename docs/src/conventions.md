# Development conventions

LineCableModels follows the SciML formatter style. Run

```julia
using JuliaFormatter
format(".")
```

before committing. The vendored GetDP frontend under
`ext/fem/getdp_frontend/` is a provenance-pinned compatibility snapshot; its
documented loader and namespace adaptations are excluded from formatting.

Versions follow [Semantic Versioning](https://semver.org/). Public behavior is
kept compatible within a minor release; deprecations must provide a migration
path before removal.

Commit subjects use scoped Conventional Commits, start with a lowercase
description, and remain within 72 characters. For example:

```text
fix(fem): resolve getdp from the configured executable
```

Changes must include tests at the lowest useful level. Core tests must not load
optional integrations. CairoMakie and FEM/GetDP are verified in separate gates.
Examples in docstrings should be executable and self-contained; examples that
require fixtures, user interfaces, or external executables belong in integration
documentation instead.

Physical quantities state their SI units. A docstring for implemented
physically meaningful mathematics includes a `# Notes` section with the
equation and definitions used by the code. This requirement follows the
implementation, not a naming prefix.
