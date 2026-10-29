# Contributing

LineCableModels is a research package whose structure is enforced by tests. A contribution
is welcome when it keeps that structure. This page gives the principles, the environment,
the checks and the form of a pull request. The developer section of the documentation
describes the architecture, the rules and the conventions in full.

## Principles

- A change is specified before it is written. The specification lists the types, the
  signatures, the files and the tests it adds, and the sources of the physics it ports, with
  their equations and references.
- A change adds only the names, tests, checks and files that it requires. What does not
  fit the architecture is discussed first, never worked around.
- A new name is a concept. It is proposed with its meaning and its terms, confirmed by the
  maintainers, and entered in `test/architecture/names.toml` in the same change as its
  definition.
- Every test calls feature code and checks a feature result. A tolerance, an assertion or a
  guard is never loosened to make a change pass.
- Numbers come before decisions. A design choice that rests on a measurement states the
  measurement, and a physical result is validated against the frozen reference cases under
  `test/references/`.
- Every fact a change relies on is written in the pull request, never left in a conversation.
- A new behavior is a method of an existing verb on an owned type before it is a new
  function. The verbs of the package are few and closed, and its types are open.
- A design choice rests on a measurement and on the physics. A passing check or an existing
  page is never the reason a structure stays.

## Environment

Use Julia 1.12. Clone the repository and run the tests from the package directory:

```julia
using Pkg
Pkg.test()
```

Set `LINECABLEMODELS_FEM=true` to run the finite-element cases as well. They download GetDP
on first use and take about ten minutes.

Build the documentation with `julia --project=docs docs/make.jl`. The build runs the
examples and the code of the guide pages, so it needs Gmsh and the GetDP download as well.

## Checks before a pull request

Every change passes these checks.

- `Pkg.test()`, with every guard of `test/architecture/` at zero. A guard failure is a
  finding about the change, never about the guard.
- The formatter, idempotent: `using JuliaFormatter; format(".")` in the SciML style.
- Vale 3.23.0 and cspell 10.3.6 on every maintained file, as `.github/prose/files.sh` lists
  them, with zero errors and zero unknown words. A finding is fixed by rewriting the
  sentence, never by a disable directive, and a word joins the dictionary only by a decision of
  the maintainers.
- gitlint on the commit message: `type(scope): description`, an imperative subject of at
  most 72 characters. Install gitlint and its hook once per clone:

```sh
uv tool install gitlint-core
gitlint install-hook
```

## Pull requests

Keep a pull request to one change. State what it adds, the measurements it rests on and the
names it introduces with their meanings. Add tests that check feature results. Update the
documentation when the public API changes, under the [conventions](docs/src/developers/conventions.md).

Continuous integration runs the whole suite on Linux, macOS and Windows, uploads coverage to
Codecov, and runs the writing checks on the sources and the rendered documentation. The
finite-element cases run on pull requests that are not drafts. Coverage of `src/` and `ext/`
stays at or above 95 % on every pull request, through tests that check feature results.

This project follows the [code of conduct](CODE_OF_CONDUCT.md).
