# Implementation verification

The [implementation table](assimilation.tsv) accounts for all 144 formula
records: 54 implemented records or approximation variants, 42 verified
existing records, 26 equivalent witnesses, and 22 deferred records. There
are no pending candidate records. These are record counts, not counts of
distinct publications or registry identifiers.

The [implementation table](assimilation.tsv) gives each record's scope and processing state.
[Formula equivalence and selection](deduplication.md) records the retained
evaluators, compatibility aliases, approximation choices, and restrictions.
Deferred formulations require field discretization, an independent modal
solution, unsupported geometry, or constitutive inputs outside the current
coaxial representation. They have not been substituted with unrelated
analytical formulas.

## Checks performed

| Check | Result |
| --- | --- |
| Default package suite, before the final validation-ownership correction | 23,402 passed, two validation-architecture failures, and two skipped GetDP tests. All numerical formula tests passed. |
| Follow-up after the ownership correction | 2,116 passed: 113 validation-architecture checks, 198 common-pipe checks, and 1,805 complete survey-to-registry checks. Both earlier failures are resolved. |
| Package quality suite | 3,662 passed. |
| Benchmark toolkit, including PSCAD mappings | 751 passed. No external PSCAD simulation was launched. |
| Earth-formula benchmark catalogue | All 59 canonical registrations enumerated; source paths and hashes resolved; six catalogue assertions passed. |
| Survey citations | All 112 cited keys resolve against the 123-record bibliography. |
| Documentation | Documenter completed its doctests, citation expansion, cross-reference checks, exported-docstring checks, and HTML rendering. |
| Patch whitespace | `git diff --check` passed. |

The final ownership correction moved common-pipe checks directly into the
`CableBlueprint` validator and replaced temporary mutation with direct
checks. It changed neither formula evaluation nor matrix assembly. The
affected tests were rerun; the complete default suite was not repeated
after this final correction.

The two skipped tests require an installed GetDP executable for a real
multifrequency solve and frozen numerical comparisons. They are not counted
as passes. Numerical checks here do not constitute a new external FEM
benchmark campaign or an independent reproduction of every publication's
figures.

The documentation build used an isolated repository copy because the build
regenerates tutorial pages, case pages, and plot assets. Its stale gauntlet
dependency manifest was resolved in that copy only. Following the
API-reference correction, Documenter reused the successfully generated
pages and assets. The remaining warnings concern PNG fallbacks for large
example outputs and the search-index warning threshold; neither prevented
rendering.

## Reproduction

Run the numerical suite and quality checks in their test environment:

```sh
julia --project=test --startup-file=no test/runtests.jl
julia --project=test --startup-file=no test/runtests.jl tag:quality
julia --project=test --startup-file=no test/runtests.jl input_validation_architecture.jl literature_common_pipe.jl literature_metadata.jl
```

The benchmark toolkit uses its separate environment:

```sh
julia --project=test/gauntlet --startup-file=no test/runtests.jl tag:gauntlet_toolkit
```

Instantiate the documentation and gauntlet environments before building
the manual. Use a disposable working copy if the generated documentation
files should not change in the integration tree:

```sh
julia --project=docs --startup-file=no docs/make.jl
```

## Integration state

Changes remain uncommitted on `release/v0.2.0-literature-survey`. No merge
or publication was performed. Nine redundant or misattributed formula
files were removed from the working tree; their tracked versions remain
recoverable from Git, and their former selections have explicit
compatibility mappings described in the equivalence document.
