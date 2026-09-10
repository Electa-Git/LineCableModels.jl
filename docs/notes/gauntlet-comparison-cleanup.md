**Gauntlet comparison cleanup — 2026-09-09**

This profile applies the native-first semantic-economy and dispatch-template playbooks to `release/v0.2.0-pscad-ext`. Its starting point is HEAD `2d694a2444f52942e9d14cd5ff835760295ff532` plus the uncommitted worktree reviewed in the compliance audit. The pre-edit source snapshot and file hashes were retained in `/tmp/lcm-audit-fix-before`; moved files are not counted as new modules or new concepts.

The calculation definitions remain authoritative. Physical models, formula tags and indexed self/mutual methods, numerical types, frequency samples, material laws and raw solver results are preserved. LCM, PSCAD and FEM remain independent calculations; the reference determines comparison direction and normalization. Numerical-zero traces retain their absolute differences and the explanation for unavailable relative RMS.

`BenchmarkDefinition.comparison_settings` is a NamedTuple. It contains `quantities`, `statistics`, `bands`, `normalizations`, `atol`, `fundamental`, `harmonics` and `unsupported`. Scalar line-parameter comparisons default to Z/Y, value statistics, all stored frequencies and reference-RMS normalization. Moment comparisons explicitly select R/L/C/G means and standard deviations; their currently supported full-band reference-RMS calculation is retained. Unsupported moment settings are rejected rather than ignored. `LineParametersPolicy`, `UQMomentPolicy` and `benchmark_calculation` are removed without resolving aliases.

`run_benchmark` validates its declaration, computes or recovers each operand, obtains the result representation through dispatch, calculates the selected comparisons once, and returns that comparison. Saving retains those same error objects and settings. `compare_saved` is the separately requested calculation over completed operands; it uses the same comparison methods and the same writer. Neither operation constructs an additional comparison with default settings.

| Construct or method family | Owner and meaning | Admission and callers |
|---|---|---|
| `record_benchmark` | Gauntlet writes completed comparisons and their operand checksums atomically; it verifies an existing record before reuse. | The sole new generic. It separates an existing file-writing algorithm from RMS calculation and is used by both `run_benchmark` and `compare_saved`. It does not create another result type or execution sequence. |
| `validate(compare; ...)` | Engine checks RMS controls without accessing calculated values. | Method on existing `validate`; reused by direct comparison and benchmark preflight. Rejects invalid band selectors, bounds, tolerances, normalization and unsupported-quantity explanations. |
| `validate(::BenchmarkDefinition)` | Gauntlet checks declared comparison settings, input coordinates and timing controls before calculation or file creation. | Method on existing `validate`; used by individual, campaign and saved-result actions. Result-dependent checks remain in `compare`. |
| Formula `NamedTuple` conversions | Each formula owner exposes its actual equation selection, parameters, hooks and options; earth formulas also expose assumptions and explicit equivalent-earth reduction. | Native constructors replace structural probing. Callables and numerical values retain their identity. Missing required conversion methods fail; they do not produce empty records. |
| PSCAD `NamedTuple` conversion | PSCAD exposes the requested/native selections and their technical descriptions. | Replaces its separate `formulation_record` generic and Gauntlet's unchanged cross-module forward. Native setting values are unchanged. |
| Parametric-result `NamedTuple` conversion | ParametricBuilder exposes its formulation, ordered values, axes and details. | Native conversion supplies the existing result record without copying arrays. Consumers no longer reach into its fields to read axes. |
| Result `semantic_sha256` methods | Gauntlet defines the numerical checksum representation once for scalar results, parametric results and moments. | Methods on the existing hash operation, shared by saving, reading and reporting. The scalar and parametric checksum encoding remains compatible with existing saved calculations. |
| `validate(read_calculation, result, metadata)` and result `select` methods | Gauntlet verifies saved numerical values and selects a declared point for reporting/plotting. | Methods on existing generics replace the report's central result-type switch. Scalar data are read through `observe`, `frequencies`, `basis` and `domain`. |
| Public result accessors | Engine selects the requested Z/Y/R/X/L/G/B/C observation. | Existing public names and results remain supported. Their methods now select an observable explicitly rather than forwarding unchanged to another convenience function. |

`repository_revision` names the recorded commit and dirty state. `records.jl` owns calculation/formula/source records. `_write_value` names the deterministic byte encoding used by the existing checksum operation; its encoding is unchanged. These are renamings, not added responsibilities. The one-caller `_benchmark_metadata` repacker is removed; the benchmark action constructs its retained metadata directly. Source evidence, resolved geometry, parameters and hooks remain retained.

The old on-disk `comparison_policy` key is decoded only by saved-file readers. Its presence in checksummed retained benchmarks establishes the compatibility requirement. Readers convert it to `comparison_settings` in memory; active APIs and new files do not emit the old key. Selection dictionaries and scalar labels are likewise converted to a named record at the reader, without discarding their contents. Plot labels use the explicit calculation IDs retained with each analysis.

The semantic-economy check now includes Gauntlet and checks public forwarding functions. No function-name suppression was added. Native interface methods and methods selecting another dispatch axis remain valid. The alternative-result test demonstrates that reporting and numerical checks do not depend on `LineParameters` field layout. Invalid-request tests observe calculation counts and directory creation; the one-comparison test observes actual calls and compares returned, persisted and tabulated values.

The native PSCAD log now says the compile call returned. The measured duration continues to describe the compile call, separately from the output-readiness wait and whole-call wall time. It no longer announces that the calculation has completed before output verification.

The inventory loads LineCableModels and Gauntlet under Julia 1.12.7 in both the
pre-edit snapshot and the completed worktree. It traverses their owned modules,
counts named generic/type bindings separately, and counts authored methods,
including methods extending Base/Core operations. Imported dependency names and
anonymous compiler-generated names are excluded. This is a binding inventory;
the method count is not a test count.

| Measured item | Before | After | Change |
|---|---:|---:|---:|
| Owned generic bindings | 882 | 879 | −3 |
| Owned type bindings | 180 | 178 | −2 |
| Owned modules | 32 | 32 | 0 |
| Authored methods | 3221 | 3239 | +18 |

The added bindings are `record_benchmark` and the two renamings
`_write_value` and `repository_revision`. Removed bindings are
`PSCAD.formulation_record`, `_benchmark_metadata`, `_canonical_write`,
`benchmark_calculation`, `comparison_policy_record` and `repository_provenance`.
The two removed types are the former comparison-setting wrappers. No registry,
alias, module, execution context or formula-selection mechanism is added.
The existing execution and saved-reanalysis entry points remain; each now
calculates its requested RMS once, with one shared persistence operation.
The earlier execution path calculated RMS twice, or three times when saving.

Gauntlet now participates in ExplicitImports checks for missing/stale imports,
true ownership and qualified access. Its public-access annotation check records
specific native APIs used for instrumentation detection, artifact hashes,
restoring recorded packages and BLAS measurement. Those Base/Pkg/BLAS APIs lack
public annotations; they are not package-owned compatibility wrappers. There is
no exception for private access to a LineCableModels owner and no suppression
of the forwarding detector.

The inventory script and TOML records are retained at
`/tmp/lcm-definition-inventory.jl`, `/tmp/lcm-inventory-before.toml` and
`/tmp/lcm-inventory-after.toml`. The before/after measurements refer to the
pre-existing dirty worktree, so the preceding PSCAD migration is not falsely
counted as newly introduced architecture.

Validation completed on this worktree:

| Check | Result | Retained log |
|---|---|---|
| Ordinary unit/integration selection, including numerical regressions and local PSCAD/FEM execution fixtures | 11,249 passed | `/tmp/lcm-fix-package-tests.log` |
| Quality selection, including ambiguity, formula ownership, semantic economy, imports and Gauntlet | 1,102 passed | `/tmp/lcm-fix-quality.log` |
| Gauntlet benchmark, recovery, UQ, result-space and reporting tests | 871 passed across the broader run and final corrected-file rerun | `/tmp/lcm-fix-gauntlet.log`, `/tmp/lcm-fix-gauntlet-final-items.log` |
| Native Makie figures, matrix overlays, controls and export | 176 passed | `/tmp/lcm-fix-native-makie.log` |
| Actual retained benchmark through the Gauntlet plotting overload | Four matrix cells, two sources, SVG export and unchanged report values verified | `/tmp/lcm-fix-benchmark-overlay.log` |
| Package doctests | Passed | `/tmp/lcm-fix-doctest.log` |
| Documentation build, executable examples and generated case catalogue | Passed; Documenter reports PNG fallbacks for large HTML examples and a search-index size warning | `/tmp/lcm-fix-docs.log` |
| Whitespace and retired-name checks | `git diff --check` passed; no active old API bindings or `.cov` files outside retained work artifacts | Source inspection |

The Gauntlet total counts each final assertion once: 645 checks in unchanged
passing files from the broader run, plus 226 checks in the final
`comparison_execution_tests.jl` and `toolkit_tests.jl` rerun. The broader run had
stopped at a stale scalar fixture supplying acceptance limits that the previous
runner silently ignored. That fixture now requests comparisons only; a separate
behavioral test verifies that unsupported acceptance limits fail with zero
calculations. The added selected-band test returns RMS 1.0 for 10–100 Hz,
observes exactly one comparison call, and verifies identical returned, saved and
tabulated values. These are execution tests, not literature or Markdown assertions.

Reproduction commands from the worktree root:

```sh
julia --project=test --startup-file=no test/runtests.jl
julia --project=test --startup-file=no -e 'using TestItemRunner, LineCableModels; TestItemRunner.run_tests(pkgdir(LineCableModels); filter=item->:quality in item.tags, verbose=true, failfast=true)'
julia --project=gauntlet --startup-file=no test/gauntlet/runtests.jl
julia --project=gauntlet --startup-file=no -e 'using TestItemRunner, LineCableModels; TestItemRunner.run_tests(pkgdir(LineCableModels); filter=item->endswith(item.filename,"toolkit_tests.jl") || endswith(item.filename,"comparison_execution_tests.jl"), verbose=true, failfast=true)'
LINECABLEMODELS_TEST_PLOTTING=true julia --project=test/visual --startup-file=no -e 'using TestItemRunner, LineCableModels; TestItemRunner.run_tests(pkgdir(LineCableModels); filter=item->endswith(item.filename,"native_makie.jl"), verbose=true, failfast=true)'
julia --project=docs --startup-file=no docs/doctest.jl
julia --project=docs --startup-file=no docs/make.jl
git diff --check
```

No live PSCAD station calculation was repeated for this cleanup. The native
formula settings and equations are unchanged; PSCAD-specific edits expose its
record through `NamedTuple`, update callers and correct the compile-call progress
message. All pre-edit scientific formula implementation files under
`src/**/formulas/` retain their bytes. The numerical checks demonstrate preserved
software behavior; they do not designate any backend as physical truth.

No commit or push was made.
