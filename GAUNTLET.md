# Gauntlet implementation ledger

This file is the persistent source of truth for the PSCAD validation dataset and
Gauntlet harness. Inventory rows are updated with evidence as work proceeds;
completed evidence is not silently discarded.

## Baseline

| Item | Value | Status |
|---|---|---|
| Branch | `release/v0.2.0` | verified |
| Commit | `475d953f5c3cc57e0ee76d3a2374f24ba06189d3` | verified |
| Worktree | clean before this file was created | verified |
| Root package tests | 2,880/2,880 passing in 316.9 s | verified |
| Tutorial 3 performance | 17.47 ms minimum, 19.03 ms median, 2,098,496 bytes, 14,056 allocations | verified |
| Representative solver performance | 61-frequency, two-cable EMT: 17.50 ms median, 17.12 ms minimum, 2,099,104 bytes | verified |

Gauntlet is a separate Julia application. It may depend on LineCableModels but
must not add dependencies, exports, or validation-only entry points to the root
package.

## Datasource alignment checkpoint

The first Gauntlet implementation completed the PSCAD campaign but left its
dataset loader and native reconstruction path coupled to PSCAD records. This
checkpoint aligns the owned API before any artifact is published:

- `Datasource` is the abstract origin of external benchmark evidence;
- `PSCAD` is the first complete datasource and owns its parser, ingestion,
  decoding, loading, reconstruction, and calibration specializations;
- `FEM` is registered as a datasource provision and fails explicitly as not
  implemented until a finite-element format is selected;
- `Dataset` is the source-neutral lazy collection of cases and rejections;
- `ingest`, `load`, and `decode` are the short dispatch verbs;
- datasource symbols are registered through `Val`, with `:pscad` supported
  end to end and `:fem` reserved;
- the normalized dataset index declares its datasource explicitly;
- the previous source-type and collection names, including the old suite field,
  are removed without aliases because no Gauntlet API or artifact has been
  released.

The comparison, modal-alignment, tolerance, reporting, and performance layers
remain datasource-neutral. New datasource adapters must normalize evidence into
the existing typed `Reference` surface; source-specific records and file-format
rules remain below their datasource directory.

## PSCAD capture repair checkpoint

This amendment began from clean branch checkpoint `11bba25f` after the initial
Gauntlet implementation was committed.

The original campaign index classified 180 successful cases as sparse because
the harvester watched PSCAD's temporary directory while those detailed files
were emitted into the process working directory. The solver configuration was
not sparse: `Output=YES` was already set. The repaired runner captured the
missing quartet from both output locations.

The ignored amendment at `runs/pendingcases` contains 180 identity-matched
records: 72 coaxial, 72 overhead, and 36 pipe-type. Gauntlet independently
verified 720/720 `_zm`, `_zp`, `_ym`, and `_yp` files, square matrix dimensions,
artifact hashes, and the shared 101-point range from `1e-3` to `1e7` Hz. The
Julia harvester observes PSCAD's temporary directory, the process working
directory, and the generated project directory as one completion set, so this
output-location split cannot force a false timeout. The
repair report SHA-256 is
`2a9c01631e951e6b52804416ff1628f881b9bda932c558caa773995d6532873c`;
the source verification record is
`ade67569822c6cd5431950081ddc466aa5e46de81ac59183fa5850f705068fc0`.

The amended dataset therefore contains 869 detailed successes and zero sparse
successes. All 72 pipe references are retained with `Deferred` fidelity. They
load and remain searchable, but have no native problem/formulation and are not
computed or included in performance or tolerance calibration until native pipe
physics exists.

## Original dataset

The ignored source dataset was inspected at
`benchmarks/PSCAD-reference-v1`. Its authoritative index is
`dataset_manifest.json`, schema version 1, SHA-256
`fc5c78cecf673ecd5272ef992181bcca8c8a853602db534d43815d1e36931b65`.

| Item | Count |
|---|---:|
| Declared records | 1,104 |
| Unique effective inputs | 1,068 |
| Duplicate aliases | 36 |
| Successful records | 905 |
| PSCAD solver rejections | 199 |
| Excluded non-frequency-dependent definitions | 1 |
| Unique successful cases after alias removal | 869 |
| Detailed full-frequency cases in original capture | 689 |
| Missing detailed captures repaired by amendment | 180 |
| Detailed full-frequency cases after amendment | 869 |
| Sparse successful cases after amendment | 0 |
| Coaxial cases | 288 |
| Overhead cases | 240 |
| Mixed overhead/cable cases | 269 |
| Pipe-type cases | 72 |

The source tree occupies approximately 2.1 GiB. The dominant raw classes are
18,628 `.out` files (1,206,127,086 bytes), 1,410 `.clo` files
(336,688,996 bytes), 2,246 `.bakx` files (287,389,408 bytes), 1,164 `.pscx`
files (184,969,730 bytes), and 499 `.tlo` files (43,731,712 bytes). The
`runs/` tree accounts for approximately 2.0 GiB, `diagnostics/` for 83 MiB,
and `donors/` for 27 MiB.

Configuration evidence:

| File | SHA-256 |
|---|---|
| `benchmark.yml` | `1b28f25089ac47152be978ca94f0156e77905dda8adee72887f3f7d671717535` |
| `official_sources.yml` | `700b8017e9de02476fa146b3b707adca0680e01f9a8db0de5e07b242410fa593` |

Disposable parser evidence:

| File | SHA-256 | Disposition |
|---|---|---|
| `import_pul_ZY.m` | `e1a0835e6290270b1b26cb3bb8b1acae6a733d543706ab3bac651c6f07448cb1` | port useful parsing rules, then delete |
| `import_Yc_H.m` | `1de66ab5ec2fb9107c335e385d6336dce3c27d64cb0ef12bd3f15337392f5209` | port useful parsing rules, then delete |
| `importFitted_Yc_H.m` | `f462d9c38cbc8761622cf53ae909f5bcb10ac8bd3a395ad7601a19ee331ba2b2` | port useful fit evaluation, then delete |
| `importFitted_1ph_Yc_H.m` | `fd9b41949f5dbd2ba341370845a5e7bdedb73fce7b4399ef466ce6d64dd7f757` | port useful one-phase fit evaluation, then delete |
| `tlo_clo_files.pdf` | `ad098f1fc80f240499da8f780e94cdfb97d65708066cf351f4f12eadcd161526` | retain name/hash provenance only, then delete |

## Artifact disposition

| Artifact class | Raw artifact | Normalized artifact | Maintained tree | Status |
|---|---|---|---|---|
| Authoritative dataset index | retain | normalized index | counts and digest in this ledger | verified |
| Canonical `.cli`/`.tli` inputs | retain once per effective case | typed source records | representative parser fixtures only | verified |
| Required ordinary/detailed `.out` | retain once | primitive arrays and explicit fields | representative parser fixtures only | verified |
| `.clo`/`.tlo` fit files | retain once | poles, residues, constants, delays, range | representative parser fixtures only | verified; 21 exact-prefix recoveries |
| Successful alias records | alias map only | alias map | one exercised smoke alias | verified |
| PSCAD solver rejections | concise terminal evidence | 199 rejection records | one smoke rejection | verified |
| Excluded definition | concise evidence | one exclusion record | index evidence | verified |
| PSCAD projects and `.bakx` files | remove | none | none | excluded from staged artifacts |
| Donor downloads | remove after provenance capture | none | pinned source descriptors only | excluded from staged artifacts |
| Routine logs and stale diagnostics | remove | none | none | excluded from staged artifacts |
| MATLAB and parser PDF | remove after verified Julia port | none | hashes in this ledger | port verified; disposable source deleted |
| Julia harvester | source/configuration only | none | `gauntlet/harvest/pscad` | owned pure-Julia app |
| Python source, caches, and generated runs | remove | none | none | no project-owned Python retained |

The original ignored dataset was deleted only after both staged archives, the
tracked smoke dataset, two deterministic normalization runs, extraction tests,
hash verification, and full-dataset accounting succeeded.

The maintained harvester is a Julia Pkg app named `linecablebenchmark`. Static
commands (`pending`, `verify`, `ingest`, `catalog`, and `sources`) load no Python
runtime. Live `inspect`, `case`, and `batch` commands enter a pinned, isolated
environment and activate `GauntletPythonCallExt`, whose only job is to call the
external `mhi.pscad` automation API. Configuration is TOML and all planning,
identity, hashing, capture, verification, and staging code is Julia. No
compatibility launcher or project-owned `.py` file remains.

Final local artifact identities:

| Artifact | Tree SHA-256 | Archive SHA-256 | Files | Uncompressed bytes |
|---|---|---|---:|---:|
| `pscad-raw-v1.tar.gz` | `cd72caef3905cd6847b47f4d20948d787f779cb393b43443bedd6abbb3255a07` | `475ee211cda031d09a0882f5d9bf18848456757c7bbf0bba45a20f75e429eb43` | 16,485 | 1,306,876,339 |
| `pscad-normalized-v1.tar.gz` | `4f3b195b9de3baf80d12089bb0f369b3a2d658d0fa73fa20bdcd7b0bcaedc794` | `f3516d641bf4e6810cc5ac49bc569b827f8dfa73afd72d25a1f1578b01d61294` | 1,071 | 349,208,379 |

Independent final amendment ingestions produced identical raw and normalized
tree hashes. Both archives were extracted into temporary storage,
rehash-verified, and the extracted normalized index was opened as a `Dataset`.
The datasource-aligned normalized archive was rebuilt from that verified tree;
all 869 cases and 199 rejections were loaded again, and archive extraction
reproduced the aligned tree hash above.

## Normalized schema and units

The normalized artifact has an explicit version and contains only primitive
arrays and fields. JLD2 is a container, not an opaque Julia object snapshot.

| Quantity | Canonical representation |
|---|---|
| Frequency | Hz |
| Series impedance | Ω/m |
| Shunt admittance | S/m |
| Angles | radians |
| Port ordering | explicit identifiers stored with every matrix family |
| Missing output | absent field; never a fabricated array |
| Reduction | `Retained`, `Reduced`, or `NoReduction` |
| Fidelity | `Exact`, `Approximate`, `Deferred`, or `Rejected` |

Detailed cases consume every available phase, sequence, modal, calculated,
and fitted channel. The amended dataset has no sparse successes; missing
channels remain absent rather than inferred from fitted propagation data.

Final normalized channel counts are 869 phase references, 568 sequence
references, 869 modal references, 868 vector fits, and 96 terminal-response
records. Case
`de4f15dc8a3dbe15e94065c7c6ff7d81b098d480f5cc643fb4f5bd5499cee4d3`
contains the single source fit file with no coefficient payload and therefore
records the fit as absent.

## Reconstruction and formulation ledger

| Family | Native materialization | Initial fidelity rule | Status |
|---|---|---|---|
| Coaxial | existing materials, cable DSL, positions, Earth, problem, formulation | approximate pending formulation equivalence | 288/288 materialized |
| Overhead | existing materialized conductor/system/problem path | approximate pending earth formulation equivalence | 240/240 materialized |
| Mixed | existing heterogeneous materialized system path | diagnostic approximation | 269/269 materialized |
| Pipe type | no substitute materialization | deferred until native pipe physics exists | 72/72 references retained |

Saad, Wedepohl, Deri–Semlyen, Ametani, and Lucca names are never presented as
implemented LineCableModels formulations unless genuine numerical equivalence
is demonstrated. Every approximation is represented by an `Assumption` and is
excluded from pass/fail tolerance calibration.

No PSCAD case is currently promoted to `Exact`. This is deliberate: the live
comparisons show useful agreement but do not yet demonstrate complete geometry
and formulation equivalence. Materialization, computation, port alignment, and
comparison exceptions remain hard failures for implemented families; only a
successfully evaluated numerical discrepancy may be diagnostic. A `Deferred`
pipe trial reports one `Unavailable` implementation check without invoking
`compute!`.

## Numerical comparison

Every matrix comparison records maximum absolute error, scaled maximum
relative error, Frobenius relative error by frequency, worst frequency, worst
port pair, symmetry, reciprocity, positive-real/passivity diagnostics, and
alignment failures. Paired cases additionally record Kron consistency.

Modal comparison performs one-to-one correlation matching, phase alignment,
degenerate-subspace comparison, modal-value comparison after matching, and
phase-domain reconstruction. Raw PSCAD mode order, sign, or complex phase is
not treated as a physical discrepancy. Native modal evidence is decomposed
directly at each frequency for comparison; Gauntlet does not alter or bypass
the production `compute!` result. Sequence transforms operate independently on
each recorded three-phase circuit and preserve retained ground-wire channels.

Tolerance calibration is deterministic: four SHA buckets calibrate and one
bucket validates; thresholds use source-resolution absolute floors, a
factor-of-two margin, 1–2–5 rounding, and a maximum gating relative tolerance
of `1e-3`. Thresholds will be frozen in `gauntlet/config/tolerances.toml` and
will never recalibrate during ordinary validation.

The source-internal vector-fit calibration produced one frozen implemented
class within the `1e-3` cap:

| Quantity/family | Calibration cases | Held-out cases | Calibration maximum | Held-out maximum | Frozen relative tolerance |
|---|---:|---:|---:|---:|---:|
| Yc / overhead | 146 | 22 | `2.46178e-5` | `9.52933e-6` | `5e-5` |

Pipe calibration is intentionally absent while the family is deferred.
Coaxial and mixed Yc coefficients and every H family contain source-internal
outliers above the cap. They remain explicit diagnostics rather than causing
tolerance inflation. Ordinary 60 Hz output is not itself a detailed-grid
sample; a declared polar interpolation on log frequency cross-check covers 552
cases (137 lack a usable ordinary matrix), with median/max residuals of
`5.83057e-4`/`1.82367e-3` for Z and `4.73084e-7`/`1.67454e-5` for Y.

## Performance

The pinned Julia 1.12 Linux baseline records warmed minimum/median time, bytes,
allocations, frequency count, conductor count, and points per second. The
allocation and byte gates allow 5% regression. Timing is a reported trend, not
a failing gate. PSCAD elapsed time is provenance and is never used as a direct
speed ratio.

| Case | Baseline | Current | Status |
|---|---|---|---|
| Tutorial 3 EMT | 17.47/19.03 ms min/median; 2,098,496 bytes; 14,056 allocations | 19.39/19.84 ms; 2,098,496 bytes; 14,056 allocations | pass |
| Overhead | 108.06/111.91 ms; 11,304,744 bytes; 57,485 allocations | 99.16/102.99 ms; 11,306,888 bytes; 57,485 allocations | pass |
| Mixed | 75.38/80.99 ms; 6,639,136 bytes; 69,814 allocations | 66.13/71.54 ms; 6,639,136 bytes; 69,814 allocations | pass |
| Coaxial, retained | 21.49/22.29 ms; 2,895,624 bytes; 26,524 allocations | 18.81/19.22 ms; 2,895,656 bytes; 26,524 allocations | pass |
| Coaxial, no reduction | 21.44/22.18 ms; 2,895,624 bytes; 26,524 allocations | 18.08/19.42 ms; 2,895,656 bytes; 26,524 allocations | pass |

Every representative allocation count is unchanged; byte ratios range from
`1.0` to `1.00019`, within the frozen five-percent gate. Timings are trends and
are not used as correctness gates.

## Verification ledger

| Gate | Evidence | Status |
|---|---|---|
| Datasource alignment | `Datasource`, `Dataset`, symbol dispatch, external fixture datasource, PSCAD `ingest`/`load`/`decode`, FEM not-implemented paths, and generic-loader isolation | 24/24 focused assertions passed |
| Datasource-aligned full dataset | all 869 PSCAD cases and 199 rejections loaded through dataset dispatch; normalized archive extraction reproduced `4f3b195b9de3baf80d12089bb0f369b3a2d658d0fa73fa20bdcd7b0bcaedc794` | verified |
| Root package baseline | 2,880/2,880 tests passed from baseline SHA | verified |
| Parser unit tests | input, ordinary, detailed, terminal, fit, malformed/truncated/checksum, modal assignment, port order, Kron, reports | all application tests passed |
| Sixteen-case tracked smoke suite | fixed IDs in `gauntlet/config/smoke.toml`; amended tree `ddbf13f6a85cb6e01b9c745912c896c694cd2a323b070c1053486d8b93344047`; 14 acceptable and two deferred pipe trials, zero failures | passed |
| 797 implemented canonical materializations | performed during each ingestion pass; 72 pipe references deferred | verified |
| 797 implemented canonical computations | full extracted-artifact run; pipe references not computed | verified |
| 199 rejection records | indexed and lazily loadable | verified |
| Detailed/sparse channel accounting | 869 detailed; zero sparse after verified amendment | verified |
| Deterministic normalization | identical v9/v10 raw and normalized trees | verified |
| Raw archive extraction and digest | tree and archive hashes above | verified |
| Normalized archive extraction and digest | tree and archive hashes above; `Dataset` opened | verified |
| Full local artifact suite | 1,068 trials and 12,392 comparisons; 996 acceptable, 72 deferred pipe trials, zero failures | passed |
| Allocation regression suite | five warmed implemented cases; unchanged allocation counts and at most 1.00019 byte ratio; pipe excluded until implemented | passed |
| Gauntlet application tests | parser, comparison, smoke, structured reports, archive identity, and optional full accounting; 115/115 default and 122/122 artifact-enabled assertions | passed |
| Harvester boundary | exact `linecablebenchmark` Pkg app; static command remained Python-free; isolated live command activated PythonCall and reached the external `mhi.pscad` import | passed |
| Strict documentation | Literate generation, doctests, cross-references, document checks, and HTML rendering | passed |
| Root package suite | 2,880/2,880 in 5m07s | passed |
| Core-only boundary | 32/32 | passed |
| Aqua | ambiguities, undefined exports, stale dependencies, compat, piracy, persistent tasks | passed |
| SciML formatting | Gauntlet source, tests, scripts, and Literate guide | passed |
| Cairo visual suite | 458/458, including all existing goldens | passed |
| GL gallery smoke | 18 panels built without native windows | passed |
| Clean installation | empty depot; Makie extension remained unloaded | passed |
| Literate documentation | `docs/literate/gauntlet.jl` uses the tracked typed records; strict doctests and HTML build | passed |
| Original dataset removal | deleted after all archive, extraction, full-suite, and smoke gates; recoverable from the raw archive or the user's backup | complete |

Report summaries are not failure counts. Implemented native executions must
succeed; `Approximate`, `Deferred`, and `ReferenceRejected` evidence remains
explicit. The CLI exits nonzero for any `Fail` verdict.

## Delivery restrictions

Work remains on `release/v0.2.0`. Gauntlet does not upload artifacts, create a
release, push, merge, rebase, force-push, tag, register the package, or mutate
`main`. A real downloadable `Artifacts.toml` binding is added only after the
user provides the uploaded artifact URLs.
