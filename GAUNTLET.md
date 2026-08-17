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
| Detailed full-frequency cases | 689 |
| Sparse successful cases | 180 |
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
| Python harvester | source/configuration only | none | `gauntlet/harvest/pscad` | verified slim boundary |
| Python caches and generated runs | remove | none | none | excluded from maintained harvester |

The original ignored dataset was deleted only after both staged archives, the
tracked smoke dataset, two deterministic normalization runs, extraction tests,
hash verification, and full-dataset accounting succeeded.

Final local artifact identities:

| Artifact | Tree SHA-256 | Archive SHA-256 | Files | Uncompressed bytes |
|---|---|---|---:|---:|
| `pscad-raw-v1.tar.gz` | `55bb88973ee3bbb183cea1342d45ea3da5493123a968232477a25af3e21d1140` | `cb78519befb85de37ed1773c5d597f40694294d4bb077124ded515f27f4eb7cd` | 13,821 | 1,151,403,828 |
| `pscad-normalized-v1.tar.gz` | `957397f155673fd54ac3fd361a6cc60815eed01e6d1fa29e8ac62bb1ae20f6fd` | `05d173f623de420cd34e98fd56c8767aa75f9dd7feb5c049dc991c213b57896c` | 1,071 | 321,801,296 |

Independent `ingest-v9` and `ingest-v10` runs produced identical raw and
normalized tree hashes. Both archives were extracted into temporary storage,
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
| Fidelity | `Exact`, `Approximate`, or `Rejected` |

Detailed cases consume every available phase, sequence, modal, calculated,
and fitted channel. Sparse cases compare ordinary 60 Hz Z/Y and sequence data
and evaluate fitted Yc/H over their available range. Full-frequency Z/Y is not
inferred from sparse fitted propagation data.

Final normalized channel counts are 869 phase references, 568 sequence
references, 689 modal references, 868 vector fits, and 72 terminal-response
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
| Pipe type | closest valid native common-enclosure/coaxial representation | diagnostic approximation | 72/72 materialized |

Saad, Wedepohl, Deri–Semlyen, Ametani, and Lucca names are never presented as
implemented LineCableModels formulations unless genuine numerical equivalence
is demonstrated. Every approximation is represented by an `Assumption` and is
excluded from pass/fail tolerance calibration.

No PSCAD case is currently promoted to `Exact`. This is deliberate: the live
comparisons show useful agreement but do not yet demonstrate complete geometry
and formulation equivalence. Materialization, computation, port alignment, and
comparison exceptions remain hard failures regardless of fidelity; only a
successfully evaluated numerical discrepancy may be diagnostic.

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

The source-internal vector-fit calibration produced two frozen classes within
the `1e-3` cap:

| Quantity/family | Calibration cases | Held-out cases | Calibration maximum | Held-out maximum | Frozen relative tolerance |
|---|---:|---:|---:|---:|---:|
| Yc / overhead | 146 | 22 | `2.46178e-5` | `9.52933e-6` | `5e-5` |
| Yc / pipe | 28 | 8 | `1.80072e-4` | `1.81587e-4` | `5e-4` |

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
| Tutorial 3 EMT | 17.47/19.03 ms min/median; 2,098,496 bytes; 14,056 allocations | 17.26/18.63 ms; 2,098,496 bytes; 14,056 allocations | pass |
| Overhead | 108.06/111.91 ms; 11,304,744 bytes; 57,485 allocations | 96.78/100.32 ms; 11,306,952 bytes; 57,485 allocations | pass |
| Mixed | 75.38/80.99 ms; 6,639,136 bytes; 69,814 allocations | 67.26/71.63 ms; 6,639,136 bytes; 69,814 allocations | pass |
| Pipe | 30.95/32.01 ms; 3,802,096 bytes; 36,761 allocations | 27.01/28.59 ms; 3,802,128 bytes; 36,761 allocations | pass |
| Coaxial, retained | 21.49/22.29 ms; 2,895,624 bytes; 26,524 allocations | 19.38/20.19 ms; 2,895,688 bytes; 26,524 allocations | pass |
| Coaxial, no reduction | 21.44/22.18 ms; 2,895,624 bytes; 26,524 allocations | 19.93/20.53 ms; 2,895,656 bytes; 26,524 allocations | pass |

Every representative allocation count is unchanged; byte ratios range from
`1.0` to `1.00020`, within the frozen five-percent gate. Timings are trends and
are not used as correctness gates.

## Verification ledger

| Gate | Evidence | Status |
|---|---|---|
| Datasource alignment | `Datasource`, `Dataset`, symbol dispatch, external fixture datasource, PSCAD `ingest`/`load`/`decode`, FEM not-implemented paths, and generic-loader isolation | 24/24 focused assertions passed |
| Datasource-aligned full dataset | all 869 PSCAD cases and 199 rejections loaded through the new dataset dispatch; normalized archive extraction reproduced `957397f155673fd54ac3fd361a6cc60815eed01e6d1fa29e8ac62bb1ae20f6fd` | verified |
| Root package baseline | 2,880/2,880 tests passed from baseline SHA | verified |
| Parser unit tests | input, ordinary, detailed, terminal, fit, malformed/truncated/checksum, modal assignment, port order, Kron, reports | all application tests passed |
| Sixteen-case tracked smoke suite | fixed IDs in `gauntlet/config/smoke.toml`; datasource-aligned tree `b696494f6ffc51274ec541e5c85349e4480a76d18cd32c41b9f3a75a28d82355` | passed before and after source-dataset deletion; zero hard failures |
| 869 successful canonical materializations | performed during each ingestion pass | verified |
| 869 successful canonical computations | full extracted-artifact run | verified |
| 199 rejection records | indexed and lazily loadable | verified |
| Detailed/sparse channel accounting | 689 detailed; 180 sparse; per-channel counts recorded above | verified |
| Deterministic normalization | identical v9/v10 raw and normalized trees | verified |
| Raw archive extraction and digest | tree and archive hashes above | verified |
| Normalized archive extraction and digest | tree and archive hashes above; `Dataset` opened | verified |
| Full local artifact suite | 1,068 canonical trials: 869 computed successes and 199 source rejections; 10,715 diagnostic, 1,022 pass, 360 unavailable, 199 source-rejected comparison rows, zero fail rows | passed |
| Allocation regression suite | six warmed cases; byte and allocation ratios below 1.05 | passed |
| Gauntlet application tests | parser, comparison, smoke, reports, archive identity, opt-in full accounting | passed |
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

`Report(smoke, 11/16 acceptable trials)` and
`Report(full, 888/1068 acceptable trials)` are not failure counts. Every native
execution succeeded. The lower summary count reflects intentionally
non-gating `Approximate`, `Unavailable`, and `ReferenceRejected` evidence; the
CLI exits nonzero for any `Fail` verdict, and both runs exited zero.

## Delivery restrictions

Work remains on `release/v0.2.0`. Gauntlet does not upload artifacts, create a
release, push, merge, rebase, force-push, tag, register the package, or mutate
`main`. A real downloadable `Artifacts.toml` binding is added only after the
user provides the uploaded artifact URLs.
