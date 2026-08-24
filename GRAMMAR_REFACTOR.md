# Unified LineCableModels Grammar refactor

## Goal

Reconcile LineCableModels with the package-local Unified CableHosting Grammar
without changing scientific equations, accepted numerical values, released
v0.1 behavior, or the immutable PSCAD benchmark collection.

## Baseline

- Branch: `release/v0.2.0`
- Source commit: `abfa45b056636592c5e919d9a0a39d74feebced4`
- Released compatibility baseline: tag `v0.1.0`, commit
  `26faf942324cbcad5e7a7173de4bdaa464358a4f`
- Gauntlet binding: `gauntlet_pscad_v1_0_0`
- Gauntlet tree: `ec9bd7850a744dad14a9547ad4c0e0555363459a`
- Engine integration fixture: `test/fixtures/reference/datamodel_engine_baseline.json`
- Independent analytical fixture:
  `test/fixtures/reference/coaxial_capacitance.toml`
- Pre-existing user change: `RELEASE_AUDIT.md`. It is excluded from this
  refactor and must remain unstaged and unmodified.
- Sibling repositories are outside this goal. Their written grammar is used
  only as a conformance reference.

## Immutable Gauntlet evidence

The accepted artifact is not recorded, replaced, rebound, or edited during
this refactor. Before executable case files change, their exact baseline bytes
are retained as inert `.jl.source` files under
`test/gauntlet/reference_sources/gauntlet_pscad_v1_0_0/`. Stored case hashes
are checked against those bytes. Refactored executable cases are checked
separately by semantic formulation, problem structure, frequency vector,
terminal order, basis, domain, dimensions, and numerical result.

Cases:

- `benchmark_132kV_630mm2_flathor_pscad`
- `benchmark_18kV_1000mm2_trifoil_pscad`
- `benchmark_380kV_2000mm2_flatver_pscad`
- `benchmark_525kV_1600mm2_bipole_pscad`
- `benchmark_640kV_2000mm2_bipole_pscad`
- `benchmark_solid_1000mm2_single_pscad`
- `benchmark_two_bare_wires_pscad`

## Released compatibility ledger

The method-level audit is performed against `v0.1.0` before a current public
name is removed. Released material, cable, geometry, collection, import,
export, DataFrame, plotting, property-forwarding, basis, domain, and return
behavior is retained. The following retirement paths remain actionable:

- FEM and sector support points to branch `legacy/fem-sector` at commit
  `b75dd2723f90a83ec090b20605ea42af57f4a9c3`.
- Unversioned legacy JSON points to the last migration-capable commit
  `a71bdfe1ac832f27a0c88b1d02596194aac46ec7`.
- Versioned material and cable persistence remains supported.
- Plot backend `force` remains an accepted no-op source-compatibility keyword.

The transitional v0.2 calculation, parameterization, and result names were
confirmed absent from v0.1.0. They were therefore replaced directly, without
compatibility aliases. The released `CableConstants` behavior and Gauntlet
evidence remain protected independently of those ownership changes.

### Method-level v0.1 audit

The released names below were inspected directly at tag `v0.1.0`. None is
removed or newly redirected by this grammar-only change. Where the baseline
at `abfa45b0` had already normalized an interface, the established v0.2
replacement and its tests remain unchanged rather than acquiring a second
adapter in this refactor.

| Released v0.1 surface | Baseline v0.2 behavior retained by this goal | Evidence |
| --- | --- | --- |
| `Material(...)`; `MaterialsLibrary(...)`; `store_materials_library!`, `remove_materials_library!`, `save_materials_library`, `list_materials_library`, `get_material` | Material construction, dictionary operations, `DataFrame`, versioned save/load, and displays remain owned by `Materials` | material unit, library, persistence, and DataFrame tests |
| `EarthLayer`, `EarthModel`, `addto_earth_model!`, `earth_data`; `CPEarth`, `EnforceLayer` | Static earth construction, `add!`, `DataFrame`, display, constant-property evaluation, and enforced-layer dispatch remain available through their established v0.2 owners | EarthProps and Engine formulation tests |
| `WireArray`, `Strip`, `Tubular`, `Semicon`, `Insulator`, `Conductor`, `addto_conductor!` | Materialized circular/rectangular strands, strip, tubular, semiconductor, insulator, conductor-group construction, equivalents, and same-scalar `add!` remain covered; the earlier approved `WireArray` retirement is not reopened here | DataModel constructor and refactor fixture tests |
| `CableComponent`, `NominalData`, `CableDesign`, `addto_design!`, `design_data`, `design_preview` | Materialized design construction, property forwarding, DataFrame views, indexing, eager equivalents, and optional Makie `preview` remain unchanged | DataModel, DataFrame, and plotting tests |
| `CablesLibrary`, `save_cables_library`, `store_cables_library!`, `remove_cables_library!`, `get_design`, `list_cables_library` | Dictionary operations, JLS trusted-input behavior, versioned JSON, DataFrame conversion, and display remain unchanged | persistence and library tests |
| `CableDef`, `LineCableSystem`, `addto_system!`, `system_data`, `system_preview`, `trifoil_formation`, `flat_formation` | Materialized positions/systems, derived counts, formations, DataFrame output, overlap validation, and optional Makie `preview` remain unchanged | system, geometry, DataFrame, and plotting tests |
| `calc_equivalent_alpha`, `calc_parallel_equivalent`, `calc_helical_params`, resistance/GMR/GMD/inductance/dielectric kernels | The previously approved v0.2 mathematical names retain the same equations and fixed numerical fixtures; this grammar refactor changes no kernel | BaseParams tests, `datamodel_engine_baseline.json`, and `coaxial_capacitance.toml` |
| `export_pscad_lcp` | The established v0.2 `export_data(:pscad, ...)` implementation and PSCAD Gauntlet snapshots remain unchanged | exporter tests and immutable Gauntlet suite |
| released `show`, DataFrame, property-forwarding, formation return shapes, physical units | Existing v0.2 tests are retained and expanded only for Grammar identity and observables | root and visual test suites |

The v0.1 utility constants and uncertainty helpers were already superseded at
the audited v0.2 baseline by owner-local constants, `nominal`,
`standard_uncertainty`, Grid uncertainty declarations, and optional extension
methods. This goal does not reintroduce global constants or a `Utils` module.

## Converged ownership ledger

| Owner | Final responsibility | Evidence |
| --- | --- | --- |
| `Grammar` | shared problem, formulation, and result roots; `compute`, `observables`, `primitives`, and `preprocess` generics | owner identity and no-fallback guards |
| `Engine` | analytical formulation, immutable solver input, primitive calculation, traces, and benchmarks | formulation and integration fixtures |
| `ParametricBuilder` | Grid/Gridspace definitions, deterministic traversal, composite results, and `CalculationManifest` | ordering, details-schema, and deterministic-identity tests |
| `UQ` | direct propagation, conditional Monte Carlo, statistics, samples, and histogram densities | replay, covariance, and extension tests |
| `PlotBuilder` | definition-driven `PlotRecipe` construction through one fixed stage sequence | architecture, Cairo, and backend-isolation tests |
| `ImportExport` | explicit format dispatch and path-or-throw write boundaries | persistence and exporter tests |

`primitives` and `preprocess` are protected future provisions. The package
defines no method for either generic, so unsupported calculation orderings fail
by ordinary Julia dispatch. Every composite result retains its details
dictionary and manifest independently of those reserved operations.

## Result and manifest freeze

- `ParametricResult` and `LinearErrorResult` contain `formulation`, primitive
  `values`, matching `space`, and `details::Dict{Symbol,NamedTuple}`.
- `MonteCarloResult.stats` has exactly `mean`, `std`, `min`, `q05`, `median`,
  `q95`, `max`, and `n`.
- Stable detail keys are `:failures`, `:samples`, `:histograms`, `:random`, and
  `:manifest`.
- Composite result primitive parameter `T` is always `CableConstants` or
  `LineParameters`; composite results do not nest.
- Manifest serialization freezes only after numerical parity is established.

## Verification ledger

| Gate | Baseline | Final |
| --- | --- | --- |
| Root package suite | 3,146/3,146 | 3,311/3,311 |
| Seven-case Gauntlet snapshot suite | 65/65 | 65/65; accepted arrays and sources unchanged |
| Gauntlet toolkit tests | 183/183 | 183/183 |
| DataModel/Engine integration fixture | included in root: 68/68 | 68/68 |
| Independent coaxial-capacitance fixture | included in root baseline | passes in Engine formulation tests |
| Gridspace and construction tests | pending | 593/593 focused ParametricBuilder tests |
| Deterministic and uncertainty tests | pending | 26/26 manifest and 152/152 result-container tests |
| Measurements/Distributions extensions | pending | 62/62 optional-adapter tests in the root suite |
| Persistence and compatibility | pending | 43/43 persistence, 225/225 exporters, and 26/26 cable constants |
| Plotting and visual smoke tests | pending | 470/470 visual tests; isolated GLMakie 6/6 and WGLMakie 12/12 activation tests |
| Documentation and doctests | pending | documentation build and 1/1 doctest pass |
| Formatting, Aqua, clean installation | pending | JuliaFormatter clean; Aqua clean; isolated core 32/32; empty-depot install/import passes |
| Source-amended production coverage | pending | 6,757/7,111 lines (95.02%) |

Comparable performance remains within the accepted per-case limits: median
time ratio `1.20`, bytes ratio `1.05`, and allocation ratio `1.05`. Evidence
from a non-comparable environment is advisory only.

The seven-case snapshot suite passed these unchanged limits in one complete
run. An earlier batch contained one transient 380 kV timing outlier; its
immediate isolated rerun and the final complete replay both passed. No
tolerance, accepted result, snapshot, or artifact binding was modified.

The Cairo suite exercised the real renderer and golden images. Dedicated
isolated checks also activated the GLMakie and WGLMakie extensions without
rendering a display-dependent scene. The existing CI and manual graphical
boundaries remain unchanged.

## Final integrity evidence

- `test/gauntlet/Artifacts.toml` remains byte-unmodified and binds
  `gauntlet_pscad_v1_0_0` to
  `ec9bd7850a744dad14a9547ad4c0e0555363459a`.
- All seven archived source files match the `case_sha256` values accepted by
  the immutable snapshot reader.
- The maintained tree contains no removed unreleased names except deliberate
  absence guards and the test-local decoder for the immutable v1 formulation
  metadata record.
- Production sources contain no dependency or reference to CableHosting or
  PowerImpedance.
- A clean installation from an empty depot resolved, precompiled, and imported
  LineCableModels without activating a Makie extension.
- `git diff --check` is clean. The pre-existing user-owned
  `RELEASE_AUDIT.md` modification remains untouched and unstaged.

## Execution status

All nine planned commits are complete. The refactor changes
ownership, naming, dispatch, and result organization only. No scientific
equation, physical constant, accepted numerical value, persistence version,
Gauntlet tolerance, or performance limit changed. No compatibility alias was
introduced for an unpublished name, and no sibling repository was modified.

## Change-control rule

Stop before changing a scientific equation, physical constant, numerical
tolerance, accepted artifact, released behavior, or sibling repository. Such a
change requires separate review and is not hidden in this ownership refactor.
