# DataModel–Engine refactor ledger

This file is the persistent source of truth for the DataModel-to-Engine refactor on
`release/v0.2.0`. Source changes must match an inventory row below. If an unlisted
public behavior or numerical change is discovered, record it here and pause for user
review before implementation.

## Baseline

- Branch: `release/v0.2.0`
- Commit: `a71bdfe1ac832f27a0c88b1d02596194aac46ec7`
- Worktree at capture: clean
- Julia: 1.12
- Fresh default-suite result: 2,962/2,962 passing in 4m12.3s
- Historical count in the accepted plan: 2,623; superseded by the fresh result above
- Numerical fixtures: `test/fixtures/reference/datamodel_engine_baseline.json`,
  captured by executing the immutable baseline snapshot in a temporary checkout;
  includes every eager component/group equivalent in the MV fixture, CableConstants,
  reduced/unreduced results for both legacy `ideal_transposition` option states,
  and all six primitive trace tensors. Fixture keys name the option state because
  the baseline implemented its meaning inconsistently between Z and Y.
- Representative EMT benchmark: three-phase MV fixture, 50/500 Hz, phase domain,
  reduction and temperature correction enabled, trace storage disabled; 11 warm runs,
  median 1.654097 ms and 180,000 bytes

Status values are `pending`, `in progress`, `verified`, and `blocked`.

## Non-negotiable design rules

1. Preserve the eager chain `Material → cable parts → groups → CableDesign →
   LineCableSystem → LineParametersProblem → Formulation → compute!`.
2. Preserve all inventoried capabilities and default numerical results except the
   approved ideal-transposition correction.
3. Material `T0` is reference data. Operating temperature exists only on the line
   problem and is applied locally inside `compute!`.
4. Natural Julia promotion replaces package scalar unions and coercion machinery.
5. Validation is nonmutating and throws native exceptions.
6. Internal mismatches are fixed at the owned type, constructor, protocol, dispatch,
   or module boundary.
7. Work only on `release/v0.2.0`; do not merge, push, rebase, force-push, tag,
   register, or touch `main`.

## Shim and Adapter Rule

A shim is legitimate only when it adapts an external API, framework, tool,
compatibility boundary, or runtime shape that the project does not control.

A legitimate shim or adapter must satisfy most of the following:

1. It adapts an external interface or runtime model.
2. It prevents external framework/tool quirks from leaking into project-owned code.
3. It is narrow and named after the boundary it adapts.
4. It has a stable contract.
5. It does not encode unrelated business rules.
6. It does not silently mutate inputs.
7. It makes testing easier by isolating the external seam.
8. Removing it would force the external quirk into multiple places.

A shim is illegitimate when any of the following are true:

1. It adapts between two internal project shapes that should instead be made coherent.
2. It exists because the project’s own types, constructors, protocols, or dispatch
   rules are weak.
3. It mutates arguments before passing them along to make incompatible internal APIs
   fit.
4. It hides branching that should be explicit dispatch.
5. It preserves a bad internal API in code that is not externally released or
   compatibility-bound.
6. It is vaguely named with prepare, normalize, massage, convert, handle, patch, or
   similar.
7. It makes debugging harder by quietly transforming inputs.
8. It is used to avoid defining the correct value object, protocol, adapter, or
   function-module boundary.

In a project written from scratch, internal shims are suspicious by default.

Do not call lazy glue code a shim.

Use a shim or adapter for external mismatch isolation.

Fix the type, protocol, dispatch, constructor, or boundary when the mismatch is
internal.

Prefer extending an existing concept-owned public surface over creating a new
abstraction. New code must attach to the correct owner. Do not introduce parallel
abstractions to avoid touching the right concept. Do not place neutral shared
contracts inside a runtime boundary.

## Public API migration ledger

| Current API | Result | Status | Evidence |
|---|---|---:|---|
| `set_verbosity!` | removed; use `ComputeOptions(verbosity=...)` | verified | source and migrated tests |
| `validate!` and scattered assertions | root `validate(x)` plus declarative rules | verified | source and migrated tests |
| `MaxFill()` | removed; use `maxfill(Type, ...)` | verified | source and migrated tests |
| `WireArray` | removed | verified | source and migrated tests |
| macro-generated builder parameter types | documented, unexported Specs | verified | source and migrated tests |
| `@gridspace` | retained | verified | source and migrated tests |
| `@relax` | retained only as unused future provision | verified | source and migrated tests |
| operating-temperature cable-part constructors | removed; set temperature on system problem | verified | source and migrated tests |
| `EMTWorkspace` and `_compute_with_workspace` | removed; use `compute!(...; inspect=true)` | verified | source and migrated tests |
| primitive-matrix storage option | removed; use `EMTTrace` | verified | source and migrated tests |
| file logging | removed; scoped console logging only | verified | source and migrated tests |
| unversioned/legacy JSON loading | deprecation notice and throwing tombstone | verified | source and migrated tests |
| maintained JLS/JSON/ATP/PSCAD/TRALIN/XLSX | retained | verified | source and migrated tests |
| `equivalent`, `nonsensify`, `vdeparse` | retained | verified | source and migrated tests |
| PlotBuilder grammar and visuals | retained | verified | source and migrated tests |

## Mathematical rename ledger

All replacements are breaking prerelease renames without aliases. Existing formulas,
mathematical notes, units, and citations remain intact.

| Current | Replacement | Status |
|---|---|---:|
| `calc_equivalent_alpha` | `equivalent_alpha` | verified |
| `calc_parallel_equivalent` | `parallel` | verified |
| `calc_helical_params` | `helix` | verified |
| `calc_strip_resistance` | `strip_resistance` | verified |
| `calc_temperature_correction` | `temperature_factor` | verified |
| `calc_tubular_resistance` | `tubular_resistance` | verified |
| `calc_tubular_inductance` | `tubular_inductance` | verified |
| `calc_circstrands_coords` | `wire_coordinates` | verified |
| `calc_inductance_trifoil` | `trifoil_inductance` | verified |
| `calc_circstrands_gmr` | `strand_gmr` | verified |
| `calc_tubular_gmr` | `tubular_gmr` | verified |
| `calc_equivalent_mu` | `equivalent_mu` | verified |
| `calc_shunt_capacitance` | `shunt_capacitance` | verified |
| `calc_shunt_conductance` | `shunt_conductance` | verified |
| `calc_equivalent_gmr` | `equivalent_gmr` | verified |
| `calc_gmd` | `gmd` | verified |
| `calc_solenoid_correction` | `solenoid_factor` | verified |
| `calc_equivalent_rho` | `equivalent_rho` | verified |
| `calc_equivalent_eps` | `equivalent_eps` | verified |
| `calc_equivalent_lossfact` | `loss_tangent` | verified |
| `calc_sigma_lossfact` | `equivalent_conductivity` | verified |
| `_calc_horz_sep!` | `horizontal_separation!` | verified |
| `_calc_gamma` | `gamma` | verified |
| `_calc_transformation_matrix_LM` | `levenberg_transform` | verified |
| `_calc_modal_quantities` | `modal_quantities` | verified |

## Source inventory

### Root and shared infrastructure

| Path | Disposition | Status | Evidence |
|---|---|---:|---|
| `Project.toml` | remove obsolete dependencies after use audit | verified | source and migrated tests |
| `src/LineCableModels.jl` | replace Commons/Utils wiring and exports | verified | source and migrated tests |
| `src/interfaces.jl` | add genuine root generics and domain tags | verified | source and migrated tests |
| `src/docstrings.jl` | move METHODLIST support; preserve CI path hiding | verified | source and migrated tests |
| `src/retired.jl` | consolidate FEM, sector, and legacy-JSON tombstones | verified | source and migrated tests |
| `src/commons/Commons.jl` | delete | verified | source and migrated tests |
| `src/commons/consts.jl` | delete | verified | source and migrated tests |
| `src/commons/docstringextension.jl` | move to root `docstrings.jl` | verified | source and migrated tests |
| `src/commons/retired.jl` | move to root `retired.jl` | verified | source and migrated tests |
| `src/utils/Utils.jl` | delete after moving owned concepts | verified | source and migrated tests |
| `src/utils/commondeps.jl` | delete | verified | source and migrated tests |
| `src/utils/logging.jl` | delete file/global logging | verified | source and migrated tests |
| `src/utils/macros.jl` | delete | verified | source and migrated tests |
| `src/utils/typecoercion.jl` | delete | verified | source and migrated tests |
| `src/validation/Validation.jl` | root `validate` plus declarative rule namespace | verified | source and migrated tests |
| `src/validation/rules.jl` | retain noncoercive primitive rules | verified | source and migrated tests |
| `src/validation/applyrules.jl` | nonmutating native-exception evaluation | verified | source and migrated tests |
| `src/unithandler/UnitHandler.jl` | retain; import root accessors | verified | source and migrated tests |
| `src/plotbuilder/PlotBuilder.jl` | preserve grammar; replace Utils imports | verified | source and migrated tests |
| `src/plotbuilder/types.jl` | specialize root validation | verified | source and migrated tests |
| `src/plotbuilder/grammar.jl` | retain behavior | verified | source and migrated tests |
| `src/plotbuilder/backendhandler/BackendHandler.jl` | retain | verified | source and migrated tests |
| `src/plotbuilder/uicomponents/UIComponents.jl` | retain visuals; replace shared accessors | verified | source and migrated tests |

### Materials and DataModel

| Path | Disposition | Status | Evidence |
|---|---|---:|---|
| `src/materials/Materials.jl` | natural Real types, promotion, validation; retain T0/alpha | verified | source and migrated tests |
| `src/materials/base.jl` | retain behavior; own explicit conversions | verified | source and migrated tests |
| `src/materials/dataframe.jl` | retain | verified | source and migrated tests |
| `src/materials/materialslibrary.jl` | dictionary semantics; reject duplicate `add!` | verified | source and migrated tests |
| `src/materials/typecoercion.jl` | delete | verified | source and migrated tests |
| `src/datamodel/DataModel.jl` | remove macro/coercion plumbing; retain materialized API | verified | source and migrated tests |
| `src/datamodel/types.jl` | retain hierarchy; remove catch-all interception | verified | source and migrated tests |
| `src/datamodel/radii.jl` | delete; radial intent belongs to Specs | verified | source and migrated tests |
| `src/datamodel/macros.jl` | delete | verified | source and migrated tests |
| `src/datamodel/validation.jl` | delete constructor-default glue | verified | source and migrated tests |
| `src/datamodel/retired.jl` | move tombstones to root | verified | source and migrated tests |
| `src/datamodel/strands_handler.jl` | own functional `maxfill` dispatch | verified | source and migrated tests |
| `src/datamodel/typecoercion.jl` | delete | verified | source and migrated tests |
| `src/datamodel/helpers.jl` | replace with `geometry.jl` | verified | source and migrated tests |
| `src/datamodel/geometry.jl` | add formations, overlap geometry, `outer_radius` | verified | source and migrated tests |
| `src/datamodel/io.jl` | retain display; RectStrands-only shape forwarding | verified | source and migrated tests |
| `src/datamodel/plotspecs.jl` | retain visuals; adapt owned fields/accessors | verified | source and migrated tests |
| `src/datamodel/baseparams/BaseParams.jl` | explicit exports; remove coercion machinery | verified | source and migrated tests |
| `src/datamodel/baseparams/geometry.jl` | add geometry formula owner | verified | source and migrated tests |
| `src/datamodel/baseparams/resistance.jl` | add resistance/temperature formula owner | verified | source and migrated tests |
| `src/datamodel/baseparams/inductance.jl` | add inductance/permeability owner | verified | source and migrated tests |
| `src/datamodel/baseparams/dielectrics.jl` | add dielectric formula owner | verified | source and migrated tests |
| `src/datamodel/circstrands.jl` | remove operating temperature; base-state formulas | verified | source and migrated tests |
| `src/datamodel/rectstrands.jl` | remove operating temperature; retain geometry | verified | source and migrated tests |
| `src/datamodel/strip.jl` | exact radii and base-state formulas | verified | source and migrated tests |
| `src/datamodel/tubular.jl` | exact radii and base-state formulas | verified | source and migrated tests |
| `src/datamodel/insulator.jl` | exact radii; correct conductance units | verified | source and migrated tests |
| `src/datamodel/semicon.jl` | exact radii; correct conductance units | verified | source and migrated tests |
| `src/datamodel/conductorgroup.jl` | strict insertion, common T0, fix undefined `Tg` | verified | source and migrated tests |
| `src/datamodel/conductorgroup/base.jl` | fold `eltype`, display, and conversion into `conductorgroup.jl`; delete split file | verified | source and migrated tests |
| `src/datamodel/insulatorgroup.jl` | explicit reference frequency; strict insertion | verified | source and migrated tests |
| `src/datamodel/insulatorgroup/base.jl` | fold `eltype`, display, and conversion into `insulatorgroup.jl`; delete split file | verified | source and migrated tests |
| `src/datamodel/nominaldata.jl` | naturally promote present optional values | verified | source and migrated tests |
| `src/datamodel/nominaldata/base.jl` | fold `eltype`, display, and conversion into `nominaldata.jl`; delete split file | verified | source and migrated tests |
| `src/datamodel/cablecomponent.jl` | preserve eager equivalent flattening at reference state | verified | source and migrated tests |
| `src/datamodel/cablecomponent/base.jl` | retain display/eltype; own conversion | verified | source and migrated tests |
| `src/datamodel/cabledesign.jl` | strict insertion, common T0, preserve transforms | verified | source and migrated tests |
| `src/datamodel/cabledesign/base.jl` | retain display/indexing; fix absent nominal data | verified | source and migrated tests |
| `src/datamodel/cabledesign/cableconstants.jl` | preserve base-temperature results | verified | source and migrated tests |
| `src/datamodel/cabledesign/dataframe.jl` | retain CableConstants delegation | verified | source and migrated tests |
| `src/datamodel/cableslibrary.jl` | retain storage and `vdeparse`; dictionary semantics | verified | source and migrated tests |
| `src/datamodel/cableslibrary/base.jl` | retain operations/concise display | verified | source and migrated tests |
| `src/datamodel/cableslibrary/dataframe.jl` | retain | verified | source and migrated tests |
| `src/datamodel/cableslibrary/vdeparse.jl` | retain mappings, parsing, and exports | verified | source and migrated tests |
| `src/datamodel/linecablesystem.jl` | strict insertion, validation, derived counts | verified | source and migrated tests |
| `src/datamodel/linecablesystem/base.jl` | fold display, `eltype`, counts, and conversion into `linecablesystem.jl`; delete split file | verified | source and migrated tests |
| `src/datamodel/linecablesystem/dataframe.jl` | retain with derived counts | verified | source and migrated tests |

### EarthProps and ParametricBuilder

| Path | Disposition | Status | Evidence |
|---|---|---:|---|
| `src/earthprops/EarthProps.jl` | static EarthLayer/EarthModel only | verified | source and migrated tests |
| `src/earthprops/base.jl` | retain display/eltype; own conversion | verified | source and migrated tests |
| `src/earthprops/dataframe.jl` | retain static summary | verified | source and migrated tests |
| `src/earthprops/fdprops.jl` | move frequency formulation into Engine | verified | source and migrated tests |
| `src/earthprops/typecoercion.jl` | delete | verified | source and migrated tests |
| `src/parametricbuilder/ParametricBuilder.jl` | retain DSL, GridSpace, and quarantined `@relax` | verified | source and migrated tests |
| `src/parametricbuilder/gridspace/grid.jl` | retain semantics; root nominal accessors | verified | source and migrated tests |
| `src/parametricbuilder/gridspace/gridspace.jl` | retain configuration/materialization grammar | verified | source and migrated tests |
| `src/parametricbuilder/gridspace/macros.jl` | retain `@gridspace`; quarantine `@relax`/`recast` | verified | source and migrated tests |
| `src/parametricbuilder/materialspec.jl` | explicit unexported MaterialSpec | verified | source and migrated tests |
| `src/parametricbuilder/cablebuilderspec.jl` | explicit PartSpec and dispatch materialization | verified | source and migrated tests |
| `src/parametricbuilder/positionspec.jl` | rename PositionBuilder to PositionSpec | verified | source and migrated tests |
| `src/parametricbuilder/systembuilderspec.jl` | EarthSpec/SystemSpec; sole operating temperature intent | verified | source and migrated tests |
| `src/parametricbuilder/wirepatterns/WirePatterns.jl` | typed patterns, shared maxfill, WireEstimate | verified | source and migrated tests |

### Engine

| Path | Disposition | Status | Evidence |
|---|---|---:|---|
| `src/engine/Engine.jl` | explicit imports/exports; expose trace/property types | verified | source and migrated tests |
| `src/engine/types.jl` | natural types; parametric formulations | verified | source and migrated tests |
| `src/engine/lineparamopts.jl` | typed verbosity/output basis; remove storage/logfile | verified | source and migrated tests |
| `src/engine/problemdefs.jl` | declarative validation; retain problem grammar | verified | source and migrated tests |
| `src/engine/lineparams.jl` | preserve containers/accessors/basis/domain | verified | source and migrated tests |
| `src/engine/base.jl` | remove workspace methods; add EMTTrace display | verified | source and migrated tests |
| `src/engine/dataframe.jl` | retain; root nominal/uncertainty accessors | verified | source and migrated tests |
| `src/engine/plotspecs.jl` | retain visuals; update imports | verified | source and migrated tests |
| `src/engine/helpers.jl` | delete after owner-local moves | verified | source and migrated tests |
| `src/engine/workspace.jl` | delete | verified | source and migrated tests |
| `src/engine/input.jl` | add immutable flattened EMTInput | verified | source and migrated tests |
| `src/engine/trace.jl` | add EMTTrace | verified | source and migrated tests |
| `src/engine/earthkernels.jl` | add shared Engine earth-return kernels | verified | source and migrated tests |
| `src/engine/earthproperties/EarthProperties.jl` | add Engine formulation namespace | verified | source and migrated tests |
| `src/engine/earthproperties/constant.jl` | add Engine-owned CPEarth | verified | source and migrated tests |
| `src/engine/solver.jl` | typed sequence, temperature, inspection, transposition fix | verified | source and migrated tests |
| `src/engine/reduction.jl` | retain algorithms; generic types and proven allocation fixes | verified | source and migrated tests |
| `src/engine/retired.jl` | move to root tombstones | verified | source and migrated tests |
| `src/engine/internalimpedance/InternalImpedance.jl` | retain API; replace shared imports | verified | source and migrated tests |
| `src/engine/internalimpedance/scaledbessel.jl` | preserve equations; generic/local rho | verified | source and migrated tests |
| `src/engine/insulationimpedance/InsulationImpedance.jl` | retain API; replace imports | verified | source and migrated tests |
| `src/engine/insulationimpedance/lossless.jl` | preserve target-typed formula | verified | source and migrated tests |
| `src/engine/earthimpedance/EarthImpedance.jl` | Engine kernels; scoped QuadGK boundary | verified | source and migrated tests |
| `src/engine/earthimpedance/base.jl` | retain | verified | source and migrated tests |
| `src/engine/earthimpedance/homogeneous.jl` | preserve methods; generic quadrature | verified | source and migrated tests |
| `src/engine/insulationadmittance/InsulationAdmittance.jl` | owner-local constants/conductivity | verified | source and migrated tests |
| `src/engine/insulationadmittance/base.jl` | retain | verified | source and migrated tests |
| `src/engine/insulationadmittance/lossless.jl` | preserve | verified | source and migrated tests |
| `src/engine/insulationadmittance/parallelrc.jl` | physical-layer broadband calculation | verified | source and migrated tests |
| `src/engine/earthadmittance/EarthAdmittance.jl` | Engine kernels; scoped QuadGK boundary | verified | source and migrated tests |
| `src/engine/earthadmittance/base.jl` | retain | verified | source and migrated tests |
| `src/engine/earthadmittance/homogeneous.jl` | preserve methods; generic quadrature | verified | source and migrated tests |
| `src/engine/ehem/EHEM.jl` | retain dispatch; consume evaluated earth properties | verified | source and migrated tests |
| `src/engine/ehem/enforcelayer.jl` | preserve behavior; native validation | verified | source and migrated tests |
| `src/engine/transforms/Transforms.jl` | own reciprocity/transposition diagnostics | verified | source and migrated tests |
| `src/engine/transforms/fortescue.jl` | preserve results | verified | source and migrated tests |
| `src/engine/transforms/eiglevenberg.jl` | preserve algorithm; scoped NLsolve verbosity | verified | source and migrated tests |

### Downstream and extensions

| Path | Disposition | Status | Evidence |
|---|---|---:|---|
| `src/computation/Computation.jl` | retain policies/results; ordinary compute path | verified | source and migrated tests |
| `src/computation/dataframe.jl` | retain | verified | source and migrated tests |
| `src/computation/plotspecs.jl` | retain | verified | source and migrated tests |
| `src/importexport/ImportExport.jl` | own external format adapters and schemas | verified | source and migrated tests |
| `src/importexport/serialize.jl` | explicit versioned type tags | verified | source and migrated tests |
| `src/importexport/deserialize.jl` | registry parsing; retired legacy detection | verified | source and migrated tests |
| `src/importexport/materialslibrary.jl` | supported-schema atomic save/load | verified | source and migrated tests |
| `src/importexport/cableslibrary.jl` | supported JSON/JLS; atomic loading | verified | source and migrated tests |
| `src/importexport/pscad.jl` | public entry point for static emission | verified | source and migrated tests |
| `src/importexport/pscad/schema.jl` | add minimal structure/master bindings | verified | source and migrated tests |
| `src/importexport/pscad/project.jl` | add EzXML component/project emitter | verified | source and migrated tests |
| `src/importexport/atp.jl` | retain output; static Earth/renames | verified | source and migrated tests |
| `src/importexport/tralin.jl` | retain output and formatting | verified | source and migrated tests |
| `src/importexport/xlsx.jl` | retain output; owner-local diagonal check | verified | source and migrated tests |
| `ext/LineCableModelsMeasurementsExt.jl` | root accessors, persistence, kernels, MC | verified | source and migrated tests |
| `ext/LineCableModelsDistributionsExt.jl` | retain | verified | source and migrated tests |
| `ext/LineCableModelsMakieExt.jl` | retain | verified | source and migrated tests |
| `ext/LineCableModelsCairoMakieExt.jl` | retain | verified | source and migrated tests |
| `ext/LineCableModelsGLMakieExt.jl` | retain | verified | source and migrated tests |
| `ext/LineCableModelsWGLMakieExt.jl` | retain | verified | source and migrated tests |

### Empty directories

| Path | Disposition | Status |
|---|---|---:|
| `ext/fem/getdp_frontend` | remove | verified |
| `integration/fem` | remove | verified |
| `src/engine/fem` | remove | verified |
| `src/plotbuilder/plotuicomponents` | remove | verified |
| `src/uncertainbessels` | remove | verified |
| `src/uq/plotspecs` | remove | verified |

## Numerical parity and intentional deltas

| Case | Baseline artifact | Expected result | Status | Evidence |
|---|---|---|---:|---|
| eager cable-part/group formulas | immutable JSON fixture at baseline SHA | unchanged at material T0 | verified | `datamodel_engine_refactor.jl` component comparisons |
| CableConstants | immutable JSON fixture at baseline SHA | unchanged | verified | R/L/C scaled comparisons |
| maintained EMT formulations | immutable JSON fixture plus formulation suites | unchanged | verified | default, reduced, unreduced, modal, earth, and insulation tests |
| reductions and modal transforms | immutable JSON fixture plus transform suites | unchanged | verified | primitive and result tensor comparisons |
| Earth impedance/admittance paths | immutable JSON fixture plus formulation suites | unchanged | verified | homogeneous earth and CPEarth tests |
| supported external formats | baseline examples and donor structures | semantically unchanged | verified | round trips and semantic ATP/PSCAD/TRALIN/XLSX tests |
| ideal transposition | crossed baseline Z/Y tensors | intentionally correct Z/Y condition in all paths | verified | approved delta isolated in `datamodel_engine_refactor.jl` |
| non-base operating temperature | baseline reference-state fixture | correction moves from eager objects into compute | verified | cold/hot/repeated/disabled end-to-end cases |

Additional numerical deltas require user approval before implementation.

## Performance ledger

| Case | Baseline | Limit | Result | Status |
|---|---:|---:|---:|---:|
| representative Float64 EMT median runtime | 1.654097 ms | ≤ 1.819507 ms | 1.558481 ms | verified |
| representative Float64 EMT allocations | 180,000 bytes | ≤ 189,000 bytes | 169,216 bytes | verified |
| `inspect=false` trace tensors | current workspace stores optional tensors | allocate none | zero-byte `_trace_buffers(..., false)` | verified |

## Test migration ledger

No legacy test disappears without an observable replacement.

| Existing area | Replacement evidence | Status |
|---|---|---:|
| Commons domain/tombstone/docstring tests | root interfaces/retired/docstrings tests | verified |
| Utils logging tests | scoped Engine logging/global logger tests | verified |
| Utils macros/coercion tests | promotion, convert, `@relax`, architecture guards | verified |
| Validation rule tests | native exceptions and declarative owner rules | verified |
| cable-part constructor suites | exact radii, promotion, base-state temperature | verified |
| group/design/system suites | strict mutation, atomic failure, shared T0, derived counts | verified |
| EarthProps suite | static Earth plus Engine frequency evaluation | verified |
| Engine integration | parity fixtures, trace, temperature, transposition | verified |
| ParametricBuilder suites | explicit Specs, whole-intent promotion, GridSpace parity | verified |
| persistence/export suites | versioned atomic JSON/JLS and semantic exporters | verified |
| plotting suites/goldens | unchanged grammar and visuals | verified |

## Architecture guards

The completed suite must reject Commons, Utils, `helpers.jl`, internal shim modules,
scalar union aliases, coercion machinery, typed-constant wrappers, removed constructor
macros, `calc_` outside tombstones, executable `@assert`, dynamic JSON evaluation,
operating-temperature fields through CableDesign, global logger mutation, and use of
`@relax`/`recast` outside its quarantined implementation and dedicated tests.

The architecture guard in `test/unit/core/architecture.jl` and a final repository
scan verify these exclusions. `@relax` and its private reconstruction machinery occur
only in their implementation, export declaration, and dedicated GridSpace/core tests.

## Verification summary

- Default migrated package suite: 2,860/2,860 passing. The lower raw count than the
  2,962-test baseline reflects removal of obsolete Commons/Utils/API assertions and
  consolidation into the replacement rows above; no inventoried behavior was dropped.
- Core-only import and tests: 32/32 passing; core import loads no plotting backend.
- Aqua: 11/11 passing.
- Cairo visual gate: 446/446 passing. Only `line_zy_polar.png` was refreshed; its
  three-pixel legend clearance was introduced at the immutable baseline before this
  refactor, so the stale reference did not represent a new visual change here.
- Native GL and WGL smoke galleries: 18/18 panels built for each backend. The gallery
  fixture now uses the static `EarthModel` constructor.
- Documentation: doctest 1/1 and strict Documenter build passing without `warnonly`.
- SciML formatter: `src`, `ext`, `test`, `integration`, `examples`, and `docs` pass
  `overwrite=false`.
- Production coverage: 6,203/6,494 lines (95.52%), above the 95% gate.
- Clean installation: succeeded from an empty depot; `using LineCableModels` left all
  Makie backend extensions unloaded.
- Performance: 11 warm runs followed by 11 measured runs using the baseline case and
  protocol. Runtime improved by 5.8%; allocations improved by 6.0%.
- Dependency and architecture scans find none of the forbidden internal modules,
  coercion aliases, constructor macros, `calc_` implementations, executable source
  assertions, dynamic deserialization, or global logger mutation.

## Delivery sequence

1. `docs(refactor): inventory data model and engine changes`
2. `refactor(core): replace common and utility modules`
3. `refactor(datamodel): use base-state eager constructors`
4. `refactor(engine): apply system temperature in compute`
5. `refactor(io): harden persistence and pscad export`
6. `test(refactor): verify numerical and type parity`
7. `docs(refactor): document the materialized model pipeline`

These commits remain local until the user explicitly authorizes a remote action.
