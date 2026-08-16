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
- Numerical fixtures: pending capture before numerical source changes
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
| `set_verbosity!` | removed; use `ComputeOptions(verbosity=...)` | pending | — |
| `validate!` and scattered assertions | root `validate(x)` plus declarative rules | pending | — |
| `MaxFill()` | removed; use `maxfill(Type, ...)` | pending | — |
| `WireArray` | removed | pending | — |
| macro-generated builder parameter types | documented, unexported Specs | pending | — |
| `@gridspace` | retained | pending | — |
| `@relax` | retained only as unused future provision | pending | — |
| operating-temperature cable-part constructors | removed; set temperature on system problem | pending | — |
| `EMTWorkspace` and `_compute_with_workspace` | removed; use `compute!(...; inspect=true)` | pending | — |
| primitive-matrix storage option | removed; use `EMTTrace` | pending | — |
| file logging | removed; scoped console logging only | pending | — |
| unversioned/legacy JSON loading | deprecation notice and throwing tombstone | pending | — |
| maintained JLS/JSON/ATP/PSCAD/TRALIN/XLSX | retained | pending | — |
| `equivalent`, `nonsensify`, `vdeparse` | retained | pending | — |
| PlotBuilder grammar and visuals | retained | pending | — |

## Mathematical rename ledger

All replacements are breaking prerelease renames without aliases. Existing formulas,
mathematical notes, units, and citations remain intact.

| Current | Replacement | Status |
|---|---|---:|
| `calc_equivalent_alpha` | `equivalent_alpha` | pending |
| `calc_parallel_equivalent` | `parallel` | pending |
| `calc_helical_params` | `helix` | pending |
| `calc_strip_resistance` | `strip_resistance` | pending |
| `calc_temperature_correction` | `temperature_factor` | pending |
| `calc_tubular_resistance` | `tubular_resistance` | pending |
| `calc_tubular_inductance` | `tubular_inductance` | pending |
| `calc_circstrands_coords` | `wire_coordinates` | pending |
| `calc_inductance_trifoil` | `trifoil_inductance` | pending |
| `calc_circstrands_gmr` | `strand_gmr` | pending |
| `calc_tubular_gmr` | `tubular_gmr` | pending |
| `calc_equivalent_mu` | `equivalent_mu` | pending |
| `calc_shunt_capacitance` | `shunt_capacitance` | pending |
| `calc_shunt_conductance` | `shunt_conductance` | pending |
| `calc_equivalent_gmr` | `equivalent_gmr` | pending |
| `calc_gmd` | `gmd` | pending |
| `calc_solenoid_correction` | `solenoid_factor` | pending |
| `calc_equivalent_rho` | `equivalent_rho` | pending |
| `calc_equivalent_eps` | `equivalent_eps` | pending |
| `calc_equivalent_lossfact` | `loss_tangent` | pending |
| `calc_sigma_lossfact` | `equivalent_conductivity` | pending |
| `_calc_horz_sep!` | `horizontal_separation!` | pending |
| `_calc_gamma` | `gamma` | pending |
| `_calc_transformation_matrix_LM` | `levenberg_transform` | pending |
| `_calc_modal_quantities` | `modal_quantities` | pending |

## Source inventory

### Root and shared infrastructure

| Path | Disposition | Status | Evidence |
|---|---|---:|---|
| `Project.toml` | remove obsolete dependencies after use audit | pending | — |
| `src/LineCableModels.jl` | replace Commons/Utils wiring and exports | pending | — |
| `src/interfaces.jl` | add genuine root generics and domain tags | pending | — |
| `src/docstrings.jl` | move METHODLIST support; preserve CI path hiding | pending | — |
| `src/retired.jl` | consolidate FEM, sector, and legacy-JSON tombstones | pending | — |
| `src/commons/Commons.jl` | delete | pending | — |
| `src/commons/consts.jl` | delete | pending | — |
| `src/commons/docstringextension.jl` | move to root `docstrings.jl` | pending | — |
| `src/commons/retired.jl` | move to root `retired.jl` | pending | — |
| `src/utils/Utils.jl` | delete after moving owned concepts | pending | — |
| `src/utils/commondeps.jl` | delete | pending | — |
| `src/utils/logging.jl` | delete file/global logging | pending | — |
| `src/utils/macros.jl` | delete | pending | — |
| `src/utils/typecoercion.jl` | delete | pending | — |
| `src/validation/Validation.jl` | root `validate` plus declarative rule namespace | pending | — |
| `src/validation/rules.jl` | retain noncoercive primitive rules | pending | — |
| `src/validation/applyrules.jl` | nonmutating native-exception evaluation | pending | — |
| `src/unithandler/UnitHandler.jl` | retain; import root accessors | pending | — |
| `src/plotbuilder/PlotBuilder.jl` | preserve grammar; replace Utils imports | pending | — |
| `src/plotbuilder/types.jl` | specialize root validation | pending | — |
| `src/plotbuilder/grammar.jl` | retain behavior | pending | — |
| `src/plotbuilder/backendhandler/BackendHandler.jl` | retain | pending | — |
| `src/plotbuilder/uicomponents/UIComponents.jl` | retain visuals; replace shared accessors | pending | — |

### Materials and DataModel

| Path | Disposition | Status | Evidence |
|---|---|---:|---|
| `src/materials/Materials.jl` | natural Real types, promotion, validation; retain T0/alpha | pending | — |
| `src/materials/base.jl` | retain behavior; own explicit conversions | pending | — |
| `src/materials/dataframe.jl` | retain | pending | — |
| `src/materials/materialslibrary.jl` | dictionary semantics; reject duplicate `add!` | pending | — |
| `src/materials/typecoercion.jl` | delete | pending | — |
| `src/datamodel/DataModel.jl` | remove macro/coercion plumbing; retain materialized API | pending | — |
| `src/datamodel/types.jl` | retain hierarchy; remove catch-all interception | pending | — |
| `src/datamodel/radii.jl` | delete; radial intent belongs to Specs | pending | — |
| `src/datamodel/macros.jl` | delete | pending | — |
| `src/datamodel/validation.jl` | delete constructor-default glue | pending | — |
| `src/datamodel/retired.jl` | move tombstones to root | pending | — |
| `src/datamodel/strands_handler.jl` | own functional `maxfill` dispatch | pending | — |
| `src/datamodel/typecoercion.jl` | delete | pending | — |
| `src/datamodel/helpers.jl` | replace with `geometry.jl` | pending | — |
| `src/datamodel/geometry.jl` | add formations, overlap geometry, `outer_radius` | pending | — |
| `src/datamodel/io.jl` | retain display; RectStrands-only shape forwarding | pending | — |
| `src/datamodel/plotspecs.jl` | retain visuals; adapt owned fields/accessors | pending | — |
| `src/datamodel/baseparams/BaseParams.jl` | explicit exports; remove coercion machinery | pending | — |
| `src/datamodel/baseparams/geometry.jl` | add geometry formula owner | pending | — |
| `src/datamodel/baseparams/resistance.jl` | add resistance/temperature formula owner | pending | — |
| `src/datamodel/baseparams/inductance.jl` | add inductance/permeability owner | pending | — |
| `src/datamodel/baseparams/dielectrics.jl` | add dielectric formula owner | pending | — |
| `src/datamodel/circstrands.jl` | remove operating temperature; base-state formulas | pending | — |
| `src/datamodel/rectstrands.jl` | remove operating temperature; retain geometry | pending | — |
| `src/datamodel/strip.jl` | exact radii and base-state formulas | pending | — |
| `src/datamodel/tubular.jl` | exact radii and base-state formulas | pending | — |
| `src/datamodel/insulator.jl` | exact radii; correct conductance units | pending | — |
| `src/datamodel/semicon.jl` | exact radii; correct conductance units | pending | — |
| `src/datamodel/conductorgroup.jl` | strict insertion, common T0, fix undefined `Tg` | pending | — |
| `src/datamodel/conductorgroup/base.jl` | retain `eltype`; own conversion | pending | — |
| `src/datamodel/insulatorgroup.jl` | explicit reference frequency; strict insertion | pending | — |
| `src/datamodel/insulatorgroup/base.jl` | retain `eltype`; own conversion | pending | — |
| `src/datamodel/nominaldata.jl` | naturally promote present optional values | pending | — |
| `src/datamodel/nominaldata/base.jl` | retain `eltype`; own conversion | pending | — |
| `src/datamodel/cablecomponent.jl` | preserve eager equivalent flattening at reference state | pending | — |
| `src/datamodel/cablecomponent/base.jl` | retain display/eltype; own conversion | pending | — |
| `src/datamodel/cabledesign.jl` | strict insertion, common T0, preserve transforms | pending | — |
| `src/datamodel/cabledesign/base.jl` | retain display/indexing; fix absent nominal data | pending | — |
| `src/datamodel/cabledesign/cableconstants.jl` | preserve base-temperature results | pending | — |
| `src/datamodel/cabledesign/dataframe.jl` | retain CableConstants delegation | pending | — |
| `src/datamodel/cableslibrary.jl` | retain storage and `vdeparse`; dictionary semantics | pending | — |
| `src/datamodel/cableslibrary/base.jl` | retain operations/concise display | pending | — |
| `src/datamodel/cableslibrary/dataframe.jl` | retain | pending | — |
| `src/datamodel/cableslibrary/vdeparse.jl` | retain mappings, parsing, and exports | pending | — |
| `src/datamodel/linecablesystem.jl` | strict insertion, validation, derived counts | pending | — |
| `src/datamodel/linecablesystem/base.jl` | retain display/eltype; derived counts/conversion | pending | — |
| `src/datamodel/linecablesystem/dataframe.jl` | retain with derived counts | pending | — |

### EarthProps and ParametricBuilder

| Path | Disposition | Status | Evidence |
|---|---|---:|---|
| `src/earthprops/EarthProps.jl` | static EarthLayer/EarthModel only | pending | — |
| `src/earthprops/base.jl` | retain display/eltype; own conversion | pending | — |
| `src/earthprops/dataframe.jl` | retain static summary | pending | — |
| `src/earthprops/fdprops.jl` | move frequency formulation into Engine | pending | — |
| `src/earthprops/typecoercion.jl` | delete | pending | — |
| `src/parametricbuilder/ParametricBuilder.jl` | retain DSL, GridSpace, and quarantined `@relax` | pending | — |
| `src/parametricbuilder/gridspace/grid.jl` | retain semantics; root nominal accessors | pending | — |
| `src/parametricbuilder/gridspace/gridspace.jl` | retain configuration/materialization grammar | pending | — |
| `src/parametricbuilder/gridspace/macros.jl` | retain `@gridspace`; quarantine `@relax`/`recast` | pending | — |
| `src/parametricbuilder/materialspec.jl` | explicit unexported MaterialSpec | pending | — |
| `src/parametricbuilder/cablebuilderspec.jl` | explicit PartSpec and dispatch materialization | pending | — |
| `src/parametricbuilder/positionspec.jl` | rename PositionBuilder to PositionSpec | pending | — |
| `src/parametricbuilder/systembuilderspec.jl` | EarthSpec/SystemSpec; sole operating temperature intent | pending | — |
| `src/parametricbuilder/wirepatterns/WirePatterns.jl` | typed patterns, shared maxfill, WireEstimate | pending | — |

### Engine

| Path | Disposition | Status | Evidence |
|---|---|---:|---|
| `src/engine/Engine.jl` | explicit imports/exports; expose trace/property types | pending | — |
| `src/engine/types.jl` | natural types; parametric formulations | pending | — |
| `src/engine/lineparamopts.jl` | typed verbosity/output basis; remove storage/logfile | pending | — |
| `src/engine/problemdefs.jl` | declarative validation; retain problem grammar | pending | — |
| `src/engine/lineparams.jl` | preserve containers/accessors/basis/domain | pending | — |
| `src/engine/base.jl` | remove workspace methods; add EMTTrace display | pending | — |
| `src/engine/dataframe.jl` | retain; root nominal/uncertainty accessors | pending | — |
| `src/engine/plotspecs.jl` | retain visuals; update imports | pending | — |
| `src/engine/helpers.jl` | delete after owner-local moves | pending | — |
| `src/engine/workspace.jl` | delete | pending | — |
| `src/engine/input.jl` | add immutable flattened EMTInput | pending | — |
| `src/engine/trace.jl` | add EMTTrace | pending | — |
| `src/engine/earthkernels.jl` | add shared Engine earth-return kernels | pending | — |
| `src/engine/earthproperties/EarthProperties.jl` | add Engine formulation namespace | pending | — |
| `src/engine/earthproperties/constant.jl` | add Engine-owned CPEarth | pending | — |
| `src/engine/solver.jl` | typed sequence, temperature, inspection, transposition fix | pending | — |
| `src/engine/reduction.jl` | retain algorithms; generic types and proven allocation fixes | pending | — |
| `src/engine/retired.jl` | move to root tombstones | pending | — |
| `src/engine/internalimpedance/InternalImpedance.jl` | retain API; replace shared imports | pending | — |
| `src/engine/internalimpedance/scaledbessel.jl` | preserve equations; generic/local rho | pending | — |
| `src/engine/insulationimpedance/InsulationImpedance.jl` | retain API; replace imports | pending | — |
| `src/engine/insulationimpedance/lossless.jl` | preserve target-typed formula | pending | — |
| `src/engine/earthimpedance/EarthImpedance.jl` | Engine kernels; scoped QuadGK boundary | pending | — |
| `src/engine/earthimpedance/base.jl` | retain | pending | — |
| `src/engine/earthimpedance/homogeneous.jl` | preserve methods; generic quadrature | pending | — |
| `src/engine/insulationadmittance/InsulationAdmittance.jl` | owner-local constants/conductivity | pending | — |
| `src/engine/insulationadmittance/base.jl` | retain | pending | — |
| `src/engine/insulationadmittance/lossless.jl` | preserve | pending | — |
| `src/engine/insulationadmittance/parallelrc.jl` | physical-layer broadband calculation | pending | — |
| `src/engine/earthadmittance/EarthAdmittance.jl` | Engine kernels; scoped QuadGK boundary | pending | — |
| `src/engine/earthadmittance/base.jl` | retain | pending | — |
| `src/engine/earthadmittance/homogeneous.jl` | preserve methods; generic quadrature | pending | — |
| `src/engine/ehem/EHEM.jl` | retain dispatch; consume evaluated earth properties | pending | — |
| `src/engine/ehem/enforcelayer.jl` | preserve behavior; native validation | pending | — |
| `src/engine/transforms/Transforms.jl` | own reciprocity/transposition diagnostics | pending | — |
| `src/engine/transforms/fortescue.jl` | preserve results | pending | — |
| `src/engine/transforms/eiglevenberg.jl` | preserve algorithm; scoped NLsolve verbosity | pending | — |

### Downstream and extensions

| Path | Disposition | Status | Evidence |
|---|---|---:|---|
| `src/computation/Computation.jl` | retain policies/results; ordinary compute path | pending | — |
| `src/computation/dataframe.jl` | retain | pending | — |
| `src/computation/plotspecs.jl` | retain | pending | — |
| `src/importexport/ImportExport.jl` | own external format adapters and schemas | pending | — |
| `src/importexport/serialize.jl` | explicit versioned type tags | pending | — |
| `src/importexport/deserialize.jl` | registry parsing; retired legacy detection | pending | — |
| `src/importexport/materialslibrary.jl` | supported-schema atomic save/load | pending | — |
| `src/importexport/cableslibrary.jl` | supported JSON/JLS; atomic loading | pending | — |
| `src/importexport/pscad.jl` | public entry point for static emission | pending | — |
| `src/importexport/pscad/schema.jl` | add minimal structure/master bindings | pending | — |
| `src/importexport/pscad/project.jl` | add EzXML component/project emitter | pending | — |
| `src/importexport/atp.jl` | retain output; static Earth/renames | pending | — |
| `src/importexport/tralin.jl` | retain output and formatting | pending | — |
| `src/importexport/xlsx.jl` | retain output; owner-local diagonal check | pending | — |
| `ext/LineCableModelsMeasurementsExt.jl` | root accessors, persistence, kernels, MC | pending | — |
| `ext/LineCableModelsDistributionsExt.jl` | retain | pending | — |
| `ext/LineCableModelsMakieExt.jl` | retain | pending | — |
| `ext/LineCableModelsCairoMakieExt.jl` | retain | pending | — |
| `ext/LineCableModelsGLMakieExt.jl` | retain | pending | — |
| `ext/LineCableModelsWGLMakieExt.jl` | retain | pending | — |

### Empty directories

| Path | Disposition | Status |
|---|---|---:|
| `ext/fem/getdp_frontend` | remove | pending |
| `integration/fem` | remove | pending |
| `src/engine/fem` | remove | pending |
| `src/plotbuilder/plotuicomponents` | remove | pending |
| `src/uncertainbessels` | remove | pending |
| `src/uq/plotspecs` | remove | pending |

## Numerical parity and intentional deltas

| Case | Baseline artifact | Expected result | Status | Evidence |
|---|---|---|---:|---|
| eager cable-part/group formulas | pending | unchanged at material T0 | pending | — |
| CableConstants | pending | unchanged | pending | — |
| maintained EMT formulations | pending | unchanged | pending | — |
| reductions and modal transforms | pending | unchanged | pending | — |
| Earth impedance/admittance paths | pending | unchanged | pending | — |
| supported external formats | pending | semantically unchanged | pending | — |
| ideal transposition | pending | intentionally correct Z/Y condition in all paths | approved | accepted plan |
| non-base operating temperature | pending | correction moves from eager objects into compute | approved | accepted plan |

Additional numerical deltas require user approval before implementation.

## Performance ledger

| Case | Baseline | Limit | Result | Status |
|---|---:|---:|---:|---:|
| representative Float64 EMT median runtime | 1.654097 ms | ≤ 1.819507 ms | — | pending |
| representative Float64 EMT allocations | 180,000 bytes | ≤ 189,000 bytes | — | pending |
| `inspect=false` trace tensors | current workspace stores optional tensors | allocate none | — | pending |

## Test migration ledger

No legacy test disappears without an observable replacement.

| Existing area | Replacement evidence | Status |
|---|---|---:|
| Commons domain/tombstone/docstring tests | root interfaces/retired/docstrings tests | pending |
| Utils logging tests | scoped Engine logging/global logger tests | pending |
| Utils macros/coercion tests | promotion, convert, `@relax`, architecture guards | pending |
| Validation rule tests | native exceptions and declarative owner rules | pending |
| cable-part constructor suites | exact radii, promotion, base-state temperature | pending |
| group/design/system suites | strict mutation, atomic failure, shared T0, derived counts | pending |
| EarthProps suite | static Earth plus Engine frequency evaluation | pending |
| Engine integration | parity fixtures, trace, temperature, transposition | pending |
| ParametricBuilder suites | explicit Specs, whole-intent promotion, GridSpace parity | pending |
| persistence/export suites | versioned atomic JSON/JLS and semantic exporters | pending |
| plotting suites/goldens | unchanged grammar and visuals | pending |

## Architecture guards

The completed suite must reject Commons, Utils, `helpers.jl`, internal shim modules,
scalar union aliases, coercion machinery, typed-constant wrappers, removed constructor
macros, `calc_` outside tombstones, executable `@assert`, dynamic JSON evaluation,
operating-temperature fields through CableDesign, global logger mutation, and use of
`@relax`/`recast` outside its quarantined implementation and dedicated tests.

## Delivery sequence

1. `docs(refactor): inventory data model and engine changes`
2. `refactor(core): replace common and utility modules`
3. `refactor(datamodel): use base-state eager constructors`
4. `refactor(engine): apply system temperature in compute`
5. `refactor(io): harden persistence and pscad export`
6. `test(refactor): verify numerical and type parity`
7. `docs(refactor): document the materialized model pipeline`

These commits remain local until the user explicitly authorizes a remote action.
