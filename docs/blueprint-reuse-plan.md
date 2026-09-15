# Explicit blueprint reuse — deferred implementation plan

Status: **DEFERRED — NOT AUTHORIZED FOR IMPLEMENTATION**.

Recorded: 2026-09-15.

Architectural debt cleanup takes priority over this optimization. This document
records a future implementation plan, not an instruction to start coding. Do
not add scaffolding, options, provenance fields, cache objects, dependencies or
benchmark changes now. Completing prerequisite cleanup does not automatically
activate this plan; implementation requires an explicit subsequent request.

The existing [shunt refactor record](shunt-model-refactor.md) describes delivered
behavior. Until this plan is implemented, ordinary compute calls continue to
construct fresh blueprints. All new API examples below are proposals.

## 1. Objective and scope

Allow a caller to recover completed coaxial cable blueprints from a traced
calculation and explicitly supply them to later compatible calculations.
Successful injection bypasses blueprint construction, including local boundary
solves, without changing the selected physical model or numerical results.

The reusable object remains `CableBlueprint`. It owns completed conductor
equivalents, dielectric layers, local assembly/terminal indexing, selected
boundary C/P coefficients and construction outcomes. No parallel preparation
object or public preparation call is introduced.

Cover both coaxial consumers, `LineParameters` and `CableConstants`, including
their scalar and formulation-collection entry points. Preserve the ordinary
no-injection path and existing within-call sharing.

Out of scope:

- Automatic global caches, eviction policies, process-wide lookup tables,
  partial cache completion or searching a pool for a matching design.
- Disk caches, portable checkpoints, cross-session reuse or compatibility
  across package/source-code revisions.
- Reuse of workspaces, earth responses, frequency arrays, reductions or dense
  boundary factorizations through this interface.
- FEM mesh/solver caching or accepting coaxial blueprints in FEM options.
- Changes to boundary equations, resolution, tolerances, audit requirements,
  fallback policy, material laws, or the default annular approximation.
- Cold-compilation optimization or dependency/environment changes.
- A new blueprint-only tracing mode. Existing trace memory costs will be
  documented and measured before considering another retention control.

## 2. Prerequisites and activation gate

The architecture cleanup is a separate workstream. This plan does not claim
to enumerate or authorize all of that debt. Before implementation:

- [ ] The user considers the relevant architectural cleanup complete and
  explicitly authorizes this deferred optimization.
- [ ] Reassess the then-current code; do not assume today's files, field names,
  numerical ownership or formula identifiers survived cleanup unchanged.
- [ ] Confirm that blueprint construction remains frequency-independent and
  that workspaces only consume, not modify, its completed numerical data.
- [ ] Stabilize ownership of resolved physical state, local formulation
  identity, uncertainty-aware comparisons and result-versus-runtime metadata.
- [ ] Agree how the existing owners represent an immutable construction-input
  snapshot. Reuse that representation if available; do not build a separate
  cache/provenance framework to compensate for unresolved ownership.
- [ ] Establish passing baseline checks for the touched owners. Recheck the
  previously reported Aqua ambiguities; do not assume they remain present or
  use this optimization to hide unrelated failures.
- [ ] Re-measure the remaining boundary construction cost. If cleanup already
  removes the practical need, reassess the benefit before adding API surface.

The behavior specified below is the implementation target. Concrete provenance
field layout and internal method placement must follow the cleaned-up owners;
that prerequisite is deliberately not settled by inventing new types here.

## 3. Dependency contract

The current blueprint is static for a completed cable design, its scalar
representation and its selected local construction formulas/settings. It is
not identified by a system name or cable name alone, nor is it independent of
every formulation choice.

Current construction entry points are in
[blueprint.jl](../src/engine/blueprint.jl) and
[blueprint_shunt.jl](../src/engine/blueprint_shunt.jl).
[lineinput](../src/engine/input.jl) subsequently adds system placement,
connections and frequency coordinates. Reconfirm these boundaries at activation.

| Input changed between calls | Reuse policy |
|---|---|
| Frequency grid or longitudinal propagation constants | Allowed if scalar representation remains compatible; rebuild frequency-dependent inputs. |
| Operating temperature or conductor temperature-correction selection | Allowed under the current blueprint contract; reevaluate runtime corrections. |
| Earth properties, equivalent-earth choices or earth-return formulas | Allowed; rebuild earth-dependent work. |
| External placement, spacing or rotation of unchanged cable designs | Allowed; rebuild system geometry, clearance handling and terminal placement. |
| Connections, grounding, reductions, transposition, line length or output basis | Allowed if the retained design terminals are unchanged; rebuild downstream work. |
| Series-impedance formulas evaluated only at runtime | Allowed after normal formulation validation. |
| Cable-internal geometry, material state, assembly topology or retained terminals | Reject stale blueprints. |
| Shunt model | Require an exact compatible construction selection. |
| Boundary resolution, integration controls, audit or fallback policy | Require matching normalized construction settings. |
| Insulation/semiconductor laws | Runtime-only for annular blueprints; boundary blueprints require matching admitted lossless selections. |
| Scalar type, precision context affecting construction, or uncertainty-source dependencies | Require an exact supported match; do not silently convert a cached numerical payload. |

External placement is distinct from changing positions of assemblies inside a
cable design. Only the former lies outside the local blueprint.

The current boundary model admits unmodified lossless laws, including their
default aliases, so completed C/P coefficients are independent of operating
frequency and temperature. A future dispersive or temperature-dependent
boundary model must revisit this contract; it cannot silently inherit permission
to cache completed C/P for all operating conditions.

## 4. Proposed public behavior

### 4.1 Retain through the existing result

Keep `compute`'s ordinary result type and return shape. With `trace=true`, retain
the completed blueprint collection under `details(result).trace.blueprints`:

```julia
# Proposed API; not currently implemented.
formulation = Formulation(shunt_model=:boundary)
result = @time compute(problem, formulation; options=(trace=true,))
blueprints = details(result).trace.blueprints

next_result = @time compute(next_problem, formulation;
    options=(blueprints=blueprints,))
```

Use the plural `blueprints` consistently: a line system has one entry per
selected cable design in `problem.system.designs` order. For `CableConstants`,
the same contract uses a one-element collection:

```julia
# Proposed API; not currently implemented.
local_formulation = CableConstantsFormulation(shunt_model=:boundary)
local_result = @time compute(local_problem, local_formulation;
    options=(trace=true,))
local_blueprints = details(local_result).trace.blueprints

next_local_result = @time compute(next_local_problem, local_formulation;
    options=(blueprints=local_blueprints,))
```

`CableConstants` currently accepts no compute options; add retention/injection
to its existing computation-options owner rather than inheriting unrelated line
or FEM execution controls. Existing line trace matrices remain available.

### 4.2 Injection rules

- `blueprints=nothing`, including omission, means ordinary construction.
- Accept one ordered vector of compatible `CableBlueprint{T}` values, not a
  heterogeneous cache object, dictionary of IDs, partial vector or callback.
- Validate collection length, order, scalar representation and construction
  compatibility before any selected calculation begins. The initial contract
  does not automatically reorder entries or match renamed designs.
- When every supplied blueprint is compatible, use those completed values
  directly. Do not call `flatten`, radial lowering, boundary geometry extraction,
  material descriptor evaluation for the boundary solve, or its numerical solver.
- When any supplied blueprint is incompatible, throw an informative
  `ArgumentError` or structural `DimensionMismatch`. Do not silently rebuild,
  drop the option, downgrade the model, or fill in missing entries.
- Normal problem/formulation validation remains mandatory. Injection is not
  permission to bypass an unsupported pipe model, invalid geometry, earth
  constraints or other ordinary input errors.
- Injection does not imply tracing. With `trace=false`, use the supplied values
  for the call without retaining a blueprint collection in the result.
- With both options enabled, the new trace exposes the blueprints actually used.
  Record whether they were supplied in a simple trace field such as
  `blueprints_reused`; do not introduce a telemetry subsystem.
- The option is execution state: do not put it in physical formulation identity,
  quantity descriptions, comparison identity or scientific model records.

### 4.3 Formulation collections

Each result's trace retains the collection used for that result's formulation.
With no injection, preserve existing sharing between compatible local selections.

With injection, the one supplied collection must be compatible with every
formulation in the call. Validate all of them before emitting results or invoking
`on_result`. Earth/series/reduction variants can share it; incompatible local
shunt selections cannot. For incompatible variants, the caller makes separate
calls with their own traced collections. Do not add a nested cache-routing API.

Preserve result order, callback behavior and stable result types. A compatible
one-design collection can be reused between the two coaxial consumers; general
cross-system partial matching is not part of this first implementation.

## 5. Construction provenance and validation

The current blueprint contains completed numerical data but not enough
construction provenance for safe injection. Add the minimum source identity
beside that data, inside the same owner and object.

The implementation must satisfy these requirements:

1. Capture the physical design state actually consumed by construction:
   geometry, materials, local placement/topology and terminal ownership, plus
   any declaration fields still consumed after architectural cleanup. Include
   scalar representation and the relevant normalized local formula settings.
2. Retain an owned snapshot of those construction inputs, not a live reference
   to mutable source arrays. Mutating the source later must not mutate the
   recorded basis against which compatibility is checked.
3. Compare current physical state against that snapshot without rebuilding its
   numerical equivalents. Use exact state comparison, preserving uncertainty
   source identity and correlations, rather than approximate numerical equality.
4. Do not use only `cable_id`, `system_id`, object identity, terminal counts,
   rounded values, nominal uncertainty statistics or an unchecked hash as proof.
   Equal equivalent annuli also do not prove equal open-wire/tape geometry.
5. Reuse the owner's existing state-comparison methods. Keep the dependency
   selection shared with ordinary blueprint construction so the two paths cannot
   drift. An owner-local validation method on an existing generic is preferable
   to a new validation/cache vocabulary.
6. Retain only construction inputs, not a second full problem, earth model,
   frequency grid, full formulation or workspace. Measure snapshot memory and
   validation cost, including their overhead on the no-reuse path.

Compare normalized choices conservatively. Initially, aliases may require the
same selected identity even when their equations coincide. Reuse across aliases
is permitted only if the cleaned-up formula owner already defines their canonical
equivalence and requested/effective metadata remains truthful. Do not infer
equivalence by evaluating a formula or special-casing names in cache code.

Resolution, tolerance, audit and fallback settings are part of construction
compatibility. Do not infer that an old calculation satisfies a new request
because its tolerance looks tighter or its matrix looks sufficiently close.

Cached fallback remains fallback. Retain its reason and effective model, and
reissue the existing user-visible fallback warning when using that replacement.
A strict boundary request must reject a blueprint constructed under an
incompatible fallback policy. An audited request must not reuse an unaudited
blueprint as though the independent checks ran.

Construction diagnostics describe the retained construction, not work performed
again in the current call. Preserve those diagnostics; use the trace reuse flag
to distinguish injection. Do not claim that the old `solves` count represents
fresh solves on a reused call or erase its construction history by setting all
diagnostics to zero.

## 6. Lifetime, mutation and persistence

- Treat completed blueprint arrays as read-only. Immutable Julia structs do not
  make their contained vectors and matrices immutable.
- Audit all consumers: coefficient scatter/remapping may share C/P matrices,
  but must never factorize or modify them in place. Each calculation still owns
  separate mutable workspace buffers and outputs.
- Snapshot the outer supplied collection for call-local ownership if needed;
  do not copy dense solver state or make a second public prepared object.
- Caller mutation of a retained blueprint, or concurrent source/blueprint
  mutation during a calculation, is unsupported. Do not add per-call checksums
  of every coefficient array to police this contract.
- Holding only the recovered collection must not retain the original result,
  frequency trace arrays, workspace, callbacks or dense boundary factorization.
  Sharing numerical matrices among equivalent blueprint entries is permitted.
- With `trace=false` and no external owner, temporary blueprint/provenance data
  must become collectible after the calculation; no hidden global owner remains.
- Scope reuse to the same live implementation and scalar/precision context.
  Discard blueprints after source changes, including relevant Revise updates, or
  package upgrades. Do not add a process-global revision counter as part of this
  feature. Older or foreign blueprints lacking valid provenance are rejected.

Scientific serialization is not executable checkpoint persistence. Keep runtime
blueprints out of scientific exports, UQ publications and comparison records.
The current `LineParameters` codec selects scientific fields explicitly, whereas
the `CableConstants` codec serializes its details more broadly; both paths and
the Measurements shared-source codecs require review. Omit runtime trace data
deliberately without dropping physical shunt outcomes. Deserialized scientific
results need not provide a reusable blueprint and must not reconstruct one
implicitly or run a solver while loading.

## 7. Parametric studies and uncertainty

UQ currently forwards execution options into each realized problem. Adding the
keyword only at the scalar entry point is not sufficient verification.

- Apply construction compatibility to every materialized problem, not just the
  study's nominal problem or first trial.
- A supplied blueprint may serve a study varying only external earth,
  temperature, placement or other permitted runtime inputs when representation
  and source-state checks actually match for each call.
- Changed local geometry or material realizations require new construction;
  explicit stale injection fails rather than rebuilding behind the caller's back.
- Preserve the Measurements dependency graph. Equal nominal values and standard
  deviations from independent sources are not interchangeable. A deterministic
  blueprint must not remove uncertain geometry/material dependence.
- Nominal/uncertain scalar-type mismatches fail clearly; do not convert stored
  coefficients to manufacture compatibility.
- An injection mismatch is an execution-contract error, not a physical
  `DomainError` eligible for Monte Carlo rejection/resampling. It must propagate
  even when `on_error=:retry` is selected.
- Preserve the existing prohibition on automatic model fallback in UQ and the
  model-coverage check across realizations. Reused fallback cannot evade either.
- Keep trace retention separate from scientific retained-details publication;
  do not accidentally retain a blueprint per sample in every exported record.

No automatic per-realization cache manager is added. Callers must omit injection
for studies whose local construction inputs change.

## 8. Implementation sequence after activation

All items remain pending while this plan is deferred.

1. **Rebaseline and finalize ownership.** Complete the activation gate. Inventory
   blueprint dependencies and consumers after cleanup; settle the minimal owned
   provenance representation and reuse the established comparison contract.
2. **Add construction provenance and compatibility.** Extend the existing
   blueprint constructor and owner-local validation methods. Keep the numerical
   kernels unchanged. Add mismatch tests before exposing the public option.
3. **Retain completed blueprints.** Extend result trace assembly for both
   consumers without reconstructing a blueprint from `LocalCableData`, changing
   result return types, or retaining an entire workspace.
4. **Add validated injection.** Normalize `blueprints=nothing` at the appropriate
   execution owners. Branch at the current flattening boundary: construct when
   absent, validate/use when present. Keep subsequent input/workspace assembly
   common to both paths. Handle scalar and collection entry points together.
5. **Integrate UQ and record boundaries.** Verify per-realization validation,
   failure propagation, result/callback typing, and exclusion of runtime objects
   from all affected scientific serialization paths.
6. **Complete native display and documentation.** Extend bounded, read-only
   blueprint/result displays only as needed. Document C/P units, source ownership,
   validity rules, trace memory costs, same-implementation lifetime, errors and
   both solver examples. No printing or plotting may trigger construction.
7. **Verify correctness and measured benefit.** Run the focused and wider checks
   below, then benchmark with no competing verification processes. Record actual
   performance and limitations before claiming completion.

Expected owner locations, subject to the architecture cleanup:

| Responsibility | Existing owner/files |
|---|---|
| Blueprint state, construction inputs and validation | `src/engine/blueprint.jl`, `blueprint_shunt.jl` |
| Execution option normalization | `src/engine/options.jl`, `cableconstants.jl` |
| Construction/injection decision and result retention | `src/engine/lineparameters.jl`, `cableconstants.jl` |
| System remapping and mutable workspace separation | `src/engine/input.jl`, `admittance.jl` |
| Uncertainty-aware state and per-realization behavior | existing DataModel/Engine state comparison, UQ owners, Measurements extension |
| Scientific serialization boundaries | `src/importexport/uncertainty.jl`, affected extension codecs |
| Bounded inspection | `src/engine/textdisplay.jl` |
| User contract and manual timings | `docs/src/engine.md`, existing disposable benchmark script in `dev/` |

Do not create a new module tree, cache service, preparation wrapper or generic
helper layer. A new file is justified only by an actual responsibility within
the existing owner. No compatibility aliases are required for unpublished API.

## 9. Verification matrix

### Correctness and construction bypass

- [ ] Both consumers: traced and untraced ordinary calculations have equal
  scientific outputs and unchanged public result types.
- [ ] Both consumers: injected and freshly constructed calculations agree for
  annular and boundary models at the same inputs.
- [ ] Changed frequencies, temperature, earth, external placements, connections,
  reductions, transposition and output basis match fresh calculation results.
- [ ] Independent local assemblies and repeated cable designs preserve terminal
  mapping and existing matrix sharing; reordered/missing/extra supplied entries
  fail when they do not match the declared design order.
- [ ] Changed physical geometry/materials with unchanged cable ID and terminal
  counts are rejected. Mutating source arrays after construction cannot mutate
  provenance into falsely accepting stale coefficients.
- [ ] Changed model, admitted boundary law, numerical controls, fallback policy,
  audit request, scalar type or supported precision context are rejected unless
  the declared compatibility contract explicitly permits them.
- [ ] Strict, fallback, audited and unaudited outcomes remain truthful, including
  warning behavior and diagnostic preservation.
- [ ] Direct construction-entry instrumentation proves zero flattening and zero
  boundary solves for compatible injection. Validation failures also perform
  neither. Use test-owned dispatch counters, not elapsed time or source tokens.
- [ ] Compatible formulation collections share the supplied data; an incompatible
  collection fails before callbacks or partial results. Callback order is unchanged.
- [ ] One-design cross-consumer reuse works with matching local choices and type.
- [ ] Ordinary calls without injection retain within-call deduplication and make
  zero boundary solves for default annuli.

### Ownership, inference and interfaces

- [ ] Repeated/concurrent read-only reuse cannot mutate C/P, conductor/layer
  arrays or provenance. Outputs and workspaces do not share mutable buffers.
- [ ] `trace=false` retains no reusable-blueprint trace; `trace=true` identifies
  the actual supplied/constructed values and keeps all existing line trace fields.
- [ ] Keeping the collection alone does not keep frequency traces or dense
  numerical workspaces alive. Check ownership deterministically, supplemented by
  memory measurements rather than fragile GC timing assertions.
- [ ] Warm/default and injected hot paths maintain the established inference
  guarantees; optional execution state does not widen frequency-loop storage.
- [ ] Result displays remain bounded, side-effect-free and owned by the package.
- [ ] FEM rejects `blueprints`; unknown option validation remains intact.

### UQ, serialization and regressions

- [ ] Shared-source LEP reuse preserves sensitivities/covariances; independent
  equal-valued sources are rejected. Deterministic/uncertain mismatches are rejected.
- [ ] MC with unchanged local construction inputs can reuse only when each draw
  matches; changed local draws fail injection without rejection/resampling.
- [ ] Existing UQ fallback and coverage policies cannot be bypassed by injection.
- [ ] Scientific serialization round-trips outputs and physical model outcomes
  with tracing enabled, but does not export or reconstruct reusable runtime objects.
- [ ] Shared-source Measurements codecs preserve scientific correlations while
  excluding blueprint provenance and other runtime trace state.
- [ ] Retained scientific records, reporting and comparisons do not acquire
  blueprint-dependent identity or require an executable checkpoint.
- [ ] Run focused blueprint/shunt, cable constants, formulation-grid, option,
  display, serialization and UQ tests; then quality, Aqua and relevant Gauntlet
  record/dielectric tests. Run the ordinary suite and documentation checks in an
  isolated environment/worktree when necessary to preserve unrelated edits.
- [ ] Preserve independent numerical controls and existing acceptance tolerances;
  do not weaken tests to accommodate reuse.

Use the installed Julia executable and already available environments for
verification. This work does not authorize package updates, a consumer environment
switch, publishing, or regeneration of retained Gauntlet reference attempts.

## 10. Performance validation

Use the existing 18 kV trefoil and 132 kV horizontal cases, 101 frequencies from
0.1 Hz to 10 MHz, fixed formulations and fixed reduction settings. Rebaseline
after cleanup; the earlier 18 kV boundary call took approximately 11.0–11.3 s
and allocated 295.6 MiB including blueprint construction. That is motivation,
not a promised post-cleanup timing or evidence for this unimplemented feature.

Measure separately with timing macros:

1. Cold first call, reporting compilation separately.
2. Warm ordinary default annular compute without tracing.
3. Warm ordinary boundary compute without tracing.
4. First traced boundary call that exposes the blueprints.
5. Warm compatible injected compute with tracing disabled.
6. Warm compatible injected compute with tracing enabled.
7. A compatible changed frequency grid and a collection of earth-return variants.

Warm each path before repeated measurement. Record Julia/BLAS thread settings,
wall time, allocations, allocated bytes, retained blueprint/provenance size and
trace size. Distinguish cumulative allocation from retained/peak memory. If
internal validation is timed separately, report it as diagnostic evidence in
addition to, not instead of, the full public-call cost.

Performance acceptance requires:

- Zero blueprint construction/solve entries on reuse, established by tests.
- A substantial measured reduction in full warm boundary-call time and
  allocations, with validation cheaper than the avoided construction.
- No dense solver matrices/factorizations retained in the reusable payload.
- No material regression in the default no-reuse path. Compare against the
  post-cleanup baseline and investigate repeatable overhead before accepting it.
- Explicit reporting of trace/provenance memory costs and remaining cold latency.

Do not set wall-clock assertions in normal CI, promise the default annular time
for a reused boundary call, or claim that this optimization improves physical
accuracy. If validating/retaining provenance costs too much, revisit the owner
design rather than dropping correctness checks or introducing hidden caches.

## 11. Completion criteria and deferred-state discipline

After explicit activation, completion requires all contract, ownership, UQ,
serialization and performance checks above, updated user documentation, and a
record of actual verification outcomes. Any excluded suite or unrelated baseline
failure must be stated rather than reported as a passing gate.

Until then, this file is the only deliverable for this proposal. It does not
change the existing no-cross-call-reuse behavior, supersede the architecture
cleanup, or authorize implementation to begin.
