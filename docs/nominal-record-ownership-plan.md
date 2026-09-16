# Nominal ownership of options and computation details — amended execution plan

Status: **IMPLEMENTED — verification and pre-existing failures recorded below**.

Recorded: 2026-09-16.

This document replaces the supplied implementation draft. The user explicitly
authorized end-to-end execution on 2026-09-16. That authorization applies to this
ownership migration, not environment changes or other deferred plans.

Principal amendments to the supplied draft:

- Public `options=(...)` shorthand is explicitly preserved only at ingress;
  internal record contracts remain strict (section 5.1).
- Selected-equation controls have an explicit formulation owner, and passive
  inspection is separated from normalization (sections 3.1–3.2).
- Equality/hashing, saved identities, codecs, and diagnostic specialization
  bounds now have concrete requirements (sections 6–7).
- Source loading, all-owner integration coverage, and before/after performance
  verification are explicit completion gates (sections 4.1 and 9).

## 1. Outcome, admission, and non-goals

Replace the three aliases in `src/grammar/types.jl` with exactly three nominal,
immutable records, each storing one named-tuple payload:

- `FormulationOptions`: inputs owned by a selected formulation/equation.
- `ComputationOptions`: inputs owned by a computation or execution backend.
- `ComputationDetails`: supplemental output owned by a completed computation.

The implementation is complete only when producers, owner stages, storage,
consumers, optional integrations, documentation, and tests enforce those roles.
Adding three declarations without migrating the actual paths is incomplete.

### 1.1 Structural-design admission record

This is an ownership-contract refactor, not an optimization or new framework.

| Warrant item | Decision and current evidence |
|---|---|
| Names and kind | The three names above; three immutable record/wrapper types. |
| Owner | Existing `Grammar` module; declarations in `src/grammar/types.jl`. Concrete policy remains with the receiving/producing owner. |
| Required behavior | Reject role substitution at real formulation, execution, and result boundaries, including identical payload schemas. |
| Native alternatives checked | `NamedTuple` already provides payload storage and operations; aliases do not establish nominal role distinction. Native payload operations remain in use. |
| Dependency alternatives checked | No dependency-owned type supplies these package-specific roles. No dependency or replacement collection interface is needed. |
| Existing repository alternatives | Retain `formulation_options`, `computation_options`, `computation_details`, `details`, constructors, native equation dispatch, and existing declaration/serialization interfaces. |
| Semantic delta | The wrapper establishes role, not owner identity, validation state, scientific selection, deep immutability, or cache validity. |
| Current callers and method families | Engine, Earth, Materials, Transforms, ParametricBuilder, UQ, Gauntlet, result consumers, and backend extensions use these aliases or corresponding bare-tuple slots. |
| Invariant and tests | Exact top-level payload storage plus wrong-role rejection at actual owner stages; owner-specific validation and end-to-end conformance tests. |
| Why direct owner code is insufficient | Repeated key checks cannot give otherwise identical named tuples different nominal roles across these existing boundaries. |
| External boundary | Existing ImportExport and optional extensions translate their owned records into existing external representations. No new adapter layer. |
| Compatibility evidence | Preserve identified supported persisted formats and current public keyword convenience. Unreleased internal tuple contracts receive no compatibility bridges. |
| Performance evidence | An isolated Julia 1.12.7 wrapper read inferred `Float64` and allocated zero bytes; this is feasibility evidence, not an end-to-end speed claim. Existing diagnostic type bounds must survive. |
| Removal condition | Not temporary compatibility. There is one final internal record contract. |

No new public action, module, hierarchy, schema registry, owner registry,
validation framework, or module-level private generic is admitted by this plan.
Use existing owner interfaces and direct native operations. Do not force the
code into a new Template Method pipeline: tighten the already owned stages.

### 1.2 Preserved behavior

Preserve numerical algorithms, equations, units, terminal/frequency order,
scientific selection, option keys and values, merge precedence, output basis,
default values, fallback policies, and failure semantics. Representation and
role-bearing interface changes specified here are deliberate.

In particular:

- `:default` remains an explicit route to the selected concrete formulation.
  Empty options, available numerical sections, and wrapper types do not choose
  physics or replace that route.
- Preserve earth `air`/`earth`/`mixed` and internal impedance
  `inner`/`outer`/`transfer` meanings, including consumer-specific projection
  and rejection of unused explicit sections.
- Do not reintroduce callable formulation overrides or a `hooks` field. This
  document uses **owner stage/interface method** for existing dispatch seams.
  Existing execution callbacks such as `on_result` are a separate supported
  contract and are not removed by this work.
- Do not implement blueprint reuse, new caches, solver tuning, numerical
  safeguards, dependency updates, or unrelated architecture cleanup.
- Preserve unrelated worktree changes, including `dev/inspect_saved_gauntlet.jl`
  if it is still modified when implementation starts.

## 2. Reinspect the live checkout and establish the baseline

Before implementation, run read-only discovery:

```sh
git status --short
rg -n 'FormulationOptions|ComputationOptions|ComputationDetails' src ext test docs gauntlet dev
rg -n 'options::|details::|parameters::|formulation_options|computation_options|computation_details' src ext gauntlet
```

Follow actual callers, not only annotations: the current aliases make many
role-bearing arguments indistinguishable from ordinary `NamedTuple` data.
Reconfirm file locations rather than treating this plan as a source manifest.

Record, in the implementation work record, a finite migration inventory with:

1. Public ingress and the fixed owner of each explicit option group.
2. Producers, raw versus resolved flow, default/validation stages, and the
   consumers of projected numerical sections.
3. Stored options, completed supplemental output, and deliberate type-erasure
   boundaries in diagnostic storage.
4. Inspection, identity, equality, display, serialization, cache/signature,
   worker transport, report, and callback consumers.
5. The concrete test items that exercise each path, including optional owners.
6. Existing failures, unavailable external tools, and numerical/performance
   baselines captured before editing.

Use current code, maintained architecture documents, and tests as evidence.
Do not recover old execution paths from historical prompts. A zero-ambiguity
result alone is not proof of coherent ownership.

Capture representative cold-load/first-call and warm-compute measurements
under the same executable, flags, thread counts, and environment that will be
used after the change. Record `VERSION`, `Base.active_project()`, and the loaded
package source. Do not edit `Project.toml`, manifests, Julia channels, startup
files, or the user's IDE environment to make this refactor pass.

## 3. Fix the ownership map before changing signatures

Classify by responsibility, not by field name, storage location, lifetime, or
the spelling of the generic currently used to normalize/read it.

| Record or boundary | Final role | Existing owners to cover |
|---|---|---|
| High-level formulation controls, including reduction/transposition and FEM physics selection | `FormulationOptions` | Engine line-parameter/cable-constant formulation constructors and FEM formulation owner |
| Controls of a selected native equation, including quadrature, boundary discretization/audit, and modal algorithm controls | `FormulationOptions` | All admitting native formula families, their declarations, resolved selections, and equation bindings |
| Calculation/backend invocation controls, including tracing, output basis, verbosity, solver process/mesh execution controls, and modal action tolerance | `ComputationOptions` | Coaxial, FEM/Gmsh, PSCAD, Transforms |
| Controls owned by composite calculations, including LEP retention, Monte Carlo trials/seed/error policy, and Combinatorial retention | `ComputationOptions` | UQ and ParametricBuilder, even though the record is stored in an `AbstractFormulation` subtype |
| Core calculation options retained only to forward to the eventual backend | `ComputationOptions` | `ParametricProblem`, Gauntlet calculation declarations, collection execution |
| Supplemental output retained by completed calculations | `ComputationDetails` | Core, modal, parametric, uncertainty, and external-backend result owners |
| Physical `parameters`, assumptions, constitutive state, geometry, axes, primary scientific/statistical products | Existing domain representation; not one of these wrappers merely because it is a tuple | Their current domain owners |
| Nested groups inside one owned record | Ordinary named tuples/data unless passed as an independently owned record | The enclosing owner |
| Publication/comparison settings and native renderer keyword attributes | Existing presentation/report representation | ReportBuilder and plotting extensions; update access to wrapped scientific records without relabelling presentation settings |

The formula-family inventory includes frequency-dependent earth,
equivalent-homogeneous earth, earth impedance/admittance, internal impedance,
insulation impedance/admittance, semiconductor admittance, shunt geometry,
temperature-dependent material laws, modal formulations, and any other live
registered family discovered at activation. Families with no option-bearing
state, such as a currently stateless pipe formula, must not acquire an empty
stored record solely for symmetry; test their admitted/rejected input contract.

### 3.1 Resolve the equation-control mismatch

`src/grammar/formulas.jl` currently declares and normalizes selected-equation
numerical controls through `computation_options(binding::FormulaMethod, ...)`.
Under this plan these are formulation-owned controls.

Migrate the complete equation-owned method family and all its callers to the
existing `formulation_options` generic, using the existing `FormulaMethod`
identity and section dispatch. Preserve its numerical defaults, per-case
projection, configured-section bookkeeping, validation order, and error
behavior. Update explicit imports and extension guidance with those methods.

Do not leave old/new forwarding methods. Do not convert a record from one role
to another midway through construction. `FormulaDefinition.options` and the
resolved formula's corresponding option record have the same role; definition
construction remains passive. Physical `parameters` stay separate.

Top-level record transport between owners/stages uses the wrapper. Existing
section-level methods operating on the *same owner's* projected nested tuples
may continue using those tuples; this does not admit a top-level tuple bypass.
Rewrap only when establishing an actual independent role-bearing boundary.

### 3.2 Separate normalization from passive inspection

Current `formulation_options` overloads also collect retained parameters,
numerical options, and equivalent-earth declarations for descriptions and
identity. Some UQ overloads return computation-owned options under that name.
Do not preserve those mixed contracts by changing only their annotations.

The target contracts are:

- Owner-resolution methods receive and return their corresponding role.
- Any retained `formulation_options(value)` read returns only that value's
  actual `FormulationOptions`, without validation or a new option envelope.
- A read of stored computation controls, where needed, belongs to an existing
  `computation_options` method for the owning completed declaration; it must
  not masquerade as formulation options.
- Mixed scientific declarations/provenance remain ordinary declaration
  records exposed through the existing owner-defined `NamedTuple`, `pairs`,
  `formula_id`, `description`, and persistence interfaces. They are not
  execution input records and must not be relabelled as such.

Migrate identity/display consumers to the appropriate existing declaration
methods. Remove mixed-purpose retained-record fallbacks that call live option
normalizers. Do not add a parallel `*_settings`, `*_controls`, or generic
inspection framework to compensate.

Native and saved identities must retain their existing scientific content:
owner scope, routes, identifiers, physical parameters, numerical controls,
equivalent-earth ordering, and relevant composite execution policy. Reading a
saved declaration must not resolve today's defaults, construct a live solver,
load an optional backend, or turn an unavailable identity into an empty one.

## 4. Declare exactly three records

Keep declarations in `Grammar`, preserve existing root reexports, and update
explicit imports. Do not create a new module or duplicate any declaration.

```julia
struct FormulationOptions{NT<:NamedTuple}
    data::NT
    FormulationOptions(data::NamedTuple) = new{typeof(data)}(data)
end
FormulationOptions(; kwargs...) = FormulationOptions((; kwargs...))

struct ComputationOptions{NT<:NamedTuple}
    data::NT
    ComputationOptions(data::NamedTuple) = new{typeof(data)}(data)
end
ComputationOptions(; kwargs...) = ComputationOptions((; kwargs...))

struct ComputationDetails{NT<:NamedTuple}
    data::NT
    ComputationDetails(data::NamedTuple) = new{typeof(data)}(data)
end
ComputationDetails(; kwargs...) = ComputationDetails((; kwargs...))
```

Document these using active repository docstring conventions. Each record has
exactly one field. Its type parameter describes storage, not scientific or
backend identity. No per-owner/schema record structs, abstract hierarchy,
owner tags, validity flags, phase parameters, or schema metadata are allowed.

Constructor requirements:

- Accept a positional named tuple, keyword construction, and the empty keyword
  case. Do not add a redundant zero-argument method.
- Preserve the exact supplied payload type, value types, key order, nesting,
  and mutable-object identity. Do not copy, sort, normalize, promote, or coerce.
- Reject positional dictionaries, ordinary tuples, arbitrary objects, and
  owned records, including same-role records. Public ingress preserves an
  already correctly wrapped value directly instead of invoking its constructor.
- Do not provide explicitly parameterized construction such as
  `FormulationOptions{NamedTuple}(payload)` or implicit `convert` methods.
- Do not validate owner-specific content or recursively wrap nested groups.
- Construction establishes role only. The same role type may carry caller
  syntax or resolved values; it is not a proof of normalization or validity.

### 4.1 Resolve the existing declaration-order dependency

The current root `src/formulas.jl` defines `FormulaDefinition` before
`Grammar` loads; `Grammar` then imports root `FormulaMethod` and
`FormulaDefinition`. Merely changing `FormulaDefinition`'s option bound to
`FormulationOptions` creates a load-order problem because that record belongs
to `Grammar`. Include order is therefore part of this migration, not a reason
to weaken the stored role or introduce late binding.

Use this bounded ownership correction:

- Place the existing shared formula primitives `FormulaDefinition` and
  `FormulaMethod`, with their essential construction/call invariants, in the
  existing `Grammar` owner. Declare the records before types that use them.
- Keep root-owned generic declarations such as `description`, `formula`, and
  `formula_id` available before `Grammar` loads, through the existing root
  interface declaration responsibility. Keep the generic functions themselves
  single-owned; do not redeclare competing functions in `Grammar`.
- Remove the now-unnecessary parent imports of these two primitive types from
  `Grammar`. Root public names continue through explicit imports/reexports of
  the same bindings; callers retain `LineCableModels.FormulaDefinition` and
  `LineCableModels.FormulaMethod`.
- Load remaining root formula construction/description methods after their
  required grammar types and bindings exist. Update includes and imports in
  dependency order, without moving unrelated scientific implementations.

This relocates existing shared declarations; it does not admit additional
types, a new module, or another formula interface. Do not use duplicate type
definitions, `@eval`, repeated includes into an already closed module, abstract
option fields, or old/new forwarding paths to work around the dependency.
Recheck type-qualified persistence/display consumers so the internal module
placement does not silently change supported saved identifiers. Verify clean
loading, explicit imports, and public type-binding identity after the move.

## 5. Public convenience, internal enforcement, and sequencing

### 5.1 One public ingress contract

Preserve existing public named-tuple convenience, including:

```julia
compute(problem, formulation; options=(trace=true,))
compute(problem, formulation; options=ComputationOptions(trace=true))
```

Existing `options=(...)` formulation constructors similarly admit a named tuple
or the correctly owned record. Existing grouped keyword syntax may materialize
its group with `(; kwargs...)` and construct the appropriate record.

This is an explicit amendment to the original draft's blanket prohibition on
bare tuples in public `options` slots. It is one supported boundary grammar,
not a compatibility execution path. The named tuple is wrapped immediately;
both spellings enter exactly the same owner-resolution and execution path.

At each existing public ingress:

1. The positional problem/formulation/backend contract establishes the role.
2. Accept only a named tuple or the correct role. Reject wrong-role records,
   dictionaries, and arbitrary objects without inspecting keys to guess intent.
3. Construct the wrapper for a named tuple; retain an existing wrapper as-is.
4. Pass that wrapper to the existing designated owner stage.

A narrow union or local two-case representation check is permitted **only at
these enumerated public boundaries**. Do not add a generic wrapper/fallback
helper or retain duplicate positional tuple overloads on internal constructors.
Use role-appropriate empty wrapper defaults in the migrated API and storage.

Do not add convenience to unrelated APIs, and do not accept bare tuples as
completed `ComputationDetails` constructor arguments merely because options
have a public shorthand. Result/detail-bearing constructors use the record;
call sites and fixtures migrate explicitly. Codecs construct it at their owned
external decoding boundary.

Keyword types do not select Julia methods. Use one keyword ingress for each
existing positional dispatch signature, and enforce internal distinctions on
positional owner-stage arguments. Test keyword rejection separately from
missing-method rejection.

### 5.2 Internal contract and normalization lifecycle

After public ingress, role-bearing arguments, returns, stored fields, and
cross-owner forwarding use only their assigned wrappers. Internal tuple/wrapper
unions, cross-role conversions, and generic tuple fallbacks are prohibited.

Preserve the existing order:

1. Receive caller configuration at its public/definition boundary.
2. Resolve defaults and validate at the existing owning formulation or
   computation stage when that owner is known.
3. Forward the resolved record into existing preparation/execution.
4. Produce supplemental output at the existing completion boundary and store
   it on the actual completed result.

This list describes ownership obligations; it does not authorize a replacement
pipeline or additional preparation stage. Each existing action keeps its own
sequence and result type.

Normalize at the established owner invocation, not at every forwarding hop or
frequency sample. Composite execution may legitimately invoke an inner solver
per design/trial: do not use this refactor to cache or move those invocations.
Record and test the expected counts at each owner rather than imposing one
global normalization count on all composites.

`ParametricProblem` and Gauntlet can retain caller `ComputationOptions` until
the backend is known. Merely receiving a wrapper must never skip the backend's
validation. Conversely, already resolved options must not be fed back into a
caller-syntax normalizer. For example, the current engine resolves Boolean
`trace` to `Val(trace)`; wrappers must not cause duplicate normalization to
reject the resolved representation. Keep this distinction in the existing
call graph, not new flags, tagged types, key probes, or catch-and-retry paths.

Within an owner stage, use `.data`, preserve that owner's defaults, precedence,
canonical resolved key order and nesting, and reconstruct the **same role**
only when it returns revised content. Wrong keys/values fail with the relevant
owner identified before expensive execution or external process launch.

Do not add universal required-stage defaults. Missing implementations remain
interface errors; legitimate empty options for a supported owner have explicit
owner methods. Arbitrary wrapper schemas across owners do not imply arbitrary
keys are accepted by any particular owner.

### 5.3 Supplemental output and output ownership

`computation_details` produces `ComputationDetails` through the existing
result/completion owners. `details(result)` returns that owned record wherever
the existing interface promises supplemental output. Cross-owner consumers
use that semantic accessor, then `.data`; direct fields remain appropriate
inside the result owner.

Use `ComputationDetails()` for successful empty supplemental output. Preserve
failure, missing, unavailable, and incomplete states as distinct states.
Preserve trace selection, callback behavior, and all existing retained fields.
Do not introduce trace fields for the separately deferred blueprint-reuse plan.

Construction/preparation diagnostics stored on blueprints or workspaces are
not automatically completed-computation records. Their existing bounded
representation may remain nested inside the final details payload.

Immutability is shallow. Preserve existing scratch/output-buffer separation,
array identity where promised, and concurrency rules. Wrapping must not expose
a scratch buffer that later computations overwrite, copy output arrays merely
to appear immutable, or imply that mutation is now thread-safe.

## 6. Payload operations, value semantics, and consumers

### 6.1 Explicit payload access

Use native operations on the payload:

```julia
options.data.integration.rtol
get(options.data, :trace, false)
keys(options.data)
pairs(options.data)
isempty(options.data)
```

Owner-local merges/projections use those payloads and reconstruct the same role
where required. Numerical kernels continue receiving the scalar/array arguments
they actually need. If a kernel receives a full owned option record, do not
silently revert that record to unowned top-level tuple transport.

Do not add tuple impersonation through `getproperty`, `propertynames`,
`getindex`, `iterate`, `keys`, `values`, `pairs`, `get`, `haskey`, `length`,
`isempty`, or `merge`. Do not add wrapper `NamedTuple`/`convert` routes,
`unwrap`, `payload`, or formatting/accessor helpers. `.data` is sufficient.

At existing external-library boundaries, explicitly select and unpack the
arguments required by that dependency. Do not pass wrappers to APIs that do
not understand them or indiscriminately splat unrelated options.

### 6.2 Deliberate equality and hashing

Do not inherit accidental object-identity comparison. The audit's Julia 1.12.7
probe showed that independent array-valued named tuples compare structurally,
whereas the proposed bare wrapper declarations do not automatically do so.
Current formula identity/description code uses `isequal` on retained controls.

Implement the following value contract through native Base methods owned by
the three record types:

- Same-role `==` delegates to payload `==`, preserving Julia's value semantics
  (including a possible `missing` result), not wrapper or payload-type identity.
- Same-role `isequal` delegates to payload `isequal` and returns a Boolean.
- Different roles, and a record versus its raw payload, are not interchangeable
  and compare unequal. Do not use an unrestricted catch-all to enforce this.
- Hash the semantic role and the payload consistently with `isequal`. Do not
  include the concrete payload parameter in a way that distinguishes payloads
  that Julia considers `isequal`.

These native value-semantic methods are deliberate exceptions to the ban on
tuple impersonation; they do not establish a collection interface. Write the
small role-specific methods directly, without macros or a record hierarchy.

Use the existing owning source file for Base behavior; create no global bucket
or new submodule. Preserve owner/formulation identity in existing scientific
keys: a role wrapper is not a sufficient key for a calculation or blueprint.
Do not use human display or process-dependent object hashes as portable keys.

Mutable payloads are not safe dictionary keys if their content changes after
insertion. This plan does not introduce caching, deep freezing, or copying.
Preserve existing snapshot/signature ownership; callbacks and opaque values
must not silently acquire a portable identity or serialization promise.

### 6.3 Display, provenance, persistence, and adapters

Update each actual consumer rather than compensating with wrapper forwarding.
Follow existing `show` conventions for readable role-aware record display;
do not print an enormous schema type or pretend the object is a bare tuple.
Scientific descriptions and comparisons continue through their existing
owner interfaces and must preserve the native/saved identity contract in 3.2.

Cover at least:

- `src/formulas.jl`, `src/grammar/formulas.jl`, and owner declaration methods.
- `src/importexport/serialize.jl`, `deserialize.jl`, and `uncertainty.jl`.
- Core, modal, parametric, and uncertainty result read/publication paths.
- Gauntlet declarations, signatures/reuse decisions, callbacks, saved results,
  worker transport, and comparison/display consumers.
- Gmsh/FEM and PSCAD execution/result transport and optional report consumers.

Preserve supported external record shapes where the enclosing owner and field
already determine the role; do not leak a new `data` nesting level into those
formats. Decode using that known owner and explicitly construct its wrapper.
Never infer a role by recognizing payload keys or evaluating saved code.

Distinguish passive historical records from executable input: reading a saved
result does not revalidate old numerical settings against current defaults.
Admitting a saved declaration for a new computation still goes through the
normal live owner validation.

Identify supported versioned artifacts and fixtures before editing codecs.
Where an existing generic format actually supports standalone owned records
without enclosing role information, use that format's existing explicit tag
mechanism, with tests; do not introduce a parallel codec framework. Record any
necessary version decision rather than silently altering the wire contract.

Same-version worker/Julia transport must round-trip the records using the
existing transport. Do not claim unsupported cross-version executable
checkpoints are portable. Preserve existing cache invalidation semantics;
do not delete saved attempts, rewrite user artifacts, or silently recompute a
saved result to mask a decoding failure.

## 7. Concrete option storage and bounded diagnostic specialization

Use inferred container parameters such as `O<:FormulationOptions`,
`O<:ComputationOptions`, or `D<:ComputationDetails`, with fields `options::O`
and `details::D`, on the existing concrete hot paths. Reuse current parameters
instead of adding redundant schema parameters. Callers do not spell schemas.

Do not replace hot concrete storage with `Any`, an abstract wrapper field, or
internal tuple/wrapper unions. Do not encode payload values, arbitrary symbols,
thread counts, tolerances, paths, or array lengths as new type parameters.
Preserve existing bounded `Val` choices; introduce no new specialization policy.

**Exact top-level payload storage does not require recursively concrete
diagnostics.** Preserve the deliberate limits already present in:

- `Engine._retained_details` in `src/engine/lineparameters.jl`.
- `CableBlueprint.shunt_details` and `LocalCableData.shunt_details` in
  `src/engine/blueprint.jl`.
- `CableConstants.details` and relevant result-space element types.

For example, preserve this bound before wrapping:

```julia
retained = NamedTuple{(:shunt_model,), Tuple{NamedTuple}}((shunt_details,))
details = ComputationDetails(retained)
```

The wrapper retains `typeof(retained)` exactly; the nested `shunt_model` field
remains intentionally bounded. Do not reconstruct it as
`(; shunt_model=shunt_details)` and silently expose the geometry-dependent
schema to result specialization. Preserve the trace/no-trace envelopes and
their existing field meanings.

Before replacing previously unparameterized diagnostic fields, establish a
stable owner-produced envelope using the existing key contract and appropriate
nested field bounds. Preserve the supported outer grammar without introducing
per-geometry detail types or changing public field names to hide type growth.
If existing behavior has truly variable outer schemas that cannot meet these
requirements, stop and present that concrete conflict; do not invent a new
result/details type family or silently erase hot types.

The acceptance criterion is concrete numerical/option paths **and** bounded
non-numerical specialization. It is not a global ban on `NamedTuple`-typed
nested diagnostic fields. Test inference and allocations on the actual
containers and bundled/UQ paths; constructor inference alone is insufficient.

## 8. Execute one finite migration after authorization

Implement in this order, without leaving parallel contracts in the final tree:

1. Establish the inventory, ownership map, raw/resolved flow, and numerical,
   identity, persistence, and performance baselines.
2. Add the three types, constructors, native value semantics, and focused
   tests in their existing owners; resolve the declaration order in 4.1.
   Do not temporarily commit live alias bridges.
3. Migrate equation/formulation ownership and passive inspection contracts,
   including all formula declarations, native method families, and imports.
4. Migrate public ingress, execution options, stage signatures, and option
   storage for core calculations and composite owners.
5. Migrate supplemental-output production, stable diagnostic envelopes, result
   storage/read interfaces, and buffer-ownership-sensitive consumers.
6. Migrate codecs, identity/display/reporting, Gauntlet, optional backends,
   worker boundaries, maintained examples, and documentation.
7. Run focused and broad verification, inspect the complete diff, and remove
   residual obsolete routes. Finalize only with all required work accounted for.

Apply each step across its real producers and consumers, not by replacing every
`NamedTuple` match. During a local edit sequence intermediate breakage is not a
reason to add permanent bridges. Do not alter numerical bodies except for
necessary record access/projection changes.

No new `_normalize_options`, `_prepare_options`, `_resolve_formulation_options`,
`_unwrap_options`, `_wrap_details`, `_prepare_computation_details`, renamed
equivalents, or catch-all parser/normalizer is allowed. Existing owner-local
algorithms and genuine reusable stages remain; do not delete them solely
because their names contain `prepare` or `details`.

## 9. Verification and acceptance requirements

Extend the repository's existing tests and conformance infrastructure. Do not
create another test harness or production instrumentation API.

### 9.1 Record construction and native semantics

For all three types, test:

- Positional, keyword, and empty construction; exactly one field named `data`.
- Exact `fieldtype(typeof(record), :data) === typeof(payload)` and payload
  identity, including nested arrays and a payload key itself named `data`.
- Multiple schemas, preserved key order/value types/nesting, and directly
  written representative `@inferred` constructor checks.
- Rejection of dictionaries, tuples, arbitrary objects, same/wrong-role records,
  and explicitly widened constructor routes.
- Distinct nominal types; records are not named tuples.
- Independent but equal array-valued payloads compare equal within a role;
  nested values, `missing`, NaN/signed-zero behavior, and equal payloads with
  different concrete numeric types follow the underlying Base semantics.
- `isequal` implies equal hash values; cross-role/raw-payload inequality is
  retained. Do not require unequal values to have collision-free hashes.
- Array mutation remains visible through the payload; no copying/freezing or
  immutability/thread-safety promise is introduced.
- Already bounded diagnostic payloads retain their bound after construction.

### 9.2 Actual boundaries and owner semantics

Test real owner stages, not only toy wrapper dispatch:

- Enumerated public `options=(...)` and explicit-record calls produce equivalent
  resolved configurations/results and traverse the same stages.
- Internal top-level record slots reject bare tuples and wrong-role records,
  even with identical keys. Public wrong-role inputs fail before side effects.
- Defaults, overrides, merge precedence, resolved key ordering, nested section
  projection, and unsupported-key/value diagnostics are preserved per owner.
- Arbitrary construction is not mistaken for owner validation. A record valid
  for one backend is not assumed valid for another because its role matches.
- Explicit default routes remain routes; wrapper/schema/default-section checks
  do not choose scientific formulations or alter supported case coverage.
- Test-owned stage instrumentation confirms existing ordering and resolution
  counts, including composite forwarding and raw/resolved `trace` handling.
- Required implementations remain required. Native equation dispatch, not
  registry mirrors or user callable overrides, owns scientific extension.
- A test-only native formula/backend adds ordinary owner methods without edits
  to a central schema map, selector switch, or production instrumentation.

Check live and retained/saved inspection independently: identity and description
must not call option resolution, solver construction, or unavailable extensions.
Retain non-default physical parameters, numerical settings, ordered
equivalent-earth reductions, and relevant UQ/composite policy in comparisons.

### 9.3 Integration coverage matrix

| Path | Required evidence |
|---|---|
| Every live native formulation family | Owner-conformance/default/override/rejection tests; numerical fixtures for affected evaluations |
| Analytical line parameters and cable constants | Default/non-default formulation, trace on/off, actual result/detail types, units/order, repeated calls |
| Earth routes and internal surfaces | Air, buried, mixed and inner/outer/transfer routing/projection; explicit and default selection equivalence |
| Boundary shunt | Selected/effective model and fallback/error diagnostics unchanged; bounded result typing across supported geometries |
| Modal computation | Formula options versus action options; forward/inverse results and retained operators/details |
| Combinatorial, LEP, Monte Carlo | Own versus forwarded options; per-point/trial detail retention; seeded behavior and failure policy unchanged |
| Gmsh/FEM and PSCAD | Option normalization, offline adapters/workers, output details, and available executable integration tests |
| Gauntlet | Declaration/forwarding/callback paths, saved result reads, signatures/reuse behavior, comparison identity |
| ImportExport | Supported native/saved round-trips, explicit role reconstruction, numeric/array types and unavailable states |
| Reports and optional plotting/export | Reading owned details without changing presentation/native attribute tuples or forcing weak dependencies |

Keep external executables optional for core tests. Use existing offline fixtures
where available; report live tests unavailable separately, not as passes.
Verify core-only loading does not import optional Makie, Excel, or FEM packages.

### 9.4 Performance and numerical gates

Use the same fixtures, runtime, environment, flags, and threads before/after.
Run timings in functions and separate first compilation from warmed execution.

- Preserve existing numerical tolerances and allocation ceilings; do not loosen
  them or replace expected values to accommodate this representation change.
- Check direct constructor/stage inference and actual core, bundled,
  parametric, and uncertainty result inference under representative settings.
- Test that geometry-dependent nested diagnostics do not create unexpected
  result-type variants. Keep the existing legitimate trace/scalar distinctions.
- Measure warm allocation count/bytes and runtime; report first-load and
  first-compute time separately. Investigate reproducible regressions instead
  of asserting that an immutable wrapper must be free.
- Preserve configured owner output shapes and callback typing without sorting
  caller payloads in constructors or specializing on arbitrary values.
- Repeated compute calls must not require environment changes or recreate
  normalization work at new frequency-loop boundaries.

There is no promised speedup or recompilation cure. The required outcome is
role enforcement without a numerical, inference, allocation, or demonstrated
latency regression. The isolated wrapper probe is not a substitute for this
evidence, and wall-clock noise is not itself proof of regression or success.

### 9.5 Architectural checks

Extend existing scoped quality checks and ambiguity checks to establish:

- No live alias remains and no duplicate role declaration exists.
- No internal top-level tuple bypass, wrong-role fallback, or implicit role
  conversion survives. Public-only ingress allowances are enumerated and tested.
- No old equation-normalization route remains under `computation_options`.
- No passive saved-record inspection falls into a live normalizer.
- No tuple-impersonation or generic unwrapping framework was introduced.
- Native equality/hash/show methods satisfy their stated purpose; they are not
  accidentally prohibited by overly broad checks against all Base extensions.
- Concrete hot storage and deliberately bounded diagnostic fields coexist.
- No new callable formulation overrides, owner/schema registries, phase flags,
  framework helpers, dependencies, or optional-package imports enter core.
- Ordinary scientific/presentation named tuples and same-owner nested sections
  remain legal. Scope checks by role and call path, not identifier substrings.

### 9.6 Run and report the existing harness

The current entry point is `test/runtests.jl`, using
`test/support/runner.jl`. It accepts name/path substrings, `tag:...` selectors,
and `--list`; recheck at activation. The ordinary default run currently excludes
quality, Aqua, and certain external/gauntlet categories, so run relevant explicit
selectors as well. Do not invent command-line options or treat the ordinary run
as coverage of every affected owner.

Select the actual available test environment without mutating the user's daily
environment. Record the executable, project, selectors, and exact commands in
the execution report. Run package loading/precompilation checks, focused record
and owner tests, quality/Aqua/import checks, affected integration/extension
tests, and the applicable ordinary regression suite. Never claim a test run
from listing its items.

Report pre-existing failures separately with baseline evidence. Unavailable
runtime integrations remain unverified; do not mark them passed or broaden the
refactor to repair unrelated faults. Any new failure or unexplained regression
on an available required path prevents declaring implementation complete.

## 10. Documentation, final review, and handoff

Update maintained API documentation, docstrings, extension examples, tests,
Gauntlet/dev call sites, and result-access examples to describe:

- The ownership table and selected-equation versus execution distinctions.
- Public named-tuple convenience versus strict internal record contracts.
- Positional/keyword/empty record construction and explicit `.data` access.
- Owner-specific validation, passive inspection, and the raw/resolved lifecycle.
- Structural equality/hash behavior and shallow immutability limitations.
- Preserved persistence contracts and explicit reconstruction at codecs.
- Concrete hot fields and intentionally bounded diagnostic envelopes.
- Scientific extension through native methods, without callable override hooks.

Remove alias/tuple-transparency claims and stale method names. Do not present
these records as a solver optimization, an immutable cache snapshot, a schema
validator, or an owner/phase identity. Keep deferred optimization documents
deferred; update only example syntax if required by this completed contract.

Finally inspect the diff and rerun scoped discovery. Check stored fields,
defaults, imports, tuple-style operations on wrappers, legacy method routes,
public-ingress allowances, serializers, and new helper names. Preserve unrelated
edits and do not commit unless the user separately requests a commit.

The implementation handoff must report:

1. Declaration locations and migrated owning modules/interfaces.
2. Public/internal boundary decisions, equation-control migration, and removal
   of mixed inspection/normalization behavior.
3. Equality/persistence treatment and preserved diagnostic specialization bounds.
4. Removed obsolete paths, with no internal compatibility bridges retained.
5. Exact verification commands/results, before/after performance evidence,
   pre-existing failures, and unavailable or blocked integrations.
6. Any remaining concrete limitation; no unexecuted test is reported as passing.

Execution is complete only when these contracts are implemented and verified
across the affected current paths.

## Implementation and verification record — 2026-09-16

### Delivered contract

- `src/grammar/types.jl` owns the three immutable nominal records, their exact
  payload constructors, and the relocated `FormulaDefinition`/`FormulaMethod`.
  Root public bindings still identify the same types. `src/grammar/base.jl`
  supplies role-specific equality, `isequal`, hashing and display only.
- Native equation controls now use `formulation_options` and
  `FormulationOptions` across the material, earth, impedance/admittance, shunt
  and modal families. Backend and composite execution controls use
  `ComputationOptions`; supplemental results use `ComputationDetails`.
- Public `options=(...)` convenience and whole-option grids remain supported.
  Internal role-bearing inputs reject bare tuples and other roles. Constructors
  establish ownership, not validation status; existing owner stages still own
  defaults, validation and resolution order.
- Passive native/saved inspection uses the existing owner-scoped `pairs`,
  `NamedTuple`, identity and description methods. It neither normalizes stored
  controls nor reconstructs equivalent-earth declarations. UQ execution controls
  are no longer exposed as formulation-owned controls.
- Engine, Transforms, ParametricBuilder, UQ, Gmsh/FEM, offline PSCAD, Gauntlet,
  reports, plotting consumers and codecs use the corresponding records.
  Known saved envelopes still contain plain payloads, without a new `.data`
  wire layer. Standalone portable records use the existing tagged codec.
- Concrete result/option fields coexist with bounded shunt diagnostics.
  Known native decoders restore that diagnostic bound, including retained
  LEP/MC details and saved Gauntlet calculations. Scientific arrays, primary
  statistics, physical parameters and presentation settings remain native data.
- No tuple-forwarding methods, implicit role conversions, alias bridges,
  owner/schema registries, wrapper helper framework, new dependencies, or
  blueprint-reuse optimization were introduced. Structural-design guidance kept
  the changes in the existing owning stages and codecs.

### Verification

Julia 1.12.7, one thread, `--startup-file=no`. Iterative runs used
`--compiled-modules=existing`. The ordinary test command was:

```sh
JULIA_LOAD_PATH='/home/amartins/Documents/KUL/LineCableModels/test:/home/amartins/Documents/KUL/LineCableModels:@v1.12:@stdlib' julia --startup-file=no --compiled-modules=existing test/runtests.jl
```

For tests needing writable caches/artifacts, only the process environment changed:

```sh
JULIA_DEPOT_PATH='/tmp/lcm-owned-record-depot.YMVUek:/home/amartins/.julia' JULIA_PKG_OFFLINE=true JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_LOAD_PATH='@:/home/amartins/Documents/KUL/LineCableModels/test:/home/amartins/Documents/KUL/LineCableModels:@v1.12:@stdlib' julia --startup-file=no --compiled-modules=existing test/runtests.jl tag:aqua
```

Gauntlet/FEM used a separate first depot,
`/tmp/lcm-owned-record-gauntlet-depot.drP8wB`, with the same remaining settings.
All dependencies were already installed; the user's project, manifests, Julia
channels, startup files and IDE configuration were not changed.

The following are executed results, not test listings. Counts overlap between
focused runs and must not be added together as a distinct-test total.

| Verification | Result |
|---|---|
| Final `grammar/owned_records formula_contracts formula_routes` | 1,154 assertions passed, 11 items |
| `tag:quality` including explicit imports and native catalogue ownership | 1,975 assertions passed, 8 items |
| Final `tag:aqua` | All 11 checks passed, including persistent-task/precompilation hygiene |
| `tag:core_only` | 36 assertions passed, 6 items; optional packages remain unloaded |
| `formula_architecture grammar/owned_records` | 561 assertions passed, 9 items |
| `retained_earth_formulas unit/engine/transforms` | 398 assertions passed, 9 items |
| `unit/parametricbuilder/formulation_grid unit/textdisplay/workflows` | 559 assertions passed, 5 items |
| `grammar/owned_records report_result_protocol internal_shunt_measurements` | 223 assertions passed, 6 items, including bounded persistence and correlated sensitivities |
| Native/saved descriptions and real reports | 9,909 assertions passed; subsequent rerun also covers historical unknown controls/order |
| `tag:fem_numerical fem_electrodynamics` | Real GetDP bare/insulated-wire solves: 164 assertions passed |
| `optional real GetDP multi-frequency scan` | 77 assertions passed, including batched retention and reuse |
| Gauntlet recovery/checkpoints/result spaces/reference replay | Fresh `tag:gauntlet_toolkit result_space_tests numerical_reference_tests checkpoint_tests selective_execution_tests`: 201 assertions passed, 9 items |
| Gauntlet spectral reports/saved comparisons/releases | Fresh `tag:gauntlet_toolkit spectral_reporting_tests summary_tests vault_tests`: 185 assertions passed, 5 items |
| External-engine documentation examples | Executed with constructor rejection and role-return assertions |

The full ordinary sweep ran all 300 maintained items in 127 files. Its initial
20,840 passes, 11 failures and 68 errors include tests captured before migration
corrections. Every migration-related failure was covered by fresh passing
file/owner reruns; the two existing CIM failures below remain. Likewise, the
43-item Gauntlet sweep was followed by fresh recovery/reporting reruns for its
old mock-backend, tuple-constructor and saved-payload assertions. These initial
non-green sweeps are not reported as green full-suite runs.

Offline PSCAD parsing, options, constitutive export, checkpoints and result
transport were exercised. Remote licensed PSCAD was not invoked. The full
documentation site build was not run; maintained examples, documentation/import
quality checks and the Gauntlet documentation-report consumer were exercised.

### Numerical, inference and allocation evidence

The same two-frequency air, buried and mixed two-cable fixtures were checked
against the independent pre-existing Z/Y snapshots at `rtol=1e-12`. All matched.
These small regression fixtures are not the earlier 101-frequency 18 kV case.

| Fixture | Original first compute | Final first compute | Original warm | Final warm |
|---|---:|---:|---:|---:|
| Air | 48.245 s | 46.815 s | 4.153 ms | 4.159 ms |
| Buried | 21.850 s | 19.692 s | 2.934 ms | 2.878 ms |
| Mixed | 20.457 s | 20.078 s | 3.989 ms | 3.963 ms |

First computations were more than 99.9% compilation time; these timings do not
promise a general speedup. The matching package-load check was 7.78 s versus
7.67 s initially, with essentially unchanged allocation count. Normal cached
loading/precompilation and absence of optional core imports also passed.

Real scalar and bundled `compute` calls infer concrete results, and composite/UQ
result/detail fields remain concrete. Trace/bundle/boundary tests retain the
intentional nested diagnostic bound.

Full-allocation profiling of the warmed air fixture identified an initial
12-allocation/384-byte increase entirely in setup/metadata, not frequency
evaluation. The existing metadata projection now traverses the typed selected
tuple directly, and the existing equation normalizer checks unique named-tuple
keys with a tuple filter rather than allocating a set-difference container.
Resolution counts, key order and rejection behavior are unchanged. Final check:
4,362 allocations / 358,000 bytes versus 4,421 / 363,296 before the migration.
Allocation count and bytes are both lower; no inlining hint or new helper remains.

### Existing failures and limitations

- Two mixed-earth CIM tests at 50 Hz and 10 kHz with
  `Γ=1e-4*(1+im)` exhaust the matrix error-budget policy. An untouched HEAD
  snapshot reproduced the exact estimates/targets and errors. Neither the
  equations, safeguards nor tolerances were changed to hide those failures.
- The optional plotting run had 786 passes, two marker failures and three SVG
  errors. The untouched snapshot reproduced both marker failures and the same
  SVG errors in the matrix tests. SVG export still calls
  `Base.require(LineCableModels, :CairoMakie)` although CairoMakie is not a direct
  dependency. That pre-existing backend-loading defect is outside this migration;
  no dependency was added to work around it.
- A pre-existing test compared the tuple-valued earth propagation constants
  with a vector. After reproducing it on the untouched source, the assertion was
  corrected to exact tuple-to-tuple equality, without changing numerical values.

Detailed logs and the initial/final comparison probes are retained under
`/tmp/lcm-records-*`; the workspace execution journal is
`local/nominal-record-ownership/2026-09-16/execution.md`.
The unrelated staged `dev/inspect_saved_gauntlet.jl` change was preserved.
No commit was made.
