# Computational engine

The physical equations and literature sources for every built-in formulation
are catalogued in [Transmission line parameters](transmission-line-parameters.md).
This developer guide describes how those identities are discovered, routed,
and extended without duplicating their scientific definitions.

LineCableModels separates a physical problem, a formulation, and the backend
that implements that formulation. A formulation identifies the equations or
mathematical method being requested and exists independently of a backend. A
backend participates by defining the dispatch required to compute it.

Formula identities are stable literature symbols. Registered literature IDs
use `:NameYYYY`, or at most `:NameNameYYYY` for a conventional joint name.
They are local to the receiving formula family: `:Papadopoulos2010` can select
the corresponding impedance law in `earth_impedance` and the corresponding
potential-coefficient law in `earth_admittance` without `Z`, `Y`, `Internal`,
or `External` suffixes. [`formula`](@ref) provides one uniform selection
wrapper; the receiving keyword supplies that namespace:

```julia
Formulation(
    internal_impedance=formula(:Schelkunoff1934),
    insulation_impedance=formula(:Ametani1980),
    earth_impedance=formula(:Papadopoulos2010),
    insulation_admittance=formula(:Marti2001),
    semicon_admittance=formula(:Ametani2004),
    earth_admittance=formula(:Papadopoulos2010),
    earth_properties=formula(:CIGRE2019),
    equivalent_earth=formula(:Xue2021; order=:after),
)
```

Bare symbols remain the concise default-only form. Formula-specific keyword
arguments are forwarded as route or assumption overrides. `order=:before` or
`:after` places an EHEM reduction relative to material frequency dependence;
ordinary formula slots reject an `order` value.

`EarthImpedance.Formula{:ID}` carries the identity in its concrete type. Its
`routes` tuple holds the leaf self, mutual, and propagation-constant formulas;
its `assumptions` tuple holds the physical approximations shared by those
routes. Calling a formula with one frequency's earth properties returns a
concrete, formula-owned functor that aggregates the shared numerical values:

```julia
formula = EarthImpedance.Formula(:Papadopoulos2010)
functor = formula(rho, epsilon, mu, jω, nothing)
value = functor(Val(:mutual), pair)
EarthImpedance.Γ(functor)
```

There are no author-named constructors. An experiment can replace one leaf
without spelling the whole recipe again:

```julia
formula = InternalImpedance.Formula(
    :Schelkunoff1934;
    inner=my_inner_surface_formula,
)
```

`routes(formula)` and `assumptions(formula)` make the selected recipe
inspectable. `propagation(formula)` reports whether an earth formula accepts an
explicit longitudinal Γ or fixes it to zero. `LineParametersProblem(...;
Γ=values)` passes explicit frequency-aligned values through the solve loop;
zero-assumption formulas reject a nonzero value.

The workspace resolves one `EarthPair` for every upper-triangular cable
interaction. A pair carries its physical heights, separation, and source and
target earth-layer indices. The frequency functor receives that pair at its
leaf call, so a route can distinguish overhead, underground, and mixed
interactions without rebuilding the full recipe. A single `:Pollaczek1926`
entry consequently owns its overhead, underground, and mixed impedance and
potential-coefficient leaves. The shared overhead integral is also selectable
directly as `:Carson1926`. `:Ametani2009` and `:Lucca1994` retain Carson and
Pollaczek for same-medium pairs and replace the mixed leaf with the selected
author's approximation. `:MartinsBritto2024` follows the same pair-complete
rule.
Support derivations do not become public formula IDs: the surface and
penetration-depth checks used by `:Xue2018` remain named leaves in its route
tuple, while the infinite-depth result is its default. An experiment can
replace any exposed leaf without creating another registry entry.

Each built-in formula lives in one `formulas/authoryear.jl` file. Formula
modules include those files in sorted order and require each file to return its
unique `Symbol` identifier. This is source discovery only: the hot loop has no
runtime registry or dictionary lookup. A contributed file defines its routes,
assumptions, formula call, and functor leaf dispatch together.

Insulation impedance uses the same discovery rule with one complete scalar
route per formula. Insulation and semicon admittance formulas are constitutive
relations applied to the full physical `Material` at each frequency; the
Coaxial Engine owns the common annular geometry and radial series aggregation.
`InsulationImpedance.formulas()` includes `:Ametani1980`;
`InsulationAdmittance.formulas()` includes `:Gustavsen2013` and `:Marti2001`;
`SemiconAdmittance.formulas()` includes `:Ametani2004`. They are selected
uniformly through `insulation_admittance` and `semicon_admittance`. A complete
experimental constitutive law can be supplied with `formula(:Marti2001;
route=my_route)` without changing the Coaxial matrix-assembly loop.

Frequency-dependent soil laws and equivalent homogeneous-earth models belong
to separate formula families under `EarthProps`. `EarthModel` continues to
store only static layer geometry and datasheet properties.

`EarthProps.FD` owns measured and material-physics relations with the scalar
contract `constitutive(formula, material, frequency)`. The relation selected by
`Formulation(earth_properties=:AuthorYear)` is applied only to soil; `nothing`
passes static properties through exactly and air is never modified.
`EarthProps.FD.formulas()` reports the discovered Longmire–Smith, Portela,
Alipio–Visacro, Datsios–Mikropoulos, Scott, Messier, Visacro–Portela,
Visacro–Alipio, and CIGRE WG C4.33 relations.

`EarthProps.EHEM` owns reductions required by homogeneous earth-return
formulations. `formula(:Xue2021; order=:after)` first evaluates every physical
layer and then reduces it; `order=:before` reduces the static layers and applies
FD to the artificial material afterward. These resolve to independent dispatch
paths. `formula(:Layer; layer=-1)` selects the explicit bottommost-layer policy,
which is the default and is not a literature formula.
`:MartinsBritto2020` reconstructs conductivity only, whereas `:Xue2021`
reconstructs conductivity and permittivity. Both registered routes implement
their published overhead-line scope and reject underground and mixed pairs.

At each frequency the Coaxial loop evaluates EHEM once per physical
`EarthPair`, maps the resulting material to a two-medium air/earth view, and
shares it between earth impedance and earth admittance. Physical layer indices
are never renumbered or written back to the `EarthModel`.

## Earth-free cable constants

`CableConstantsProblem`, `CableConstantsFormulation`, and `CableConstants`
belong to Engine. They reuse the registered internal-impedance,
insulation-impedance, insulation-admittance, and semicon-admittance formulas,
but not `LineParametersWorkspace` or either line-parameter assembly routine.
The default bundle is:

```julia
CableConstantsFormulation(
    internal_impedance = formula(:Schelkunoff1934),
    insulation_impedance = formula(:Ametani1980),
    insulation_admittance = formula(:Marti2001),
    semicon_admittance = formula(:Ametani2004),
)
```

`DataModel.flatten(design, frequency)` supplies the canonical components.
Contiguous components sharing one radial centre form one concentric assembly.
The Engine retains its innermost terminal, grounds all outward terminals,
assembles and reduces the local N-terminal series-impedance matrix, and
combines the physical dielectric layers between the core and first grounded
terminal in radial series. Earth impedance, earth admittance, EHEM, Γ,
position, transposition, and bundle reduction never enter this workflow.

`CableConstants(design; temperature=20, frequency=1e-3)` is the convenience
entry point. The result owns `cores`, aligned `R/L/C/G` vectors, and the
evaluation frequency. A conventional coaxial cable has one row and supports
`only(constants)`.

An external scalar relation can remain outside the built-in directory. It can
either supply a complete `EarthProps.FD.Formula(:Experiment, route, assumptions)`
or extend `constitutive(relation, material::EarthMaterial, frequency)`. The
resolved relation is concrete before the frequency scan. An external EHEM can
similarly use `EarthProps.EHEM.Formula(:Experiment, route, assumptions)` and
select its ordering with `AfterFD` or `BeforeFD`.

There is no `ReferenceEarthImpedance` category: whether a backend implements a
formulation does not change the formulation's place in the scientific
vocabulary. `:DeriSemlyen1981`, `:WedepohlWilcox1973`, `:Saad1996`,
`:Ametani2009`, and `:Lucca1994` describe formulae applicable to
homogeneous-earth models. All except the explicit PSCAD-only
`:DeriSemlyen1981` vocabulary are executable by the Coaxial backend; the PSCAD
backend independently maps the identifiers it supports to its own input
fields. Backend support is not a type-hierarchy category.

PSCAD's direct numerical integration setting is different. PSCAD exposes a
numerical integration choice through the same input field that selects an
earth formula and requires separate overhead and underground variants.
`PSCADBenchmarks.DirectNumericalIntegration` therefore remains local to that
backend and is not part of `Engine.EarthImpedance`.

## Completed-result read side

`AbstractCoreResult` marks direct LineCableModels-owned computation results;
`CableConstants` and `LineParameters` are the current core result types.
`AbstractResultSpace{T}` marks completed finite collections of stored core
results. Its element type remains open so an external solver's concrete result
can be stored without subtyping a LineCableModels type. Result-space
constructors reject abstract element types and nested result-space envelopes.

Core results own their scientific extraction methods.
[`observe`](@ref) reads native numerical values through function-object
selectors:

```julia
observe(parameters, Z)
observe(parameters, L, 1, 1, Colon())
observe(parameters, Y, angle, 1, 1, Colon())
```

The public `Z`, `L`, and other laconic accessors delegate to these methods.
Consumers do not inspect `LineParameters` storage or repeat the R/X/L/G/B/C
formulae.

Direct numerical access remains available:

```julia
parameters.Z[1, 1, :]
@view parameters.Y[1, 1, :]
Z(parameters, 1, 2)
@observe parameters L[1, 2, :]
```

`@observe` expands the indexed expression to `observe(parameters, L,
1, 2, Colon())`. With no source argument it constructs the same plain request
tuple without reading a result:

```julia
request = @observe R[:, :, :]
magnitude_request = @observe (Z, abs)[:, :, :]
magnitude = @observe parameters (Z, abs)[1, 2, :]
```

The request tuple is an implementation representation, not a second selector
type. `quantity(Z, abs)` and `quantity(Z, angle)` remain the transformed
scientific identities.

[`observables`](@ref) publishes only explicitly requested values:

```julia
published = observables(
    parameters,
    (
        (frequencies, Colon()),
        @observe(R[1, 1, :]),
    );
    units = (
        LineCableModels.Units.units(:base, :hertz),
        LineCableModels.Units.units(
            :base,
            :ohm;
            per = (:kilo, :meter),
        ),
    ),
)
```

Every positional payload contains only `values`, `quantity`, and `unit`.
Publication converts and detaches `values`; it does not attach labels, result
objects, execution options, Gridspace points, or Monte Carlo context.

Line plotting accepts explicit observable requests. Its public convenience
forms expand selectors such as `Z`, `real`, and `angle` once, at the plotting
surface, then enter the same request path. PlotBuilder groups completed
observations with the qualified `Units.family(::Quantity)` metadata. Series and
shunt identities return `Val(:series)` and `Val(:shunt)` respectively. Neither
PlotBuilder nor ReportBuilder owns another quantity or family map.

The qualified `Grammar.validate_observables` method is the single entitlement
check used by direct publication and generic reports. It validates the source
declaration, request identities, and positional unit alignment.
`Grammar.unit_targets` resolves a tuple of requests to aligned `UnitExpr`
values. A unit override may be `nothing`, a metric-prefix `Symbol`, an explicit
`UnitExpr`, or a quantity-keyed collection used by an entry-point normalizer.
Line plots, Monte Carlo plots, and reports all use this path.

`LineCableModels.Units` owns `Unit`, `UnitExpr`, `Quantity`, `units`,
`quantity`, `native_unit`, `display_unit`, `scale_factor`, `label`, and
`symbol`. PlotBuilder derives scientific axes from publication payloads.
ReportBuilder derives human-facing tables through `report`. Neither consumer
owns a quantity map, physical transform, or ordinary scientific unit string.

`Quantity{Q}` is a fieldless typed identity for extension methods and internal
publication payloads. Ordinary calls use scientific selector functions:

```julia
label(R)
symbol(Z, angle)
native_unit(R, :pul)
display_unit(Z, abs, :total)
```

The selector methods delegate through `quantity`; they do not contain a second
label or unit map. An external selector adds its identity and metadata at the
Units boundary:

```julia
function profile_response end

LineCableModels.Units.quantity(::typeof(profile_response)) =
    LineCableModels.Units.Quantity{:profile_response}()

LineCableModels.Units.native_unit(
    ::LineCableModels.Units.Quantity{:profile_response},
) = LineCableModels.Units.units(:base, :ohm)

LineCableModels.Units.display_unit(
    ::LineCableModels.Units.Quantity{:profile_response},
) = LineCableModels.Units.units(:milli, :ohm)

LineCableModels.Units.label(
    ::LineCableModels.Units.Quantity{:profile_response},
) = "Profile response"

LineCableModels.Units.symbol(
    ::LineCableModels.Units.Quantity{:profile_response},
) = "u"

label(profile_response)
display_unit(profile_response)
```

`Quantity`, `Unit`, and `UnitExpr` remain qualified extension vocabulary. The
package root exports only the six metadata functions used with scientific
selectors.

Higher-order results remain containers of owned products. `result`,
`statistics`, `samples`, and `histograms` select those products; they are not
zero-argument aliases for publication. UQ reads trial results with
`observe`, while its retained statistics, samples, and histograms implement
the same selector grammar for later publication. Monte Carlo run settings and
resolved point values are read through `root_seed`, `point_seed`,
`trial_count`, `confidence`, `cdf_tolerance`, and `sampling_distribution`.
Consumers do not inspect the result or its formulation fields.

## Coaxial workspace and supplemental output

`LineCableModelsCoaxial` solves concentric coaxial assemblies. A sector,
stranded, or otherwise nonconcentric part must be represented by the equivalent
round/concentric properties owned by DataModel before it reaches this backend.
The backend then owns frequency scans, self and mutual line parameters, earth
effects, and reduction; it does not redefine cable-design equivalence or modal
coordinates.

`LineParametersWorkspace` is the coaxial backend's per-computation working
state. Its constructor adapts a completed physical system once, evaluates
earth data, constructs cable and reduction indices, and allocates the matrices
and adaptive-quadrature segment storage used by the frequency loop. The
workspace is not a result, public data model, or alternate cable
representation. Its shunt solvers consume DataModel's ordered physical
dielectric layers directly, so the analysis is independent of the frequency
used when a lossy homogeneous export representation is requested.

The workspace separates four owned concerns:

- `input`: immutable numerical input derived from the problem;
- `invariants`: reusable physical values and index maps;
- `buffers`: mutable storage reused while solving every frequency;
- `capture`: optional diagnostic matrices allocated before the loop.

The ordinary result is always `LineParameters`. Requesting
`options=(trace=true,)` retains completed diagnostic arrays under
`details(parameters).trace`; it does not select another result type.

## Modal transformations

Modal decomposition is independent of the backend that produced fully coupled
phase-domain matrices. `LineCableModels.Transforms` owns its own problem,
formulation, registered formula files, and default backend:

```julia
phase = compute(line_problem, line_formulation)
modal = compute(
    ModalTransformationProblem(phase),
    ModalTransformationFormulation(
        formula(:Chrysochos2014; convergence=1e-8),
    ),
)
rebuilt = compute(ModalTransformationProblem(modal))
```

`LineCableModelsModal` is the default backend for this workflow. Each formula
has one route and returns a `ModalOperators` value containing the complete
frequency-dependent phase-to-modal voltage and current tensors. The shared
backend applies those operators to both `Z` and `Y`.

The modal `LineParameters` result carries `ModalDomain(operators, formula)` as
its concrete domain value. The stored operators preserve the resolved mode
order, scaling, and complex phase convention and make the transformation
bidirectional without rerunning the decomposition. An operator-less modal
result cannot be constructed through the admitted `LineParameters` interface.

Built-in formula files are included deterministically and selected by symbols:

```julia
ModalTransformationFormulation(formula(:Fortescue))
ModalTransformationFormulation(formula(:Chrysochos2014; convergence=1e-8))
ModalTransformationFormulation(formula(:Fan2009; history_weight=0.3))
ModalTransformationFormulation(formula(:Wedepohl1996; convergence=1e-9))
```

`Chrysochos2014` solves independently tracked eigenpairs with a real-valued
Levenberg–Marquardt step. `Wedepohl1996` follows the corresponding complex
Newton–Raphson route. `Fan2009` postprocesses conventional eigensolutions with
optimal assignment, complex phase alignment, and Procrustes alignment of
coalescent eigenspaces.

There are no author-named formulation structs or runtime registry lookups.
The built-in route can be replaced at the call site without changing its
stable identifier:

```julia
experiment = ModalTransformationFormulation(
    :Fortescue;
    route=my_route,
)
```

A wholly new formula uses the same typed path with an explicit assumptions
tuple:

```julia
formula = Transforms.Formula(
    :Experiment,
    my_route,
    (tolerance=1e-8,),
)
experiment = ModalTransformationFormulation(formula)
```

For a contributed built-in, one `formulas/authoryear.jl` file defines
`assumptions(::Val{:ID})`, `description(::Formula{:ID})`, and the single
`(::Functor{:ID})(parameters, assumptions)` route, then returns `:ID` from the
file. Discovery happens at module load; numerical dispatch remains static.

[`ComputationDetails`](@ref) is an alias for `NamedTuple`.
[`computation_details`](@ref) reads the fixed-key details tuple owned by a
registered computation owner. There is no general method: an unregistered
owner raises `MethodError`.

Higher-order calculations normalize the inner formulation once through the
qualified `Grammar.computation_owner` method before collecting retained
records. This produces the unparameterized type token required by the existing
`Val(OwnerType)` details contract; it is not a registry.

[`ParametricResult`](@ref), [`LinearErrorResult`](@ref), and
[`MonteCarloResult`](@ref) store the concrete details tuple type. Retention is
disabled by default, so `details(result) === (;)`. The higher-order formulation
owns the retention option:

```julia
Combinatorial(formulation; options=(retain_details=true,))
LinearError(formulation; options=(retain_details=true,))
MonteCarlo(formulation; trials=100, options=(retain_details=true,))
```

Parametric and linear calculations retain `(points=records,)`, with one record
per core result. Monte Carlo retains `trials`, `failures`, and
`failure_summary`, each aligned by Gridspace point. `trials` contains one inner
computation record per accepted trial. Each failure record contains the
attempt, target trial, failure stage, realised argument tuple, error type and
message, and a bounded stack summary. Statistics, samples, histograms, seeds,
and accepted-trial counts remain dedicated result fields.

## Formulation options

[`formulation_options`](@ref) validates values that alter the mathematical
calculation represented by a formulation. Dispatch uses the formulation owner
type rather than the public construction selector:

```julia
formulation_options(Val(LineParametersFormulation), options)
```

The default line-parameter formulation owns:

- bundle and Kron reduction.
- ideal transposition.
- temperature correction.

The normalised named tuple is stored in `LineParametersFormulation.options`.
`PSCADFormulation` currently has no formulation options because its method
bundle already contains every mathematical choice it owns.

`Formulation()` constructs the default method bundle without a backend
tag. `LineParametersFormulation` owns the formulation options;
`LineCableModelsCoaxial` separately owns execution. Symbol and `Val` selectors
remain available for external backends, but there is no
`:line_cable_models` or legacy `:analytical` selector.

## Computation options

[`computation_options`](@ref) validates values belonging to one execution.
Computation options do not change the selected mathematical formulation and are not
stored in it.

The coaxial backend accepts:

```julia
(
    verbosity = (default = 0,),
    output_basis = :pul,
    trace = false,
)
```

`trace=true` preallocates diagnostic capture with the workspace and attaches
the retained matrices to `details(result).trace` after computation.

The PSCAD backend accepts:

```julia
(
    output_stem = "case_name",
    remote = remote_config,
    verbosity = (default = 0, PSCAD = 2),
    output_basis = :pul,
)
```

`remote` must be a `PSCADBenchmarks.RemoteConfig`. `output_stem` names files
created by that execution. Neither value belongs to `PSCADFormulation`.

Both option sets are ordinary `NamedTuple`s, aliased as
[`FormulationOptions`](@ref) and [`ComputationOptions`](@ref). Callers can
compose them with `merge`. Each owner rejects unknown keys and returns a
fixed-key normalised tuple. There is no general fallback and no conversion
from dictionaries, pairs, or `nothing`.

`MonteCarlo` owns a separate outer computation-option tuple. Its normalised
keys are `retain_details`, `on_error`, and `max_failures`. `on_error=:fail` is
the default and rethrows every exception. `on_error=:retry` requires
`retain_details=true` and rejects only `DomainError` realisations until the
requested accepted-trial count is reached or `max_failures` is exhausted.
Other exception types always propagate immediately.

## Gauntlet routing

`GauntletCase` coordinates two computations but does not own either backend's
keys. Its computation options form an outer tuple:

```julia
(
    output_basis = :pul,
    reference = (
        output_stem = "case_name",
        remote = remote_config,
        verbosity = (default = 0, PSCAD = 2),
    ),
    candidate = (
        verbosity = (default = 0,),
    ),
    benchmark = (
        samples = 10,
        seconds = 10.0,
    ),
)
```

The runner validates only this outer shape. The runner passes `reference` and
`candidate` to the corresponding `computation_options` methods and forces the same
`output_basis` into both so their results are comparable. Live and record runs
load a configured remote endpoint only when `reference.remote` is absent.
Snapshot runs never load remote configuration. Run mode, snapshot writing,
comparison tolerances, expected dimensions, and port ordering are not
computation options.

## Extending the engine

An external package may own a backend identity and a separate formulation
type. The backend's `compute` method normalises execution options before doing
work:

```julia
import LineCableModels:
    AbstractFormulation,
    AbstractCoreResult,
    ComputationOptions,
    FormulationOptions,
    computation_owner,
    computation_options,
    compute,
    formulation_options

struct ExternalEngine end

struct ExternalFormulation{O <: NamedTuple} <: AbstractFormulation
    options::O
end

function formulation_options(
    ::Val{ExternalFormulation},
    options::NamedTuple,
)::FormulationOptions
    isempty(options) || throw(ArgumentError("unsupported formulation option"))
    return (;)
end

function computation_options(
    ::Val{ExternalEngine},
    options::NamedTuple,
)::ComputationOptions
    unknown = filter(key -> key != :tolerance, keys(options))
    isempty(unknown) || throw(ArgumentError("unsupported computation option"))
    normalized = merge((tolerance = 1.0e-8,), options)
    normalized.tolerance > 0 || throw(ArgumentError("tolerance must be positive"))
    return (tolerance = Float64(normalized.tolerance),)
end

computation_owner(::ExternalFormulation) = ExternalEngine

function compute(
    ::ExternalEngine,
    problem,
    formulation::ExternalFormulation;
    options::NamedTuple = (;),
)
    execution = computation_options(Val(ExternalEngine), options)
    # Use `problem`, `formulation`, and `execution` here.
end
```

An external implementation does not need a dedicated options struct or private
wrapper around the two normalisation functions. If it omits either Grammar
method, Julia raises `MethodError`.

The same backend may expose supplemental output without changing the generic
higher-order result types:

```julia
import LineCableModels: ComputationDetails, computation_details

struct ExternalResult <: AbstractCoreResult
    parameters
    diagnostics::NamedTuple
    raw::Dict{String,Any}
end

function computation_details(
    ::Val{ExternalEngine},
    output::ExternalResult,
)::ComputationDetails
    return (
        diagnostics=output.diagnostics,
        raw=output.raw,
    )
end
```

An external backend that cannot modify or wrap its solver's concrete return
type may store that type directly in a result space. `AbstractCoreResult` marks
owned direct results; it is not an admission requirement for external result
payloads.

The outer keys and their types are fixed for `ExternalFormulation`. Dynamic
vendor channels remain inside the explicit `raw` leaf. ParametricBuilder and
UQ collect these records only when `retain_details=true`; they do not inspect
the fields.

## Reports and XLSX output

[`report`](@ref) executes `entitle`, `select`, `tabulate`, `illustrate`,
`encode`, `write`, and `finish`. In-memory definitions implement explicit
no-op `encode` and `write` methods and return `ReportArtifact.output ===
nothing`.

[`XLSXReportDefinition`](@ref) owns the human-facing line-parameter workbook:

```julia
using XLSX

artifact = report(
    XLSXReportDefinition(file_name="line_parameters.xlsx"),
    parameters,
)
artifact.output
```

ReportBuilder selects values through `observables`, builds one wide table with
coordinate columns followed by one column per observed quantity, and encodes a complete
[`LineCableModels.ReportBuilder.XLSXWorkbook`](@ref) containing the destination,
ordered sheet names, and final cell strings. Loading XLSX activates the package
extension that writes only this encoded description and records its path in
[`ReportArtifact`](@ref). Relative and default paths resolve from the caller's
current working directory; the package source tree is never the implicit
destination. `export_data(:xlsx, parameters; ...)` remains a thin ImportExport
convenience call that returns the same path. ImportExport owns no second
workbook implementation.
