# Computational engine

LineCableModels separates a physical problem, a formulation, and the backend
that implements that formulation. A formulation identifies the equations or
mathematical method being requested. A formulation exists independently of any backend.
A backend participates by defining the dispatch required to compute that
formulation.

For example, `EarthImpedance.Saad()` is an Engine-owned earth-impedance
formulation. The built-in analytical backend currently defines no evaluation
method for it, so analytical computation fails with an ordinary `MethodError`.
The PSCAD backend defines the corresponding PSCAD input mapping and can execute
the same formulation. There is no `ReferenceEarthImpedance` category: whether
a backend implements a formulation does not change the formulation's place in
the scientific vocabulary. `Deri`, `Wedepohl`, `Saad`, `Ametani`, and `Lucca`
describe formulae applicable to homogeneous-earth models and are direct
children of `EarthImpedanceFormulation`. Backend support is not a type-hierarchy
category.

PSCAD's direct numerical integration setting is different. PSCAD exposes a
numerical integration choice through the same input field that selects an
earth formula and requires separate overhead and underground variants.
`PSCADBenchmarks.DirectNumericalIntegration` therefore remains local to that
backend and is not part of `Engine.EarthImpedance`.

## Completed-result read side

Primitive calculation results own their scientific extraction methods.
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

[`observables`](@ref) publishes only explicitly requested values:

```julia
published = observables(
    parameters,
    (
        frequency = (frequencies, Colon()),
        resistance = (R, 1, 1, Colon()),
    );
    units = (
        resistance = LineCableModels.Units.units(
            :base,
            :ohm;
            per = (:kilo, :meter),
        ),
    ),
)
```

Every published field contains only `values`, `quantity`, and `unit`.
Publication converts and detaches `values`; it does not attach labels, result
objects, execution options, Gridspace points, or Monte Carlo context.

`LineCableModels.Units` owns `Unit`, `UnitExpr`, `QuantityTag`, `units`,
`quantity`, `native_unit`, `display_unit`, `scale_factor`, `label`, and
`symbol`. PlotBuilder derives scientific axes from publication payloads.
ReportBuilder derives human-facing tables through `report`. Neither consumer
owns a quantity map, physical transform, or ordinary scientific unit string.

Higher-order results remain containers of owned products. `result`,
`statistics`, `samples`, and `histograms` select those products; they are not
zero-argument aliases for publication. UQ reads primitive trials with
`observe`, while its retained statistics, samples, and histograms implement
the same selector grammar for later publication.

## Formulation options

[`formulation_options`](@ref) validates values that alter the mathematical
calculation represented by a formulation. Dispatch uses the formulation owner
type rather than the public construction selector:

```julia
formulation_options(Val(AnalyticalFormulation), options)
```

The analytical formulation owns:

- bundle and Kron reduction.
- ideal transposition.
- temperature correction.
- parameter or trace output selection.

The normalised named tuple is stored in `AnalyticalFormulation.options`.
`PSCADFormulation` currently has no formulation options because its method
bundle already contains every mathematical choice it owns.

`Val(:analytical)` remains the public selector accepted by `Formulation`.
`Val(AnalyticalFormulation)` is the internal owner token used by the option
grammar. The two selectors are not interchangeable.

## Computation options

[`computation_options`](@ref) validates values belonging to one execution.
Computation options do not change the selected mathematical formulation and are not
stored in it.

The analytical backend accepts:

```julia
(
    verbosity = (default = 0,),
    output_basis = :per_length,
)
```

The PSCAD backend accepts:

```julia
(
    output_stem = "case_name",
    remote = remote_config,
    verbosity = (default = 0, PSCAD = 2),
    output_basis = :per_length,
)
```

`remote` must be a `PSCADBenchmarks.RemoteConfig`. `output_stem` names files
created by that execution. Neither value belongs to `PSCADFormulation`.

Both option sets are ordinary `NamedTuple`s, aliased as
[`FormulationOptions`](@ref) and [`ComputationOptions`](@ref). Callers can
compose them with `merge`. Each owner rejects unknown keys and returns a
fixed-key normalised tuple. There is no general fallback and no conversion
from dictionaries, pairs, or `nothing`.

## Gauntlet routing

`GauntletCase` coordinates two computations but does not own either backend's
keys. Its computation options form an outer tuple:

```julia
(
    output_basis = :per_length,
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

An external backend owns a formulation type and both option methods that apply
to it. The backend's `compute` method normalises execution options before doing work:

```julia
import LineCableModels:
    AbstractFormulation,
    ComputationOptions,
    FormulationOptions,
    computation_options,
    compute,
    formulation_options

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
    ::Val{ExternalFormulation},
    options::NamedTuple,
)::ComputationOptions
    unknown = filter(key -> key != :tolerance, keys(options))
    isempty(unknown) || throw(ArgumentError("unsupported computation option"))
    normalized = merge((tolerance = 1.0e-8,), options)
    normalized.tolerance > 0 || throw(ArgumentError("tolerance must be positive"))
    return (tolerance = Float64(normalized.tolerance),)
end

function compute(problem, formulation::ExternalFormulation; options::NamedTuple = (;))
    execution = computation_options(Val(ExternalFormulation), options)
    # Use `problem`, `formulation`, and `execution` here.
end
```

An external implementation does not need a dedicated options struct or private
wrapper around the two normalisation functions. If it omits either Grammar
method, Julia raises `MethodError`.
