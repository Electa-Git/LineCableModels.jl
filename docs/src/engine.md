# Computational engine

LineCableModels separates a physical problem, a formulation, and the backend
that implements that formulation. A formulation identifies the equations or
mathematical method being requested. It exists independently of any backend.
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
children of `EarthImpedanceFormulation`; backend support is not a type-hierarchy
category.

PSCAD's direct numerical integration setting is different. PSCAD exposes a
numerical integration choice through the same input field that selects an
earth formula and requires separate overhead and underground variants.
`PSCADBenchmarks.DirectNumericalIntegration` therefore remains local to that
backend and is not part of `Engine.EarthImpedance`.

## Formulation options

[`formulation_options`](@ref) validates values that alter the mathematical
calculation represented by a formulation. Dispatch uses the formulation owner
type rather than the public construction selector:

```julia
formulation_options(Val(AnalyticalFormulation), options)
```

The analytical formulation owns:

- bundle and Kron reduction;
- ideal transposition;
- temperature correction;
- parameter or trace output selection.

The normalized named tuple is stored in `AnalyticalFormulation.options`.
`PSCADFormulation` currently has no formulation options because its method
bundle already contains every mathematical choice it owns.

`Val(:analytical)` remains the public selector accepted by `Formulation`.
`Val(AnalyticalFormulation)` is the internal owner token used by the option
grammar. They are deliberately not interchangeable.

## Computation options

[`computation_options`](@ref) validates values belonging to one execution.
These values do not change the selected mathematical formulation and are not
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
created by that execution; neither value belongs to `PSCADFormulation`.

Both option roles are ordinary `NamedTuple`s, aliased as
[`FormulationOptions`](@ref) and [`ComputationOptions`](@ref). Callers can
compose them with `merge`. Each owner rejects unknown keys and returns a
fixed-key normalized tuple. There is no general fallback and no conversion
from dictionaries, pairs, or `nothing`.

## Gauntlet routing

`GauntletCase` coordinates two computations but does not own either backend's
keys. Its computation options form an envelope:

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

The runner validates only this outer shape. It passes `reference` and
`candidate` to their respective backend normalizers and forces the same
`output_basis` into both so their results are comparable. Live and record runs
load a configured remote endpoint only when `reference.remote` is absent.
Snapshot runs never load remote configuration. Run mode, persistence authority,
comparison tolerances, expected dimensions, and port ordering are not
computation options.

## Extending the engine

An external backend owns a formulation type and both option methods that apply
to it. Its `compute` method normalizes execution options before doing work:

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

No `ExternalOptions` composite type or backend-named private normalization
function is needed. If the owner does not implement one of the Grammar
functions, that omission remains visible through ordinary Julia dispatch.
