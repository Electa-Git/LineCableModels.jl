# Extension API

Extension code adds methods to the owner that defines their meaning. The public
calculation entry points are `compute`, `observe`, `observables`,
`report`, `plot`, and `preview`.

`PlotBuilder` owns only optional plotting entry points and the live `UIPlot`
handle. Scientific owners expose observations, quantities, geometry, and
physical property ranges. The Makie extension owns request normalization,
palettes, legend grouping, layout helpers, widgets, and native rendering.

## Formula selection and descriptions

`formula_id` identifies a scientific selection. `description` supplies its
human-readable text. Registered formula types implement
`description(::Type{<:OwnedFormula{:ID}}; compact=false)`, and their instances
delegate with the same keyword. The default is the detailed scientific
description. `compact=true` is the short name used in legends. Override this
method to customize a name, not `formula_id`.

A family declares independently selectable children through
`pairs(Family.Formula; quantity=nothing)`, returning ordered `slot => family`
pairs, or an empty mapping for a scalar leaf. Constructors, retained-record
readers and formulation projections use that same owner declaration.
Internal impedance admits `inner`, `outer`, `transfer`. The external earth
families admit `air`, `earth`, `mixed`. Explicit composites retain every branch,
including default branches, in quantity-relevant legends.

Formulation owners expose ordered `(owner, route_tuple) => selection` pairs
through `pairs(source; quantity)` and `pairs(owner, retained; quantity)`.
Each selected leaf is paired with its declaration controls.
`formulation_options(selection)` reads the live selection's `FormulationOptions`.
Consumers obtain identity from `formula_id`, structure and ordering from `pairs`,
and display text from `description`. The contextual
`description(owner, selection; compact)` method describes the backend's equations,
including PSCAD's native defaults.

## User-owned equations

Select a concrete type directly in the physical slot, such as
`Formulation(earth_properties=MySoil(...))`.
`FormulaMethod(selected, operation, Val(...), ...)` calls the operation with the
selected object first. IDs are for inspection only. `:default` resolves to an
explicit implementation before physical validation or computation.

| Family | User type and owning operation |
|---|---|
| Internal impedance | `Engine.InternalImpedanceFormulation`. `InternalImpedance.internal_impedance(selected, Val(kind), functor, workspace)` |
| Insulation impedance | `Engine.InsulationImpedanceFormulation`. `InsulationImpedance.insulation_impedance(selected, r_in, r_ex, mu_r, s, parameters, options, workspace)` |
| Insulation / semicon admittivity | Corresponding `Engine.*AdmittanceFormulation`. `insulation_material` / `semicon_material(selected, material, frequency, temperature, parameters, options, workspace)` |
| Earth impedance / potential | Corresponding `Engine.Earth*Formulation`. `earth_impedance` / `earth_potential_coefficient(selected, Val(kind), Val(source), Val(target), functor, pair, workspace)` |
| Soil frequency dependence | `Earth.FrequencyDependent.FrequencyDependentFormulation`. `earth_material(selected, material, frequency, parameters, options, workspace)` |
| Temperature dependence | `Materials.TemperatureDependent.TemperatureDependentFormulation`. `temperature_resistivity(selected, material, temperature, parameters, options, workspace)` |
| Equivalent earth | `Earth.EquivalentHomogeneous.AbstractRule`. `equivalent_material(selected, Val(kind), Val(source), Val(target), rho, eps_r, mu_r, model, pair, frequency, parameters, options, workspace)` |
| Modal decomposition | `AbstractFormulation`, selected by `ModalAnalysisFormulation`. `Engine.initialize_buffers(selected, T, input, invariants, common)` and `ModalAnalysis.decompose!(selected, workspace, parameters, options)` |
| Local shunt geometry | `Engine.ShuntModelFormulation`. `Engine.internal_shunt_response(selected, design, geometry, T, material_selections, solutions, design_index)` during blueprint construction |
| Pipe applicability | `Engine.PipeImpedanceFormulation`. `Formulation(backend, selected, Val(topology))`. No analytical pipe equation is supplied. |

`parameters` are model data and `options` are formulation-owned physical choices
and numerical controls. Custom constructors validate and normalize their own `FormulationOptions`.
Execution controls use `ComputationOptions`. Completed supplemental output uses
`ComputationDetails`. Read their payloads explicitly through `.data`.
Indexed families declare numerical defaults for the actual selected type and
case with `formulation_options(::FormulaMethod{<:MyType,typeof(operation),...})`.
Earth field equations consume evaluated material properties and define their
own wave numbers and field approximations. Material-law families implement
`constitutive` with a valid material argument.
Unified's formulation options accept a prescribed longitudinal coefficient Γ
as a scalar or a frequency-aligned vector. Its precision participates in
allocation of the calculation's numerical storage. The positive-time convention
and explicit convention conversions are documented with the Unified equations.

Coaxial equations and material laws receive the defining computation workspace,
or `nothing` for a standalone evaluation that does not require it. Numerical arrays are in
`workspace.buffers`. The workspace input and bindings remain the authority for
geometry and material mappings. An allocation method is unnecessary for an algebraic equation.
For an equation that needs scratch, extend
`Engine.initialize_buffers(selected, T, input, invariants, buffers)` to return
the record extended with owned arrays only. Only selections reached by the
required indexed calls participate in initialization, before evaluating materials
or equations. Unused recipe branches remain unallocated.
The default uses only existing storage. Existing arrays may not be replaced.
Numerical formulas provision the common quadrature storage through
`Engine.initialize_buffers(Val(:quad), T, input, invariants, buffers)`. An
integration option is not a capability declaration. Cable constants uses this
same buffer-initialization method with its local numerical input and no earth invariants.

Each formula defines its complete integrand, transformations, Jacobians, branch
choices and physical subdivision hints. `SpectralIntegral` contains only that
callable. `integrate` passes it and the supplied numeric subdivision points to
QuadGK. Physical expressions and subdivision choices belong to each formula.

`ShuntModel` owns conductor and dielectric geometry extraction, numerical coefficients and the requested fallback
model. Its `blueprint_dependencies` methods identify the actual local selections
that affect those coefficients. Engine uses that dependency record for reuse
within one blueprint construction. Completed coefficient blocks are retained
after construction. The charge-system factorization and workspace are reused
only during construction.

The existing `Engine.earth_bindings` constructor binds material interactions and
output entries. A coupled formula can extend its selected-type method to require
the complete system and its two-selection method to bind compatible Z/P consumers.
`Engine.earth!` then performs the actual calculation from completed material
inputs and writes the selected destinations. Unpaired selections use its ordinary
indexed implementation. Each selection uses the Engine frequency sequence and
its computation workspace.

Internal state is constructed once per conductor and frequency through the selected
type's physical constructor and `InternalImpedance.Functor`. Shunt geometry
construction runs once before the frequency loop and returns frequency-independent
blueprint blocks and diagnostics. A consuming earth
equation explicitly admits a custom equivalent-earth rule using `validate`.

Results use result-type and unit checks for each formula family, including custom types.
Impedances are in Ω/m, material admittivities in S/m, earth potential coefficients
in m/F, and temperature-law resistivities in Ω·m. Scalar material laws return
`EarthMaterial` or a finite scalar as appropriate. Modal equations write
`workspace.Tv`, `workspace.Ti`, and `workspace.roots`. The owner constructs
`ModalOperators` and intrinsic coefficients after decomposition. Supply
`formula_id`, `description`, `NamedTuple` and
`formulation_options` methods for metadata. Serialized identities and data do
not reconstruct executable methods. See the [temperature-law example](engine.md#Cable-material-temperature-dependence).

### A Val-dispatched modal equation

The following complete one-mode example implements a custom modal equation.
The modal workspace supplies common slices, coordinate conversion scratch, the
admittance-impedance product, eigenpair history and voltage-vector scratch.
`initialize_buffers` extends that record only for additional numerical work.
It assumes a completed one-mode phase scan named `phase`, with nonzero
diagonal coefficients and known source length. The example is algebraic. It
does not replace a broadband modal model.

```julia
using LineCableModels
import LineCableModels.Engine: initialize_buffers, description
import LineCableModels.Grammar: formulation_options, FormulationOptions
import LineCableModels.ModalAnalysis: decompose!, Formula
import LineCableModels: FormulaMethod

description(::Type{<:Formula{:diagonal_example}}; compact=false) =
    compact ? "diagonal example" : "one-mode diagonal example"
formulation_options(::FormulaMethod{<:Formula{:diagonal_example},typeof(decompose!)}) =
    FormulationOptions()

function initialize_buffers(::Val{:diagonal_example}, ::Type{T}, input,
        invariants, common) where {T<:Complex}
    invariants.n == 1 || throw(DimensionMismatch("diagonal example requires one mode"))
    return merge(common, (diagonal_product=Vector{T}(undef,invariants.nf),))
end

function decompose!(::Val{:diagonal_example}, workspace,
        parameters::NamedTuple, options::FormulationOptions)
    scratch=workspace.buffers.diagonal_product
    for k in eachindex(scratch)
        z=workspace.input.Z[1,1,k]/workspace.input.root_scale
        y=workspace.input.Y[1,1,k]/workspace.input.root_scale
        scratch[k]=z*y
        root=sqrt(scratch[k])
        (real(root)<0 || (iszero(real(root)) && imag(root)<0)) && (root=-root)
        workspace.roots[1,k]=root*workspace.input.root_scale
        workspace.Tv[1,1,k]=one(root)
        workspace.Ti[1,1,k]=one(root)
        workspace.diagnostics.eigen_residual[1,k]=zero(real(root))
        workspace.diagnostics.iterations[1,k]=0
        workspace.diagnostics.converged[1,k]=true
    end
    return workspace
end

selected=ModalAnalysisFormulation(Formula(Val(:diagonal_example)))
modal=compute(ModalAnalysisProblem(phase),selected)
segment=PropagationParameters(modal)
size(gamma(modal)) == (1,length(frequencies(phase)))
size(H(segment)) == size(gamma(modal))
```

Only the selected equation's `initialize_buffers` method runs. The common
workspace supplies per-frequency normalization and coordinate scratch. The
equation owns `diagonal_product`. It writes phase-row by mode-column bases
and mode-by-frequency roots. The owner checks structural shape and finite arithmetic before computing intrinsic coefficients. It then copies the returned arrays and records diagnostics. Numerical targets are reported as warnings and facts,
without becoming result-admission rules.

## Input validation

A materialized input defines one direct `validate(::OwnedType)` method. The method
returns its argument unchanged or throws a native Julia exception identifying
the rejected field and value:

```julia
import LineCableModels: validate

struct Annulus{T <: Real}
    r_in::T
    r_ex::T

    function Annulus{T}(r_in::T, r_ex::T) where {T <: Real}
        return validate(new{T}(r_in, r_ex))
    end
end

function validate(value::Annulus)
    value.r_in >= zero(value.r_in) || throw(DomainError(
        value.r_in,
        "Annulus.r_in must be nonnegative"
    ))
    value.r_ex > value.r_in || throw(DomainError(
        value.r_ex,
        "Annulus.r_ex must be greater than r_in"
    ))
    return value
end
```

Constructors normalize their admitted grammar. The owning `validate` method
checks the completed value directly and returns it unchanged.

```@docs
LineCableModels.InputValidation
```

## Grammar, observations, and units

```@autodocs
Modules = [
    LineCableModels.Grammar,
    LineCableModels.Units,
]
Order = [:module, :constant, :type, :function, :macro]
Filter = developer_reference_entry
Public = true
Private = false
```

## Data model and calculations

```@autodocs
Modules = [
    LineCableModels.DataModel,
    LineCableModels.Engine,
    LineCableModels.ParametricBuilder,
    LineCableModels.ImportExport,
]
Order = [:module, :constant, :type, :function, :macro]
Filter = developer_reference_entry
Public = true
Private = false
```

## Report definitions

```@autodocs
Modules = [LineCableModels.ReportBuilder]
Order = [:module, :constant, :type, :function, :macro]
Filter = developer_reference_entry
Public = true
Private = false
```

## Plotting functions

```@docs
LineCableModels.PlotBuilder
```

## Index

```@index
Pages = ["extensions.md"]
Order = [:constant, :type, :function, :macro]
```
