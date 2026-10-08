# Extension API

Extension code adds methods to the owner that defines their meaning. The public
computation entry points are `compute`, `observe`, `observables`,
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

## User-owned formulas

A formula is a collection of expressions. A user-owned formula is a concrete type below its
family's formulation supertype. It stores its model `parameters` and its normalized
`options`, a `FormulationOptions` record, and the user selects it directly in its physical
slot, such as `Formulation(earth_properties=MySoil(...))`. Its methods dispatch on its own
type: the family operation for each of its expressions, `formulation_options` for the
options of each expression, and the metadata methods `formula_id`, `description` and
`NamedTuple`. A formula that checks its input or shares values adds a `Functor` method, and
one that needs arrays adds an `initialize_buffers` method. IDs are for inspection only.
`:default` resolves to an explicit implementation before physical validation or computation.

### Expressions

`Commons.Expression(formula, operation, selectors...)` names one expression of a formula: the
formula object, the family operation and the `Val` selectors of the part, such as the kind
and the source and target layers of an earth interaction, or the inner, outer or transfer
kind of an internal impedance. Calling the expression with a Functor and a workspace evaluates
`operation(formula, selectors..., functor, workspace)`. The formula object comes first, and
the selectors follow it.

A formula declares the defaults of each expression with
`formulation_options(::Expression{<:MyFormula, typeof(operation), ...})`. The family projects
the options supplied to the formula onto the expressions that consume them. Before the
frequency loop, `validate(expression)` checks that the operation has a method for the
formula and its selectors.

### The Functor

`Commons.Functor(formula, input, state)` is a formula at one evaluation point. `input`
stores the values of that point, such as the material and the frequency, with the options
of the expression evaluated there. `state` stores plain values that the expressions of the
formula share at that point, such as Schelkunoff's scaled Bessel values or Unified's
per-frequency system. Arrays come from `workspace.buffers`, never from the state. Every
family evaluates its expressions through the same path:

```julia
functor = Functor(formula, input; workspace)
value = Expression(formula, operation, selectors...)(functor, workspace)
# This evaluates operation(formula, selectors..., functor, workspace).
```

Without a method of its own, a formula gets an empty state. A family or a formula that
checks its input, shares values or reads buffers adds a method
`Functor(formula::MyFormula, input::NamedTuple; workspace)` that returns
`Functor(formula, input, state)`. `Functor(functor, extension)` keeps the formula and the
state and extends the input, such as with one conductor pair of an earth calculation.

| Family | Formulation supertype | Operation | Functor input |
|---|---|---|---|
| Internal impedance | `Engine.InternalImpedanceFormulation` | `InternalImpedance.internal_impedance(formula, Val(kind), functor, workspace)`, with `kind` one of `:inner`, `:outer`, `:transfer` | `r_in`, `r_ex`, `rho`, `mu_r`, `jω` and the `options` of that kind |
| Insulation impedance | `Engine.InsulationImpedanceFormulation` | `InsulationImpedance.insulation_impedance(formula, functor, workspace)` | `r_in`, `r_ex`, `mu_r`, `jω`, `options` |
| Insulation and semicon admittivity | `Engine.InsulationAdmittanceFormulation`, `Engine.SemiconAdmittanceFormulation` | `insulation_material` or `semicon_material(formula, functor, workspace)` | `material`, `frequency`, `temperature`, `options` |
| Earth impedance and potential | `Engine.EarthImpedanceFormulation`, `Engine.EarthAdmittanceFormulation` | `earth_impedance` or `earth_potential_coefficient(formula, Val(kind), Val(source), Val(target), functor, workspace)` | `pair`, `physical`, the pair's columns of `rho`, `epsilon` and `mu`, `thickness`, `jω`, `frequency`, `media`, `options` and the `destinations` |
| Unified earth return | `EarthImpedance.Formula{:unified}`, `EarthAdmittance.Formula{:unified}` | `EarthAdmittance.source_coefficients(formula, Val(kind), Val(source), Val(target), functor, workspace)`, which returns `(axial, potential)` | as for the earth families, with Unified's per-frequency state |
| Soil frequency dependence | `Earth.FrequencyDependent.FrequencyDependentFormulation` | `earth_material(formula, functor, workspace)` | `material`, `frequency`, `options` |
| Temperature dependence | `Materials.TemperatureDependent.TemperatureDependentFormulation` | `temperature_resistivity(formula, functor, workspace)` | `material`, `temperature`, `options` |
| Equivalent earth | `Earth.EquivalentHomogeneous.AbstractRule` | `equivalent_material(formula, Val(kind), Val(source), Val(target), functor, workspace)` | `rho`, `eps_r`, `mu_r`, `model`, the physical `pair`, `frequency`, `options` |
| Modal decomposition | `AbstractFormulation`, selected by `ModalAnalysisFormulation` | `Commons.initialize_buffers(formula, T, input, plan, common)` and `ModalAnalysis.decompose!(formula, workspace, parameters, options)` | none |
| Local shunt geometry | `Engine.ShuntModelFormulation` | `Engine.internal_shunt_response(formula, design, geometry, T, material_selections, solutions, design_index)`, during blueprint construction | none |
| Pipe applicability | `Engine.PipeImpedanceFormulation` | `validate(design, formula, backend)` admits the topology of `design` or throws. No analytical pipe equation is supplied. | none |

An earth expression returns one coefficient per destination of its calculation. It returns a
number for one destination and a tuple, ordered as the destinations, for several.

### Layered earth

The source and target layers of an earth expression select its method. Layer 1 is air,
and layer 2 is the first earth layer. Methods written for `Val{1}` and `Val{2}` declare a
formula for air and one homogeneous earth. On an earth model with more layers, such a
formula consumes its explicit `equivalent_earth` reduction, or the `:default` reduction
without one. A layer left as a generic `Val{S}` declares the expression for every layer:

```julia
function earth_impedance(formula::MyLayeredImpedance, ::Val{:self}, ::Val{S}, ::Val{S},
        functor, workspace) where {S}
    # functor.input.thickness holds the thicknesses of the layered earth.
end
```

A formula that admits a layer above 2 consumes the whole layered earth and its interfaces.
Before the frequency loop, the plan checks that the formula has an expression for every
conductor pair on the earth it consumes. The `ArgumentError` names the number of earth
layers and the highest layer that the formula admits.

`parameters` are model data and `options` are formulation-owned physical choices
and numerical controls. Custom constructors validate and normalize their own `FormulationOptions`.
Execution controls use `ComputationOptions`. Completed supplemental output uses
`ComputationDetails`. Read their payloads explicitly through `.data`.
Indexed families declare numerical defaults for the actual selected type and
case with `formulation_options(::Expression{<:MyType,typeof(operation),...})`.
Earth field equations consume evaluated material properties and define their
own wave numbers and field approximations. Material-law families implement
`constitutive` with a valid material argument.
Unified's formulation options accept a prescribed longitudinal coefficient Γ
as a scalar or a frequency-aligned vector. Its precision participates in
allocation of the computation's numerical storage. The positive-time convention
and explicit convention conversions are documented with the Unified equations.

Coaxial equations and material laws receive the defining computation workspace,
or `nothing` for a standalone evaluation that does not require it. Numerical arrays are in
`workspace.buffers`. The workspace input and plan remain the authority for geometry and
material mappings.

The workspace provides the buffers, and formulas use them. `initialize_buffers` is the
one action that builds buffers. Each equation builds its own in its method, usually as a
named tuple of preallocated arrays, without a storage type of its own. An allocation
method is unnecessary for an algebraic equation. For an equation that needs buffers, extend
`Commons.initialize_buffers(selected, T, input, plan, buffers)` to return
the record extended with owned arrays only. Only selections reached by the
required indexed calls participate in initialization, before evaluating materials
or equations. Unused recipe branches remain unallocated.
The default uses only existing storage. Existing arrays may not be replaced.
Numerical formulas provision the common quadrature storage through
`Commons.initialize_buffers(SpectralIntegral, Val(:quad), T, input, plan, buffers)`. An
integration option is not a capability declaration. Cable constants uses this
same buffer-initialization method with its local numerical input and an empty plan.

Each formula defines its complete integrand, transformations, Jacobians, branch
choices and physical subdivision hints. `SpectralIntegral` selects a Sommerfeld-type
integral over the spatial Fourier variable and contains only that callable.
`integrate(integral, Val(:quad), controls, buffers)` passes it and the supplied numeric
subdivision points to QuadGK. Physical expressions and subdivision choices belong to each formula.

`ShuntModel` owns conductor and dielectric geometry extraction, numerical coefficients and the requested fallback
model. Its `blueprint_dependencies` methods identify the actual local selections
that affect those coefficients. Engine uses that dependency record for reuse
within one blueprint construction. Completed coefficient blocks are retained
after construction. The charge-system factorization and workspace are reused
only during construction.

Conductor pairs with the same inputs share one computed value when their media agree.
By default the pairs' destination indices take part, so distinct pairs are computed
separately. A formula whose arithmetic does not read the indices declares the inputs it
reads with `same_physical_state(formula, a::EarthPair, b::EarthPair, geometry)`. A formula
that computes the whole system at once, such as Unified, adds an `Engine.EarthPlan`
constructor method on its own type and an `Engine.earth!(formula, functor, workspace)`
method that converts the coefficients of its parts into the physical matrices. Each
selection uses the Engine frequency sequence and its computation workspace.

An internal-impedance formula builds the values that its inner, outer and transfer impedances
share in its `Functor` method, once per conductor and frequency. Shunt geometry construction runs once before the
frequency loop and returns frequency-independent blueprint blocks and diagnostics. A
consuming earth formula explicitly admits a custom equivalent-earth rule using `validate`.

Results use result-type and unit checks for each formula family, including custom types.
Impedances are in Ω/m, material admittivities in S/m, earth potential coefficients
in m/F, and temperature-law resistivities in Ω·m. Scalar material laws return
`EarthMaterial` or a finite scalar as appropriate. Modal equations write
`workspace.Tv`, `workspace.Ti`, and `workspace.roots`. The owner constructs
`ModalOperators` and intrinsic coefficients after decomposition. Supply
`formula_id`, `description`, `NamedTuple` and
`formulation_options` methods for metadata. Serialized identities and data do
not reconstruct executable methods. See the [temperature-law example](engine.md#Cable-material-temperature-dependence).

### A custom modal formula

The following complete one-mode example implements a custom modal formula. Its
methods take the formula object first, `Formula{:diagonal_example}`. The modal workspace supplies common slices, coordinate conversion buffers, the
admittance-impedance product, eigenpair history and voltage-vector buffers.
`initialize_buffers` extends that record only for additional numerical work.
It assumes a completed one-mode phase scan named `phase`, with nonzero
diagonal coefficients and known source length. The example is algebraic. It
does not replace a broadband modal model.

```julia
using LineCableModels
import LineCableModels.Engine: description
import LineCableModels.Commons: formulation_options, FormulationOptions, initialize_buffers
import LineCableModels.ModalAnalysis: decompose!, Formula
import LineCableModels: Expression

description(::Type{<:Formula{:diagonal_example}}; compact=false) =
    compact ? "diagonal example" : "one-mode diagonal example"
formulation_options(::Expression{<:Formula{:diagonal_example},typeof(decompose!)}) =
    FormulationOptions()

function initialize_buffers(::Formula{:diagonal_example}, ::Type{T}, input,
        plan, common) where {T<:Complex}
    plan.n == 1 || throw(DimensionMismatch("diagonal example requires one mode"))
    return merge(common, (diagonal_product=Vector{T}(undef,plan.nf),))
end

function decompose!(::Formula{:diagonal_example}, workspace,
        parameters::NamedTuple, options::FormulationOptions)
    product=workspace.buffers.diagonal_product
    for k in eachindex(product)
        z=workspace.input.Z[1,1,k]/workspace.input.root_scale
        y=workspace.input.Y[1,1,k]/workspace.input.root_scale
        product[k]=z*y
        root=sqrt(product[k])
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

selected=ModalAnalysisFormulation(Formula(:diagonal_example))
modal=compute(ModalAnalysisProblem(phase),selected)
segment=PropagationParameters(modal)
size(gamma(modal)) == (1,length(frequencies(phase)))
size(H(segment)) == size(gamma(modal))
```

Only the selected formula's `initialize_buffers` method runs. The common
workspace supplies per-frequency normalization and coordinate buffers. The
formula owns `diagonal_product`. It writes phase-row by mode-column bases
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
checks the completed value directly and returns it unchanged. [`validate`](@ref)
documents the rules that every method follows.

## Commons, observations, and units

```@autodocs
Modules = [
    LineCableModels.Commons,
    LineCableModels.Units,
]
Order = [:module, :constant, :type, :function, :macro]
Filter = developer_reference_entry
Public = true
Private = false
```

## Data model and computations

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
