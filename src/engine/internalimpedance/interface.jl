"""
$(TYPEDEF)

Select cylindrical surface equations and their model and numerical controls.
The actual conductor geometry determines which surfaces are required.

$(TYPEDFIELDS)
"""
struct Formula{ID, P <: NamedTuple, O <: FormulationOptions, C <: Tuple} <: InternalImpedanceFormulation
    "Resolved physical/model parameters."
    parameters::P
    "Normalized numerical sections indexed by surface kind."
    options::O
    "Explicit numerical section names, checked against actual consumers."
    configured_options::C
end

"""
$(TYPEDEF)

Retain a selected formulation and shared state for one conductor and frequency.
Leaf evaluation uses that selection's equation directly.

$(TYPEDFIELDS)
"""
struct Functor{F, S, O <: FormulationOptions}
    "Selected internal-impedance formulation."
    selection::F
    "Shared conductor state at one frequency."
    state::S
    "Normalized numerical sections for the surfaces or the current leaf."
    options::O
end

"""Return the stable identifier of a selected internal-impedance formulation."""
formula_id(::Formula{ID}) where {ID} = ID

"""Evaluate a selected cylindrical surface coefficient in Ω/m."""
function internal_impedance end

"""Evaluate required cylindrical surface impedances in Ω/m."""
function surface_impedances end

"""
$(TYPEDSIGNATURES)

Construct an internal-impedance formulation. Numerical controls are projected
onto its implemented surfaces; unknown or unused sections are rejected.
Custom formulations subtype `InternalImpedanceFormulation`, supply their
shared-state constructor, and extend `internal_impedance` on their own type.
"""
Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)

function Formula(::Val{ID}; parameters::NamedTuple=(;), options::Union{NamedTuple, FormulationOptions} = FormulationOptions()) where {ID}
    options = options isa NamedTuple ? FormulationOptions(options) : options
    isempty(parameters) || throw(ArgumentError("internal impedance :$ID has no model parameters"))
    kinds = (:inner, :outer, :transfer)
    supplied = keys(options.data)
    selected = Formula{ID, typeof(parameters), typeof(options), typeof(supplied)}(
        parameters, options, supplied)
    declarations = map(kinds) do kind
        binding = FormulaMethod(selected, internal_impedance, Val(kind))
        signature = Tuple{typeof(selected), typeof(Val(kind)), Any, Any}
        which(internal_impedance, signature) === EQUATION_FALLBACK && return nothing
        formulation_options(binding)
    end
    all(isnothing, declarations) && throw(ArgumentError(
        "internal impedance :$ID has no implemented surface equations"))
    admitted = union((keys(value.data) for value in declarations if value !== nothing)...)
    isempty(setdiff(supplied, admitted)) || throw(ArgumentError(
        "unused internal-impedance numerical sections"))
    normalized = FormulationOptions(NamedTuple{kinds}(map(kinds, declarations) do kind, defaults
        defaults === nothing && return (;)
        binding = FormulaMethod(selected, internal_impedance, Val(kind))
        names = Tuple(intersect(supplied, keys(defaults.data)))
        formulation_options(binding, defaults, FormulationOptions(options.data[names])).data
    end))
    return Formula{ID, typeof(parameters), typeof(normalized), typeof(supplied)}(
        parameters, normalized, supplied)
end

Formula(selected::InternalImpedanceFormulation) = selected

function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    selection.equivalent_earth === nothing || throw(ArgumentError(
        "equivalent_earth applies only to external earth formulas"))
    return Formula(Val(ID); parameters=selection.parameters, options=selection.options)
end

function internal_impedance(selected::InternalImpedanceFormulation, ::Val{Kind},
        functor, workspace) where {Kind}
    throw(ArgumentError(
        "internal_impedance :$(formula_id(selected)): formula not implemented for kind :$Kind"))
end

const EQUATION_FALLBACK = which(internal_impedance,
    Tuple{InternalImpedanceFormulation, Val, Any, Any})

function validate(binding::FormulaMethod{S, typeof(internal_impedance), A}) where {S, A}
    which(internal_impedance, Tuple{S, A.parameters..., Any, Any}) === EQUATION_FALLBACK &&
        binding(nothing, nothing)
    return binding
end

function (functor::Functor)(::Val{Kind}, workspace=nothing) where {Kind}
    hasproperty(functor.options.data, Kind) || throw(ArgumentError(
        "internal surface :$Kind is not implemented by this selection"))
    options = getproperty(functor.options.data, Kind)
    leaf = Functor(functor.selection, functor.state, FormulationOptions(options))
    value = internal_impedance(functor.selection, Val(Kind), leaf, workspace)
    value isa Number && isfinite(value) || throw(DomainError(value,
        "internal_impedance must return a finite surface coefficient [Ω/m]"))
    return value
end

(functor::Functor)(kind::Symbol, workspace=nothing) = functor(Val(kind), workspace)

"""
$(TYPEDSIGNATURES)

Evaluate cylindrical surface coefficients in Ω/m. The wall operator is
`[inner transfer; transfer outer]` in the surface-current basis
`(-enclosed axial current, total axial current including the wall)`.
Assemblers own basis transformation and matrix placement.

# Arguments

- `formula`: One formulation or complete `inner/outer/transfer` selections.
- `r_in`, `r_ex`: Inner and outer conductor radii \\[m\\].
- `rho`: Conductor resistivity \\[Ω·m\\].
- `mu_r`: Relative permeability \\[dimensionless\\].
- `jω`: Imaginary angular frequency \\[1/s\\].
- `workspace`: Optional numerical resources passed to each selected equation.

# Returns

- NamedTuple of `inner`, `outer`, and `transfer` coefficients \\[Ω/m\\].
"""
function surface_impedances(formula::Union{InternalImpedanceFormulation,
        NamedTuple{(:inner,:outer,:transfer)}}, r_in, r_ex, rho, mu_r, jω;
        workspace=nothing)
    validate(formula, (:inner, :outer, :transfer))
    return surface_impedances(formula, Val((:inner,:outer,:transfer)),
        r_in, r_ex, rho, mu_r, jω; workspace)
end

"""Evaluate prevalidated surfaces, preparing shared state once per conductor and frequency."""
@inline function surface_impedances(formula::InternalImpedanceFormulation, ::Val{Kinds},
        r_in, r_ex, rho, mu_r, jω; workspace=nothing) where {Kinds}
    functor = formula(r_in, r_ex, rho, mu_r, jω)
    length(Kinds) == 1 && return NamedTuple{Kinds}((functor(Val(Kinds[1]),workspace),))
    length(Kinds) == 3 || throw(ArgumentError("internal surfaces require one or three kinds"))
    return NamedTuple{Kinds}((functor(Val(Kinds[1]),workspace),
        functor(Val(Kinds[2]),workspace), functor(Val(Kinds[3]),workspace)))
end

@inline function surface_impedances(selected::NamedTuple{(:inner,:outer,:transfer)},
        ::Val{Kinds}, r_in, r_ex, rho, mu_r, jω; workspace=nothing) where {Kinds}
    if length(Kinds) == 1
        return surface_impedances(selected[first(Kinds)], Val(Kinds),
            r_in, r_ex, rho, mu_r, jω; workspace)
    end
    length(Kinds) == 3 || throw(ArgumentError("internal surfaces require one or three kinds"))
    a, b, c = selected[Kinds[1]], selected[Kinds[2]], selected[Kinds[3]]
    first_functor = a(r_in, r_ex, rho, mu_r, jω)
    second_functor = b === a ? first_functor : b(r_in, r_ex, rho, mu_r, jω)
    third_functor = c === a ? first_functor : c === b ? second_functor :
                    c(r_in, r_ex, rho, mu_r, jω)
    return NamedTuple{Kinds}((first_functor(Val(Kinds[1]),workspace),
        second_functor(Val(Kinds[2]),workspace), third_functor(Val(Kinds[3]),workspace)))
end

"""Validate each required surface against its own selected formulation and controls."""
function validate(selected::NamedTuple{(:inner,:outer,:transfer)}, kinds::Tuple)
    for (kind, leaf) in pairs(selected)
        if kind in kinds
            validate(leaf, (kind,))
        else
            isempty(leaf.parameters) && isempty(leaf.configured_options) || throw(ArgumentError(
                "explicit internal $kind controls are unused by this assembly"))
        end
    end
    return selected
end

function validate(formula::InternalImpedanceFormulation, kinds::Tuple)
    foreach(kind -> validate(FormulaMethod(formula, internal_impedance, Val(kind))), kinds)
    admitted = union((keys(getproperty(formula.options.data, kind)) for kind in kinds)...)
    isempty(setdiff(formula.configured_options, admitted)) || throw(ArgumentError(
        "an explicitly configured numerical section is unused by the required internal surfaces"))
    return formula
end

"""Expose the selected identity, model parameters and numerical controls."""
Base.NamedTuple(value::Formula) = (identifier=formula_id(value),
    parameters=value.parameters, options=value.options.data, configured_options=value.configured_options)

import ...Grammar: formulation_options
description(value::Formula; compact::Bool=false) = description(typeof(value); compact)
Base.pairs(::Type{<:Formula}; quantity=nothing) = pairs((inner=Formula, outer=Formula, transfer=Formula))
formula_id(::Type{<:Formula{ID}}) where {ID} = ID
formulation_options(value::Formula) = value.options
