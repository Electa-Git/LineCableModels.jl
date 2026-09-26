"""
$(TYPEDEF)

Own one earth impedance formulation, its indexed equation declarations and explicit
controls. `assumptions` contains scientific restrictions; each interaction
is selected by its explicit kind/source-layer/target-layer equation signature.
`parameters` stores model data; `options` stores formulation choices and numerical controls.

$(TYPEDFIELDS)
"""
struct Formula{ID, A <: NamedTuple, P <: NamedTuple, O <: FormulationOptions, E} <:
       EarthImpedanceFormulation
    "Model class and supported physical medium inventory."
    assumptions::A
    "Explicit physical/model parameters."
    parameters::P
    "Formulation options; projected onto required indexed equations during initialization."
    options::O
    "Independent equivalent homogeneous-earth reduction, or nothing."
    equivalent_earth::E
end

"""
$(TYPEDEF)

Bind one indexed interaction to evaluated physical state. Formulation options and
evaluated material quantities remain distinct. The main computation workspace is passed
when evaluating an equation that needs reusable numerical buffers.

$(TYPEDFIELDS)
"""
struct Functor{B, S, O <: FormulationOptions}
    "Selected equation and its exact source/target geometry."
    binding::B
    "Evaluated material quantities and angular frequency."
    state::S
    "Normalized formulation options of the selected equation."
    options::O
end

formula_id(::Formula{ID}) where {ID} = ID
assumptions(formula::EarthImpedanceFormulation) = formula.assumptions
media(formula::EarthImpedanceFormulation) = formula.assumptions.media

"""
Declare model inventory and physical restrictions, without callable behavior.
"""
function assumptions end

function assumptions(::Val{ID}) where {ID}
    throw(ArgumentError("unknown earthimpedance formula :$ID"))
end
"""
Evaluate one source-owned earth impedance equation.
"""
function earth_impedance end

Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)

"""
$(TYPEDSIGNATURES)

Resolve a formulation's model parameters and numerical controls. Its concrete
selection type owns the indexed equation; material properties are evaluated
before equation execution.
Formula-specific arguments are normalized by their selected owner.
A missing physical case is unsupported.
"""
function Formula(::Val{ID}; parameters::NamedTuple = (;),
        options::Union{NamedTuple, FormulationOptions} = FormulationOptions(), equivalent_earth = nothing) where {ID}
    options = options isa NamedTuple ? FormulationOptions(options) : options
    options = formulation_options(Formula{ID}, options)
    parameters = earth_parameters(Val(ID), parameters)
    declared_constraints = assumptions(Val(ID))
    constraints = merge(declared_constraints, (media = Val(declared_constraints.media),))
    reduction = if equivalent_earth === nothing
        nothing
    elseif equivalent_earth isa EquivalentHomogeneous.AbstractSequence
        equivalent_earth
    elseif equivalent_earth isa FormulaDefinition
        EquivalentHomogeneous.AbstractSequence(equivalent_earth)
    else
        throw(ArgumentError("equivalent_earth must be a formula definition or explicit reduction sequence"))
    end
    reduction !== nothing && constraints.media !== Val(:homogeneous) &&
        throw(ArgumentError("a full multilayer formula cannot consume an equivalent homogeneous reduction"))
    return Formula{ID, typeof(constraints), typeof(parameters),
        typeof(options), typeof(reduction)}(
        constraints, parameters, options, reduction)
end

"""
$(TYPEDSIGNATURES)

Prepare one validated indexed interaction at fixed angular frequency. Material
vectors use physical indices for a stratified equation and exactly `(air,soil)`
for a homogeneous equation. Explicit reduction and its physical/effective
mapping are owned by the computation workspace.

# Arguments

- `resistivity`: Aligned resistivities [Ω·m].
- `permittivity`: Absolute permittivities [F/m].
- `permeability`: Absolute permeabilities [H/m].
- `jω`: Imaginary angular frequency [1/s].
- `pair`: Indexed conductor interaction; lengths [m].

# Keywords

- `thickness`: Aligned layer thicknesses [m] for a stratified model.
"""
function (formula::EarthImpedanceFormulation)(
        resistivity::AbstractVector{T}, permittivity::AbstractVector{T},
        permeability::AbstractVector{T}, jω::Complex{T}, pair::EarthPair;
        thickness = nothing, physical_pair = pair
) where {T <: Real}
    selected = validate(formula, pair)
    return formula(resistivity, permittivity, permeability, jω, pair, selected;
        thickness, physical_pair)
end

# Evaluate a declaration already bound by validate() before the frequency loop.
function (formula::EarthImpedanceFormulation)(
        resistivity::AbstractVector{T}, permittivity::AbstractVector{T},
        permeability::AbstractVector{T}, jω::Complex{T}, pair::EarthPair, selected;
        thickness = nothing, physical_pair = pair
) where {T <: Real}
    isfinite(jω) && !iszero(jω) || throw(DomainError(jω, "jω must be finite and nonzero"))
    options = selected.options
    validate(formula, resistivity, permittivity, permeability, thickness)
    thickness === nothing || validate(pair, thickness)
    layers = media(formula) === Val(:homogeneous) ? (1, 2) : eachindex(permeability)
    σ = map(layer -> conductivity(resistivity[layer]), layers)
    materials = (rho = resistivity, epsilon = permittivity,
        mu = permeability, sigma = σ)
    state = merge(materials, (; jω, thickness))
    binding = (pair = pair, physical_pair = physical_pair,
        kind = selected.kind, equation = selected.equation)
    return Functor(binding, state, options)
end

function (functor::Functor)(workspace = nothing)
    pair = functor.binding.pair
    value = functor.binding.equation(functor, pair, workspace)
    value isa Number && isfinite(value) || throw(DomainError(value,
        "earth-impedance contribution must be a finite scalar"))
    return oftype(functor.state.jω, value)
end

function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    return Formula(Val(ID); parameters = selection.parameters,
        options = selection.options, equivalent_earth = selection.equivalent_earth)
end

Formula(selected::EarthImpedanceFormulation) = selected

function FormulaMethod(formula::EarthImpedanceFormulation, pair::EarthPair)
    return FormulaMethod(formula, earth_impedance,
        Val(pair.row == pair.column ? :self : :mutual), Val.(layer_index(pair))...)
end

function earth_impedance(
        selected::EarthImpedanceFormulation, ::Val{Kind}, ::Val{S}, ::Val{T},
        functor, pair, workspace) where {Kind, S, T}
    throw(ArgumentError(
        "earth_impedance :$(formula_id(selected)) ($Kind): formula not implemented for source in layer $S and target in layer $T"))
end


"""
$(TYPEDSIGNATURES)

Expose the selected identity, scientific restrictions, model and numerical
controls, and explicit equivalent-earth reduction as a native record.
"""
function Base.NamedTuple(value::Formula)
    return (identifier = formula_id(value), assumptions = value.assumptions,
        parameters = value.parameters, options = value.options.data,
        equivalent_earth = value.equivalent_earth === nothing ? nothing :
                           NamedTuple(value.equivalent_earth))
end

# Identity-only dispatch also describes retained selections without constructors.
import ...Grammar: formulation_options
description(value::Formula; compact::Bool = false) = description(typeof(value); compact)

"""Iterate the independently selectable child slots admitted by this formula family."""
function Base.pairs(::Type{<:Formula}; quantity = nothing)
    pairs((air = Formula, earth = Formula, mixed = Formula))
end
formula_id(::Type{<:Formula{ID}}) where {ID} = ID
formulation_options(value::Formula) = value.options

# Indexed numerical controls remain deferred to their consuming equations.
formulation_options(::Type{<:Formula}, options::FormulationOptions) = options
