"""
$(TYPEDEF)

Own one earth potential coefficient formulation, its indexed equation declarations and explicit
controls. Each interaction is selected by its explicit kind and source-layer and
target-layer equation signature. The layers in the signatures of its
`earth_potential_coefficient` methods declare the media it handles: methods up to layer 2 treat air
and one homogeneous earth. `parameters` stores model data. `options` stores formulation
choices and numerical controls.

$(TYPEDFIELDS)
"""
struct Formula{ID, P <: NamedTuple, O <: FormulationOptions, E} <:
       EarthAdmittanceFormulation
    "Explicit physical model parameters."
    parameters::P
    "Formulation options. Projected onto required indexed equations during initialization."
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
    "Selected equation and its exact source-target geometry."
    binding::B
    "Evaluated material quantities and angular frequency."
    state::S
    "Normalized formulation options of the selected equation."
    options::O
end

"""
Evaluate one source-owned earth potential coefficient equation.
"""
function earth_potential_coefficient end

"""
$(TYPEDSIGNATURES)

Resolve a formulation's model parameters and numerical controls. Its concrete
selection type defines the indexed equation. Material properties are evaluated
before equation execution.
Formula-specific arguments are normalized by their selected owner.
A missing physical case is unsupported.
"""
function Formula{ID}(; parameters::NamedTuple = (;),
        options::Union{NamedTuple, FormulationOptions} = FormulationOptions(), equivalent_earth = nothing) where {ID}
    options = options isa NamedTuple ? FormulationOptions(options) : options
    options = formulation_options(Formula{ID}, options)
    isempty(parameters) ||
        throw(ArgumentError("earth formula :$ID has no configurable physical parameters"))
    ID in formulas(Formula) || throw(ArgumentError("unknown earthadmittance formula :$ID"))
    reduction = if equivalent_earth === nothing
        nothing
    elseif equivalent_earth isa EquivalentHomogeneous.AbstractSequence
        equivalent_earth
    elseif equivalent_earth isa FormulaDefinition
        EquivalentHomogeneous.AbstractSequence(equivalent_earth)
    else
        throw(ArgumentError("equivalent_earth must be a formula definition or explicit reduction sequence"))
    end
    return Formula{ID, typeof(parameters), typeof(options), typeof(reduction)}(
        parameters, options, reduction)
end

"""
$(TYPEDSIGNATURES)

Construct one validated indexed interaction at fixed angular frequency. Material
vectors list air and every physical earth layer when layer thicknesses are given,
and exactly `(air,soil)` otherwise. Explicit reduction and its physical and effective
mapping are owned by the computation workspace.

# Arguments

- `resistivity`: aligned resistivities [Ω·m].
- `permittivity`: absolute permittivities [F/m].
- `permeability`: absolute permeabilities [H/m].
- `jω`: imaginary angular frequency [1/s].
- `pair`: indexed conductor interaction. Lengths [m].

# Keywords

- `thickness`: aligned layer thicknesses [m] for a layered earth.
"""
function (formula::EarthAdmittanceFormulation)(
        resistivity::AbstractVector{T}, permittivity::AbstractVector{T},
        permeability::AbstractVector{T}, jω::Complex{T}, pair::EarthPair;
        thickness = nothing, physical_pair = pair
) where {T <: Real}
    selected = only(bindings(formula, (pair,)))
    return formula(resistivity, permittivity, permeability, jω, pair, selected;
        thickness, physical_pair)
end

# Evaluate a declaration already bound by bindings() before the frequency loop.
function (formula::EarthAdmittanceFormulation)(
        resistivity::AbstractVector{T}, permittivity::AbstractVector{T},
        permeability::AbstractVector{T}, jω::Complex{T}, pair::EarthPair, selected;
        thickness = nothing, physical_pair = pair
) where {T <: Real}
    isfinite(jω) && !iszero(jω) || throw(DomainError(jω, "jω must be finite and nonzero"))
    options = selected.options
    validate(resistivity, formula, permittivity, permeability, thickness)
    thickness === nothing || validate(pair, thickness)
    layers = thickness === nothing ? (1, 2) : eachindex(permeability)
    σ = map(layer -> conductivity(resistivity[layer]), layers)
    materials = (rho = resistivity, epsilon = permittivity,
        mu = permeability, sigma = σ)
    state = merge(materials, (; jω, thickness))
    binding = (pair = pair, physical_pair = physical_pair,
        kind = selected.kind, expression = selected.expression)
    return Functor(binding, state, options)
end

function (functor::Functor)(workspace = nothing)
    pair = functor.binding.pair
    expression = functor.binding.expression
    values = expression isa Tuple ?
             map(method -> method(functor, pair, workspace), expression) :
             (expression(functor, pair, workspace),)
    all(value -> value isa Number && isfinite(value), values) || throw(DomainError(values,
        "earth coefficients must be finite scalars"))
    converted = map(value -> oftype(functor.state.jω, value), values)
    return expression isa Tuple ? converted : only(converted)
end

function (selected::EarthAdmittanceFormulation)(materials, binding, workspace, frequency::Int)
    state = (jω = workspace.input.jω[frequency], materials, thickness = materials.thickness)
    return (coefficients = (workspace.buffers.Pearth,), state)
end

function (selected::EarthAdmittanceFormulation)(state::NamedTuple, interaction::NamedTuple, declaration)
    index = interaction.index
    return selected(@view(state.materials.rho[:, index]),
        @view(state.materials.epsilon[:, index]), @view(state.materials.mu[:, index]),
        state.jω, interaction.pair, declaration; thickness = state.thickness,
        physical_pair = interaction.physical_pair)
end

function earth!(::EarthAdmittanceFormulation, calculation, workspace)
    (potential = only(calculation.coefficients),)
end

Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)
Formula(::Val{ID}; kwargs...) where {ID} = Formula{ID}(; kwargs...)
Formula(selected::EarthAdmittanceFormulation) = selected
Formula(::Nothing) = Formula(:default)

# A recipe selects one formula per earth route. Within it, a `nothing` route supplies no
# expression. Only omitting the whole slot selects the default.
function Formula(recipe::NamedTuple)
    routes = (; pairs(Formula)...)
    all(in(keys(routes)), keys(recipe)) || throw(ArgumentError(
        "$Formula selections admit only $(join(keys(routes), ", "))"))
    names = filter(in(keys(recipe)), keys(routes))
    return NamedTuple{names}(map(name -> recipe[name] === nothing ? nothing :
        Formula(recipe[name]), names))
end

function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    return Formula{ID}(; parameters = selection.parameters,
        options = selection.options, equivalent_earth = selection.equivalent_earth)
end

formula_id(::Formula{ID}) where {ID} = ID
formula_id(::Type{<:Formula{ID}}) where {ID} = ID
# Identity-only dispatch also describes retained selections without constructors.
description(value::Formula; compact::Bool = false) = description(typeof(value); compact)
formulation_options(value::Formula) = value.options
# Indexed numerical controls remain deferred to their consuming equations.
formulation_options(::Type{<:Formula}, options::FormulationOptions) = options

"""
$(TYPEDSIGNATURES)

Expose the selected identity, model and numerical controls, and explicit
equivalent-earth reduction as a native record.
"""
function Base.NamedTuple(value::Formula)
    return (identifier = formula_id(value),
        parameters = value.parameters, options = value.options.data,
        equivalent_earth = value.equivalent_earth === nothing ? nothing :
                           NamedTuple(value.equivalent_earth))
end

"""Iterate the independently selectable child slots admitted by this formula family."""
function Base.pairs(::Type{<:Formula}; quantity = nothing)
    pairs((air = Formula, earth = Formula, mixed = Formula))
end
