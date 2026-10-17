"""
$(TYPEDEF)

Select cylindrical surface equations and their model and numerical controls.
The actual conductor geometry determines which surfaces are required.

$(TYPEDFIELDS)
"""
struct Formula{ID, P <: NamedTuple, O <: FormulationOptions} <: InternalImpedanceFormulation
    "Resolved physical model parameters."
    parameters::P
    "Normalized numerical sections indexed by surface kind."
    options::O
end

"""Evaluate a selected cylindrical surface coefficient in Ω/m."""
function internal_impedance end

"""Evaluate required cylindrical surface impedances in Ω/m."""
function surface_impedances end

"""
$(TYPEDSIGNATURES)

Construct an internal-impedance formulation. Numerical controls are projected
onto its surface equations. Unknown numerical sections are rejected.
A custom formulation subtypes `InternalImpedanceFormulation` and extends
`internal_impedance` on its own type. Its `Functor` method builds the values shared by its
surface impedances.
"""
function Formula{ID}(; parameters::NamedTuple=(;), options::Union{NamedTuple, FormulationOptions} = FormulationOptions()) where {ID}
    ID in formulas(Formula) || throw(ArgumentError("unknown internal-impedance formula :$ID"))
    options = options isa NamedTuple ? FormulationOptions(options) : options
    isempty(parameters) || throw(ArgumentError("internal impedance :$ID has no model parameters"))
    kinds = (:inner, :outer, :transfer)
    selected = Formula{ID, typeof(parameters), typeof(options)}(parameters, options)
    expressions = map(kind -> Expression(selected, internal_impedance, Val(kind)), kinds)
    projected = formulation_options(selected, expressions)
    normalized = FormulationOptions(NamedTuple{kinds}(map(expressions) do expression
        projected.options[findfirst(==(expression), projected.expressions)].data
    end))
    return Formula{ID, typeof(parameters), typeof(normalized)}(parameters, normalized)
end

"""
$(TYPEDSIGNATURES)

Evaluate cylindrical surface coefficients in Ω/m. The wall operator is
`[inner transfer; transfer outer]` in the surface-current basis
`(-enclosed axial current, total axial current including the wall)`.
Assemblers own basis transformation and matrix placement.

# Arguments

- `formula`: the selected formulation.
- `r_in`, `r_ex`: inner and outer conductor radii \\[m\\].
- `rho`: conductor resistivity \\[Ω·m\\].
- `mu_r`: relative permeability \\[dimensionless\\].
- `jω`: imaginary angular frequency \\[1/s\\].
- `workspace`: optional computation workspace passed to each selected equation.

# Returns

- NamedTuple of `outer` for a solid primitive, or `inner`, `outer`, and
  `transfer` for a tubular primitive \\[Ω/m\\].
"""
function surface_impedances(formula::InternalImpedanceFormulation, r_in, r_ex, rho, mu_r, jω;
        workspace=nothing)
    return surface_impedances(formula,
        r_in > 0 ? Val((:inner,:outer,:transfer)) : Val((:outer,)),
        r_in, r_ex, rho, mu_r, jω; workspace)
end

"""
Evaluate the surface impedances that `Kinds` names. The formula's `Functor` stores the values
shared by its surface impedances, built once per conductor and frequency. Each surface
impedance evaluates with the options of its own section.
"""
@inline function surface_impedances(formula::InternalImpedanceFormulation, ::Val{Kinds},
        r_in, r_ex, rho, mu_r, jω; workspace=nothing) where {Kinds}
    length(Kinds) in (1, 3) ||
        throw(ArgumentError("internal surfaces require one or three kinds"))
    functor = Functor(formula, (; r_in, r_ex, rho, mu_r, jω); workspace)
    values = map(map(Val, Kinds)) do kind
        options = formulation_options(formula, kind)
        value = Expression(formula, internal_impedance, kind)(
            Functor(functor, (; options)), workspace)
        value isa Number && isfinite(value) || throw(DomainError(value,
            "internal_impedance must return a finite surface coefficient [Ω/m]"))
        value
    end
    return NamedTuple{Kinds}(values)
end

Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)
Formula(::Val{ID}; kwargs...) where {ID} = Formula{ID}(; kwargs...)
Formula(selected::InternalImpedanceFormulation) = selected
Formula(::Nothing) = Formula(:default)

function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    selection.equivalent_earth === nothing || throw(ArgumentError(
        "equivalent_earth applies only to external earth formulas"))
    return Formula{ID}(; parameters=selection.parameters, options=selection.options)
end

"""Return the stable identifier of a selected internal-impedance formulation."""
formula_id(::Formula{ID}) where {ID} = ID
formula_id(::Type{<:Formula{ID}}) where {ID} = ID
description(value::Formula; compact::Bool=false) = description(typeof(value); compact)
formulation_options(value::Formula) = value.options

"""Expose the selected identity, model parameters and numerical controls."""
Base.NamedTuple(value::Formula) = (identifier=formula_id(value),
    parameters=value.parameters, options=value.options.data)

Base.pairs(::Type{<:Formula}; quantity=nothing) = pairs((;))
