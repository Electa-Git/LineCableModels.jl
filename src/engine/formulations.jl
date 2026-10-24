# Engine-owned formulation hierarchy.
"""
$(TYPEDEF)

Select the LineCableModels backend for concentric coaxial cable assemblies.

Series impedance uses the equivalent concentric representation supplied by
DataModel. The default local shunt model uses the equivalent annular layer. Explicit
`shunt_model=:boundary` resolves eligible open wires and finite tapes in a
lossless, radially layered circular shielded domain during blueprint
construction. Other geometry and material selections retain their
equivalent-coaxial treatment.
"""
struct LineCableModelsCoaxial end

"""Identify the coaxial backend without executing or configuring a computation."""
description(::Type{LineCableModelsCoaxial}; compact::Bool = false) = "coaxial"
function description(::LineCableModelsCoaxial; compact::Bool = false)
    description(LineCableModelsCoaxial; compact)
end
formula_id(::Type{LineCableModelsCoaxial}) = :coaxial
formula_id(::LineCableModelsCoaxial) = :coaxial

"""
$(TYPEDEF)

Supertype for Engine impedance formulations.
"""
abstract type AbstractImpedanceFormulation <: AbstractFormulation end
"""
Select equations for the surface impedance of conductors [Ω/m]. Concrete subtypes build, in a
`Functor` method, the values that inner, outer and transfer surface impedances share, from
the conductor dimensions and material properties. They implement
`InternalImpedance.internal_impedance` for their supported `inner`, `outer`, and `transfer`
cases.
"""
abstract type InternalImpedanceFormulation <: AbstractImpedanceFormulation end
"""
Select the pipe contribution. The topology rule of `DataModel.flatten` rejects pipe-type
cables. A pipe formulation applies to the designs that flatten.
"""
abstract type PipeImpedanceFormulation <: AbstractImpedanceFormulation end
"""
Select magnetic impedance across an insulation annulus [Ω/m]. Concrete subtypes
implement `InsulationImpedance.insulation_impedance` on their selected type.
"""
abstract type InsulationImpedanceFormulation <: AbstractImpedanceFormulation end
"""
Select earth-return impedance equations [Ω/m]. Concrete subtypes implement
indexed `EarthImpedance.earth_impedance` methods on evaluated material inputs.
Air, earth, and mixed selections refer to conductor locations.
"""
abstract type EarthImpedanceFormulation <: AbstractImpedanceFormulation end

"""Supertype for local shunt geometry and dielectric and earth admittance selections."""
abstract type AbstractAdmittanceFormulation <: AbstractFormulation end
"""
Select cable-local shunt geometry, independently of material admittivity.
Concrete subtypes implement `ShuntModel.dielectric_groups`, which the shunt stage calls
after flattening, and `ShuntModel.internal_capacitance` when they return groups to solve.
"""
abstract type ShuntModelFormulation <: AbstractAdmittanceFormulation end
"""
Select insulation admittivity [S/m]. Concrete subtypes implement
`InsulationAdmittance.insulation_material`. Radial geometry is applied separately.
"""
abstract type InsulationAdmittanceFormulation <: AbstractAdmittanceFormulation end
"""
Select semiconducting-layer admittivity [S/m]. Concrete subtypes implement
`SemiconAdmittance.semicon_material`. Radial geometry is applied separately.
"""
abstract type SemiconAdmittanceFormulation <: AbstractAdmittanceFormulation end
"""
Select earth potential-coefficient equations [m/F]. Concrete subtypes implement
indexed `EarthAdmittance.earth_potential_coefficient` methods on evaluated
material inputs. Matrix assembly converts
the potential coefficients to shunt admittance [S/m].
"""
abstract type EarthAdmittanceFormulation <: AbstractAdmittanceFormulation end

"""
Validate the admittivity [S/m] that a dielectric law returns, before radial aggregation
represents it as `Complex{T}`. A result requiring a wider scalar type than `T` is rejected
instead of silently discarding precision or uncertainty. Return `value`.
"""
function validate(value, ::Union{InsulationAdmittanceFormulation, SemiconAdmittanceFormulation},
        ::Type{T}) where {T <: Real}
    value isa Number && !(value isa Bool) && isfinite(value) || throw(DomainError(
        value, "a dielectric material law must return a finite scalar admittivity [S/m]"))
    promote_type(T, typeof(real(value))) === T || throw(ArgumentError(
        "dielectric admittivity requires scalar type $(typeof(real(value))); " *
        "use a material and problem scalar type that preserves its precision and uncertainty"))
    return value
end

# The expression of an earth formula for one interaction: its kind and its source and target
# layers.
function Expression(formula::EarthImpedanceFormulation, pair::EarthPair)
    return Expression(formula, EarthImpedance.earth_impedance,
        Val(pair.row == pair.column ? :self : :mutual), Val.(layer_index(pair))...)
end
function Expression(formula::EarthAdmittanceFormulation, pair::EarthPair)
    return Expression(formula, EarthAdmittance.earth_potential_coefficient,
        Val(pair.row == pair.column ? :self : :mutual), Val.(layer_index(pair))...)
end

"""
$(TYPEDSIGNATURES)

Build the expression of an earth recipe for `pair`: the `air` formula for layers (1, 1), the
`earth` formula for (2, 2) and the `mixed` formula for (1, 2) and (2, 1). The pair's layers
are the decided ones: a recipe describes a homogeneous earth.

# Errors

- Throws `ArgumentError` when the recipe has no formula for the pair's route, or for another
  layer pair.
"""
function Expression(recipe::NamedTuple, pair::EarthPair)
    source, target = layer_index(pair)
    route = (source, target) == (1, 1) ? :air : (source, target) == (2, 2) ? :earth :
            (source, target) in ((1, 2), (2, 1)) ? :mixed :
            throw(ArgumentError("homogeneous selection is not defined for source in layer " *
                "$source and target in layer $target"))
    selected = get(recipe, route, nothing)
    selected === nothing &&
        throw(ArgumentError("explicit equation recipe has no requested :$route case"))
    return Expression(selected, pair)
end

# Whether the formula of `expression` admits earth layer `k`. It does when a method of its
# operation accepts `Val{k}` in the source or the target position, with the `Functor` and the
# workspace that the evaluation passes, whatever the types of the other arguments. A layer
# left generic, such as `::Val{S}`, admits every layer. The earth decision, the error message
# of `validate(expression, layers)` and the standalone formula call read the layers here.
function _admits_layer(expression::Expression, k::Int)
    F = typeof(expression.selection)
    return !isempty(methods(expression.method, Tuple{F, Any, Val{k}, Any, Functor, Any})) ||
           !isempty(methods(expression.method, Tuple{F, Any, Any, Val{k}, Functor, Any}))
end

"""
$(TYPEDSIGNATURES)

Check that the formula of `expression` defines it for an earth with one entry of `layers` per
medium, before evaluating it: the layers of an earth model, or the resistivities of the
standalone formula call. The signatures of a formula's methods declare the earth layers it
handles, whatever the types of their runtime arguments. In layers that the formula admits,
`validate(expression)` checks the expression. Return `expression`.

# Errors

- Throws `ArgumentError` with the layer count N and the highest layer up to N that the
  formula admits, when the formula does not admit the source or target layer.
- Throws the `ArgumentError` of `validate(expression)` for a missing expression in layers
  that the formula admits.

When the formula admits layer N + 1, the message reads "defined beyond layer N" instead.
"""
function validate(expression::Expression{<:Union{EarthImpedanceFormulation, EarthAdmittanceFormulation}},
        layers::Union{Tuple, AbstractVector})
    selected = expression.selection
    signature = Tuple{typeof(selected), map(typeof, expression.arguments)..., Functor, Any}
    hasmethod(expression.method, signature) && return expression
    interaction(::Val{K}, ::Val{S}, ::Val{T}) where {K, S, T} = (K, S, T)
    kind, source, target = interaction(expression.arguments...)
    if _admits_layer(expression, source) && _admits_layer(expression, target)
        validate(expression)
        return expression
    end
    N = length(layers)
    defined = _admits_layer(expression, N + 1) ? "defined beyond layer $N" :
              "defined up to layer $(something(findlast(k -> _admits_layer(expression, k), 1:N), 0))"
    throw(ArgumentError(
        "the earth model has $N layers and formula " *
        ":$(formula_id(selected)) is $defined; it has no expression " *
        "for a $kind interaction from layer $source to layer $target"))
end

function validate(earth::EarthModel, formula::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation})
    earth.vertical_layers &&
        throw(ArgumentError("earth-return equations require horizontal interfaces or an explicit EquivalentHomogeneous reduction"))
    return earth
end

function validate(rho::AbstractVector,
        formula::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation},
        epsilon::AbstractVector, mu::AbstractVector, thickness)
    length(rho) == length(epsilon) == length(mu) ||
        throw(DimensionMismatch("material vectors must align"))
    all(x -> x > 0 && !isnan(x), rho) ||
        throw(DomainError(rho, "resistivities must be positive, including infinite air resistivity"))
    all(x -> isfinite(x) && !iszero(x), epsilon) && all(x -> isfinite(x) && x > 0, mu) ||
        throw(DomainError((epsilon, mu),
            "permittivities must be nonzero and finite; permeabilities positive and finite"))
    epsilon[1] > 0 || throw(DomainError(epsilon[1], "air permittivity must be positive"))
    all(>(0), epsilon) ||
        throw(DomainError(epsilon,
            "formula :$(formula_id(formula)) requires positive permittivity; the earth-material data type permits artificial negative values"))
    if thickness === nothing
        length(rho) == 2 || throw(DimensionMismatch(
            "material vectors without layer thicknesses describe air and one earth medium"))
    else
        length(thickness) == length(rho) ||
            throw(DimensionMismatch("layer thicknesses must align with materials"))
        isinf(first(thickness)) && isinf(last(thickness)) &&
        all(x -> isfinite(x) && x > 0, @view(thickness[2:(end - 1)])) ||
            throw(DomainError(thickness,
                "air and bottom half-spaces must be infinite; internal layers positive and finite"))
    end
    return rho
end

"""
$(TYPEDSIGNATURES)

Return the normalized options of the surface impedance `kind` of the internal-impedance
formula `formula`: the section of its options for this kind, or empty options.
"""
function formulation_options(formula::InternalImpedanceFormulation, ::Val{Kind}) where {Kind}
    return FormulationOptions(get(formula.options.data, Kind, (;)))
end

"""
$(TYPEDSIGNATURES)

Build the Functor of an insulation-impedance formula for one insulation annulus at one
frequency, after checking its radii, its relative permeability and `jω`. The state is empty.
"""
function Functor(formula::InsulationImpedanceFormulation, input::NamedTuple;
        workspace = nothing)
    (; r_in, r_ex, mu_r, jω) = input
    T = typeof(r_in)
    isfinite(r_in) && isfinite(r_ex) && zero(T) <= r_in <= r_ex ||
        throw(DomainError((r_in, r_ex), "insulation radii must satisfy 0 ≤ r_in ≤ r_ex [m]"))
    isfinite(mu_r) && mu_r > zero(T) || throw(DomainError(mu_r,
        "relative insulation permeability must be positive and finite"))
    isfinite(jω) || throw(DomainError(jω, "jω must be finite"))
    return Functor(formula, input, (;))
end

"""
$(TYPEDSIGNATURES)

Build the Functor of an insulation admittivity law for one material at one frequency and
temperature, after checking both. The state is empty.
"""
function Functor(formula::InsulationAdmittanceFormulation, input::NamedTuple;
        workspace = nothing)
    (; frequency, temperature) = input
    isfinite(frequency) && frequency > zero(frequency) || throw(DomainError(
        frequency, "insulation constitutive frequency must be positive and finite"))
    isfinite(temperature) || throw(DomainError(
        temperature, "insulation constitutive temperature must be finite"))
    return Functor(formula, input, (;))
end

"""
$(TYPEDSIGNATURES)

Build the Functor of a semiconductor admittivity law for one material at one frequency and
temperature, after checking both. The state is empty.
"""
function Functor(formula::SemiconAdmittanceFormulation, input::NamedTuple;
        workspace = nothing)
    (; frequency, temperature) = input
    isfinite(frequency) && frequency > zero(frequency) || throw(DomainError(
        frequency, "semicon constitutive frequency must be positive and finite"))
    isfinite(temperature) || throw(DomainError(
        temperature, "semicon constitutive temperature must be finite"))
    return Functor(formula, input, (;))
end

"""
$(TYPEDSIGNATURES)

Build the Functor of an earth calculation of the impedance formula `formula` at one
frequency, whose parts write the impedance coefficients of each conductor pair into `Zearth`
of the workspace buffers, with an empty state.
"""
function Functor(formula::EarthImpedanceFormulation, input::NamedTuple; workspace)
    return Functor(formula, merge(input, (destinations = (workspace.buffers.Zearth,),)), (;))
end

"""
$(TYPEDSIGNATURES)

Build the Functor of an earth calculation of the admittance formula `formula` at one
frequency, whose parts write the potential coefficients of each conductor pair into `Pearth`
of the workspace buffers, with an empty state.
"""
function Functor(formula::EarthAdmittanceFormulation, input::NamedTuple; workspace)
    return Functor(formula, merge(input, (destinations = (workspace.buffers.Pearth,),)), (;))
end

"""
$(TYPEDSIGNATURES)

Check the input of one conductor pair before an earth formula evaluates it: a finite nonzero
`jω`, the aligned material vectors of the pair, and, on a layered earth, the pair's heights
against the layer thicknesses. Return `input`.
"""
function validate(input::NamedTuple,
        formula::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation})
    isfinite(input.jω) && !iszero(input.jω) ||
        throw(DomainError(input.jω, "jω must be finite and nonzero"))
    validate(input.rho, formula, input.epsilon, input.mu, input.thickness)
    input.thickness === nothing || validate(input.pair, input.thickness)
    return input
end

"""
$(TYPEDSIGNATURES)

Evaluate the earth coefficient of one indexed conductor pair at fixed angular frequency.
Material vectors list air and every physical earth layer when layer thicknesses are given,
and exactly `(air,soil)` otherwise. The call runs the checks of a computation: the pair and
its geometry, the expression for its kind and layers on an earth with these media, the
formula's options for that expression, and the pair's input. Return the coefficient as
`Complex{T}`.

# Arguments

- `resistivity`: aligned resistivities [Ω·m].
- `permittivity`: absolute permittivities [F/m].
- `permeability`: absolute permeabilities [H/m].
- `jω`: imaginary angular frequency [1/s].
- `pair`: indexed conductor interaction. Lengths [m].

# Keywords

- `thickness`: aligned layer thicknesses [m] for a layered earth.
- `physical_pair`: the physical pair before an equivalent-earth reduction.
- `workspace`: the computation workspace, when the expression reads its buffers.
"""
function (formula::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation})(
        resistivity::AbstractVector{T}, permittivity::AbstractVector{T},
        permeability::AbstractVector{T}, jω::Complex{T}, pair::EarthPair;
        thickness = nothing, physical_pair = pair, workspace = nothing
) where {T <: Real}
    validate(pair)
    expression = Expression(formula, pair)
    validate(expression, resistivity)
    options = only(formulation_options(formula, (expression,)).options)
    # One pair does not need destinations. The plain constructor gives the empty state of a
    # formula that does not share values. Unified, whose state is its whole system, runs only in
    # a computation.
    functor = Functor(formula, (; jω, thickness, pair, physical = physical_pair,
        rho = resistivity, epsilon = permittivity, mu = permeability, options), (;))
    validate(functor.input, formula)
    value = expression(functor, workspace)
    value isa Number && isfinite(value) ||
        throw(DomainError((value,), "earth coefficients must be finite scalars"))
    return oftype(jω, value)
end

"""
$(TYPEDSIGNATURES)

Build the backend formulation whose identity is `tag`: `:coaxial` or `:cable_constants`.
`Formulation(::Val{tag}; kwargs...)` forwards `kwargs` to the backend's constructor. A
package extension adds the tag of its backend, as the PSCAD extension adds `:pscad` when
PythonCall is loaded.

# Errors

- Throws `ArgumentError` for an unknown `tag`.
"""
Formulation(tag::Symbol; kwargs...) = Formulation(Val(tag); kwargs...)

Formulation(::Val{tag}; kwargs...) where {tag} =
    throw(ArgumentError("unknown backend formulation :$tag. A backend held in a package " *
                        "extension loads with its trigger package."))

# An earth equation admits an equivalent-earth reduction only through its own method.
function validate(reduction::EquivalentHomogeneous.AbstractRule,
        expression::Expression{<:Union{EarthImpedanceFormulation, EarthAdmittanceFormulation}})
    throw(ArgumentError("$expression does not admit equivalent-earth reduction :$(formula_id(reduction))"))
end
