# Engine-owned formulation hierarchy.
"""
$(TYPEDEF)

Select the LineCableModels backend for concentric coaxial cable assemblies.

Series impedance uses the equivalent concentric representation supplied by
DataModel. The default local shunt model uses coaxial annuli. Explicit
`shunt_model=:boundary` resolves eligible open wires and finite tapes in a
lossless, radially layered circular shielded domain during blueprint
construction. Other geometry and material selections retain their
equivalent-coaxial treatment.
"""
struct LineCableModelsCoaxial end

"""Identify the coaxial backend without executing or configuring a calculation."""
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
Select conductor surface impedance equations [Ω/m]. Concrete subtypes prepare
shared state when called with conductor dimensions and material properties,
and implement `InternalImpedance.internal_impedance` for their supported
`inner`, `outer`, and `transfer` cases.
"""
abstract type InternalImpedanceFormulation <: AbstractImpedanceFormulation end
"""
Select the pipe contribution and its backend/topology applicability. Concrete
subtypes extend `Formulation(backend, selected, Val(topology))`; admission alone
does not supply an impedance equation.
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

"""Supertype for local shunt geometry and dielectric/earth admittance selections."""
abstract type AbstractAdmittanceFormulation <: AbstractFormulation end
"""
Select cable-local shunt geometry, independently of material admittivity.
Concrete subtypes implement `internal_shunt_response` during blueprint construction.
"""
abstract type ShuntModelFormulation <: AbstractAdmittanceFormulation end
"""
Select insulation admittivity [S/m]. Concrete subtypes implement
`InsulationAdmittance.insulation_material`; radial geometry is applied separately.
"""
abstract type InsulationAdmittanceFormulation <: AbstractAdmittanceFormulation end
"""
Select semiconducting-layer admittivity [S/m]. Concrete subtypes implement
`SemiconAdmittance.semicon_material`; radial geometry is applied separately.
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
Validate a dielectric law's admittivity [S/m] at its material boundary and
represent it as `Complex{T}` for radial aggregation. A result requiring a wider
scalar type is rejected instead of silently discarding precision or uncertainty.
"""
function validate(::Union{InsulationAdmittanceFormulation, SemiconAdmittanceFormulation},
        ::Type{T}, value) where {T <: Real}
    value isa Number && !(value isa Bool) && isfinite(value) || throw(DomainError(
        value, "a dielectric material law must return a finite scalar admittivity [S/m]"))
    promote_type(T, typeof(real(value))) === T || throw(ArgumentError(
        "dielectric admittivity requires scalar type $(typeof(real(value))); " *
        "use a material and problem scalar type that preserves its precision and uncertainty"))
    return convert(Complex{T}, value)
end

"""
Return whether an earth formulation consumes homogeneous or stratified media.
"""
function media end

"""
Resolve scalar or named selections through the child slots declared by their formula owner.
"""
function Formulation(::Type{F},
        selected) where {F <: AbstractFormulation}
    return F(selected)
end

function Formulation(::Type{F},
        selected::NamedTuple) where {F <: AbstractFormulation}
    children = (; pairs(F)...)
    all(in(keys(children)), keys(selected)) || throw(ArgumentError(
        "$F selections admit only $(join(keys(children), ", "))"))
    names = filter(in(keys(selected)), keys(children))
    # An explicit recipe stays explicit: a missing or `nothing` leaf supplies
    # no equation. Only omission of the whole family chooses its default.
    return NamedTuple{names}(map(names) do name
        value = selected[name]
        value === nothing ? nothing : children[name](value)
    end)
end

function Formulation(::Type{F}, ::Nothing) where {
        F <: Union{EarthImpedanceFormulation, EarthAdmittanceFormulation,
            InternalImpedanceFormulation}}
    return F(:default)
end

function Formulation(selected::NamedTuple, ::Val{Kind}) where {Kind}
    value = get(selected, Kind, nothing)
    value === nothing && throw(ArgumentError("explicit equation recipe has no requested :$Kind case"))
    return value
end

"""
Resolve the selected formula for exact source and target layer indices.
"""
function Formulation(
        selected::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation},
        ::Val{S}, ::Val{T}) where {S, T}
    selected
end

Formulation(selected::NamedTuple, ::Val{1}, ::Val{1}) = Formulation(selected, Val(:air))
function Formulation(selected::NamedTuple, ::Val{2}, ::Val{2})
    Formulation(selected, Val(:earth))
end
function Formulation(selected::NamedTuple, ::Val{1}, ::Val{2})
    Formulation(selected, Val(:mixed))
end
function Formulation(selected::NamedTuple, ::Val{2}, ::Val{1})
    Formulation(selected, Val(:mixed))
end

function Formulation(::NamedTuple, ::Val{S}, ::Val{T}) where {S, T}
    throw(ArgumentError(
        "homogeneous selection is not defined for source in layer $S and target in layer $T"))
end

function validate(selected::NamedTuple, earth::EarthModel)
    validate(earth)
    return selected
end

function validate(formula::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation},
        pair::EarthPair)
    return only(validate(formula, (pair,)))
end

function validate(formula::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation},
        pairs::Union{Tuple, AbstractVector{<:EarthPair}})
    equations = map(pairs) do pair
        validate(pair)
        equation = FormulaMethod(formula, pair)
        validate(pair, equation)
        equation
    end
    identities = unique(equations)
    bindings = map(identities) do equation
        defaults = formulation_options(equation)
        (equation = equation, kind = typeof(first(equation.arguments)).parameters[1],
            defaults = defaults)
    end
    admitted = union((keys(binding.defaults.data) for binding in bindings)...)
    unknown = setdiff(keys(formula.options.data), admitted)
    isempty(unknown) || throw(ArgumentError(
        "unused formulation options $(Tuple(unknown)) for required cases of :$(formula_id(formula))"))
    resolved = map(bindings) do binding
        names = Tuple(intersect(keys(formula.options.data), keys(binding.defaults.data)))
        options = formulation_options(binding.equation, binding.defaults,
            FormulationOptions(formula.options.data[names]))
        (equation = binding.equation, kind = binding.kind, options = options)
    end
    return map(equation -> resolved[findfirst(==(equation), identities)], equations)
end

# Equation-specific geometric restrictions extend the existing validation protocol.
validate(pair::EarthPair, ::FormulaMethod) = pair

function validate(formula::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation}, count::Integer)
    count in formula.assumptions.layers || throw(DimensionMismatch(
        "formula :$(formula_id(formula)) requires $(formula.assumptions.layers) media including air; received $count"))
    return formula
end

function validate(formula::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation}, earth::EarthModel)
    validate(earth)
    earth.vertical_layers &&
        throw(ArgumentError("earth-return equations require horizontal interfaces or an explicit EquivalentHomogeneous reduction"))
    validate(formula, length(earth.layers))
    return formula
end

function validate(formula::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation},
        rho::AbstractVector, epsilon::AbstractVector, mu::AbstractVector, thickness)
    length(rho) == length(epsilon) == length(mu) ||
        throw(DimensionMismatch("material vectors must align"))
    validate(formula, length(rho))
    all(x -> x > 0 && !isnan(x), rho) ||
        throw(DomainError(rho, "resistivities must be positive, including infinite air resistivity"))
    all(x -> isfinite(x) && !iszero(x), epsilon) && all(x -> isfinite(x) && x > 0, mu) ||
        throw(DomainError((epsilon, mu),
            "permittivities must be nonzero and finite; permeabilities positive and finite"))
    epsilon[1] > 0 || throw(DomainError(epsilon[1], "air permittivity must be positive"))
    restriction = formula.assumptions.permittivity
    restriction in (:positive, :nonzero) ||
        throw(ArgumentError("unknown source permittivity restriction"))
    restriction === :positive && !all(>(0), epsilon) &&
        throw(DomainError(epsilon,
            "formula :$(formula_id(formula)) requires positive permittivity; the earth-material data type permits artificial negative values"))
    if thickness === nothing
        media(formula) === Val(:stratified) && length(rho) > 2 &&
            throw(DimensionMismatch("stratified equations require aligned physical layer thicknesses"))
    else
        length(thickness) == length(rho) ||
            throw(DimensionMismatch("layer thicknesses must align with materials"))
        isinf(first(thickness)) && isinf(last(thickness)) &&
        all(x -> isfinite(x) && x > 0, @view(thickness[2:(end - 1)])) ||
            throw(DomainError(thickness,
                "air and bottom half-spaces must be infinite; internal layers positive and finite"))
        media(formula) === Val(:homogeneous) && length(thickness) != 2 &&
            throw(DimensionMismatch("homogeneous equations have no internal soil interfaces"))
    end
    return formula
end

"""
Route an explicit external formulation tag to its `Val` dispatch method.
"""
Formulation(backend::Symbol; kwargs...) = Formulation(Val(backend); kwargs...)

function validate(binding::FormulaMethod, reduction::EquivalentHomogeneous.AbstractRule)
    throw(ArgumentError("$binding does not admit equivalent-earth reduction :$(formula_id(reduction))"))
end
