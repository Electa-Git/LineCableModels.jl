"""
$(TYPEDEF)

Define one scalar line-parameter calculation over a completed cable system.
Operating temperature and analysis frequencies are fields of the problem.

$(TYPEDFIELDS)
"""
struct LineParametersProblem{
    T <: Real,
    S <: LineCableSystem{T},
    P <: Union{Nothing, Vector{Complex{T}}}
} <: AbstractProblemDefinition
    "Physical cable system."
    system::S
    "Operating temperature \\[°C\\]."
    temperature::T
    "Static earth model."
    earth_props::EarthModel{T}
    "Strictly positive, sorted analysis frequencies \\[Hz\\]."
    frequencies::Vector{T}
    "Optional longitudinal propagation constants aligned with frequency \\[1/m\\]."
    Γ::P
end

Base.eltype(::LineParametersProblem{T}) where {T} = T
Base.eltype(::Type{LineParametersProblem{T}}) where {T} = T

function _check_line_parameters_problem(problem::LineParametersProblem)
    validate(problem.system)
    validate(problem.earth_props)
    for (design, pose) in zip(problem.system.designs, problem.system.positions)
        iszero(pose.y) && throw(DomainError(
            pose.y,
            "a cable centre cannot lie on the air-earth interface"
        ))
        radius = DataModel.outer_radius(design)
        abs(pose.y) >= radius || throw(DomainError(
            pose.y,
            "the cable cross-section crosses the air-earth interface"
        ))
    end
    phases = unique(problem.system.connection_order)
    positive = filter(>(0), phases)
    isempty(positive) && throw(ArgumentError(
        "at least one conductor must be assigned to a positive phase",
    ))
    maximum(positive) <= nphases(problem.system) || throw(DomainError(
        positive,
        "a phase assignment exceeds the number of distinct positive phases"
    ))
    for design in problem.system.designs
        reference = first(design.geometry.regions).source.material.T0
        abs(problem.temperature - reference) < oftype(problem.temperature, 150) ||
            throw(DomainError(
                problem.temperature,
                "operating temperature is outside the linear resistivity model range"
            ))
    end
    isempty(problem.frequencies) && throw(ArgumentError("frequencies cannot be empty"))
    all(>(zero(eltype(problem))), problem.frequencies) || throw(DomainError(
        problem.frequencies, "frequencies must be positive"
    ))
    issorted(problem.frequencies) || throw(ArgumentError("frequencies must be sorted"))
    if problem.Γ !== nothing
        length(problem.Γ) == length(problem.frequencies) ||
            throw(DimensionMismatch(
                "longitudinal propagation constants must align with frequencies"
            ))
        all(value -> isfinite(real(value)) && isfinite(imag(value)), problem.Γ) ||
            throw(DomainError(
                problem.Γ,
                "longitudinal propagation constants must be finite"
            ))
    end
    maximum(problem.frequencies) > oftype(first(problem.frequencies), 1e8) &&
        @warn("Frequencies above 100 MHz exceed the quasi-TEM validity range.",
            max_frequency=maximum(problem.frequencies),)
    return nothing
end

function Validation.rules(::Type{<:LineParametersProblem})
    (Validation.OwnerRule(:line_parameters_problem, _check_line_parameters_problem),)
end

"""
$(TYPEDSIGNATURES)

Construct a problem after promoting the system, operating temperature, static
earth model, frequencies, and optional propagation constants to one real scalar
type.

# Keywords

- `temperature`: Operating temperature \\[°C\\].
- `earth_props`: Static earth model.
- `frequencies`: Positive sorted analysis frequencies \\[Hz\\].
- `Γ`: Optional longitudinal propagation constants aligned with `frequencies`
  \\[1/m\\]. A formula that fixes Γ to zero rejects nonzero values.
"""
function LineParametersProblem(
        system::LineCableSystem;
        temperature::Real = oftype(float(system.line_length), 20),
        earth_props::EarthModel,
        frequencies::AbstractVector{<:Real} = [oftype(float(system.line_length), 50)],
        Γ::Union{Nothing, AbstractVector{<:Number}} = nothing
)
    isempty(frequencies) && throw(ArgumentError("frequencies cannot be empty"))
    Γ !== nothing && isempty(Γ) && throw(ArgumentError("Γ cannot be empty"))
    propagation_type = Γ === nothing ? typeof(float(first(frequencies))) :
                       promote_type(
        typeof(float(real(first(Γ)))),
        typeof(float(imag(first(Γ))))
    )
    T = promote_type(
        eltype(system), typeof(float(temperature)), eltype(earth_props),
        typeof(float(first(frequencies))), propagation_type
    )
    converted_system = convert(LineCableSystem{T}, system)
    propagation = Γ === nothing ? nothing :
                  Complex{T}[convert(Complex{T}, value) for value in Γ]
    problem = LineParametersProblem{
        T,
        typeof(converted_system),
        typeof(propagation)
    }(
        converted_system,
        convert(T, float(temperature)),
        convert(EarthModel{T}, earth_props),
        T[convert(T, float(value)) for value in frequencies],
        propagation
    )
    return validate(problem)
end

"""
$(TYPEDSIGNATURES)

Construct a line-parameter problem from completed cable designs and physical
placements.

# Arguments

- `designs`: One completed design, or a collection aligned with `placements`.
- `placements`: One placement or a collection of physical placements.

# Keywords

- `connections`: Terminal-to-phase declarations.
- `environment`: Optional physical environment declaration.
- `system_id`: Stable system identifier.
- `line_length`: Physical line length in metres.
- `temperature`: Operating temperature in °C.
- `earth_props`: Static earth model.
- `frequencies`: Positive sorted analysis frequencies in Hz.
- `Γ`: Optional longitudinal propagation constants aligned with `frequencies`
  in inverse metres.
- `combine`: Rule used to combine designs and placements.

# Returns

One validated [`LineParametersProblem`](@ref) containing a completed
[`LineCableSystem`](@ref).
"""
function LineParametersProblem(
        designs::Union{CableDesign, AbstractVector{<:CableDesign}, Tuple},
        placements,
        connections,
        environment,
        system_id::AbstractString,
        line_length::Real,
        temperature::Real,
        earth_props::EarthModel,
        frequencies::AbstractVector{<:Real};
        Γ::Union{Nothing, AbstractVector{<:Number}} = nothing,
        combine::Symbol = :product
)
    system = build(
        LineCableSystem,
        designs,
        placements;
        connections,
        environment,
        system_id,
        line_length,
        combine
    )
    return LineParametersProblem(system; temperature, earth_props, frequencies, Γ)
end

"""
$(TYPEDEF)

Store the physical methods selected for a line-parameter calculation.

$(TYPEDFIELDS)
"""
struct LineParametersFormulation{M <: NamedTuple, O <: NamedTuple} <: AbstractFormulation
    methods::M
    options::O
end

function LineParametersFormulation(;
        internal_impedance::InternalImpedanceFormulation,
        insulation_impedance::InsulationImpedanceFormulation,
        earth_impedance::EarthImpedanceFormulation,
        insulation_admittance::InsulationAdmittanceFormulation,
        semicon_admittance::SemiconAdmittanceFormulation,
        earth_admittance::EarthAdmittanceFormulation,
        earth_properties,
        equivalent_earth,
        options::NamedTuple
)
    methods = (;
        internal_impedance, insulation_impedance, earth_impedance,
        insulation_admittance, semicon_admittance, earth_admittance, earth_properties,
        equivalent_earth
    )
    return LineParametersFormulation(methods, options)
end

_earth_impedance_formula(formula::EarthImpedanceFormulation) = formula
_earth_impedance_formula(identifier::Symbol) = EarthImpedance.Formula(identifier)
function _earth_impedance_formula(
        selection::FormulaSpec{ID, Order}
) where {ID, Order}
    _direct(selection, :earth_impedance)
    return EarthImpedance.Formula(Val(ID); selection.overrides...)
end

_earth_admittance_formula(formula::EarthAdmittanceFormulation) = formula
_earth_admittance_formula(identifier::Symbol) = EarthAdmittance.Formula(identifier)
function _earth_admittance_formula(
        selection::FormulaSpec{ID, Order}
) where {ID, Order}
    _direct(selection, :earth_admittance)
    return EarthAdmittance.Formula(Val(ID); selection.overrides...)
end

_internal_impedance_formula(formula::InternalImpedanceFormulation) = formula
_internal_impedance_formula(identifier::Symbol) = InternalImpedance.Formula(identifier)
function _internal_impedance_formula(
        selection::FormulaSpec{ID, Order}
) where {ID, Order}
    _direct(selection, :internal_impedance)
    return InternalImpedance.Formula(Val(ID); selection.overrides...)
end

_insulation_impedance_formula(formula::InsulationImpedanceFormulation) = formula
_insulation_impedance_formula(identifier::Symbol) = InsulationImpedance.Formula(identifier)
function _insulation_impedance_formula(
        selection::FormulaSpec{ID, Order}
) where {ID, Order}
    _direct(selection, :insulation_impedance)
    return InsulationImpedance.Formula(Val(ID); selection.overrides...)
end

_insulation_admittance_formula(formula::InsulationAdmittanceFormulation) = formula
function _insulation_admittance_formula(identifier::Symbol)
    InsulationAdmittance.Formula(identifier)
end
function _insulation_admittance_formula(
        selection::FormulaSpec{ID, Order}
) where {ID, Order}
    _direct(selection, :insulation_admittance)
    return InsulationAdmittance.Formula(Val(ID); selection.overrides...)
end

_semicon_admittance_formula(formula::SemiconAdmittanceFormulation) = formula
function _semicon_admittance_formula(identifier::Symbol)
    SemiconAdmittance.Formula(identifier)
end
function _semicon_admittance_formula(
        selection::FormulaSpec{ID, Order}
) where {ID, Order}
    _direct(selection, :semicon_admittance)
    return SemiconAdmittance.Formula(Val(ID); selection.overrides...)
end

_earth_properties_formula(::Nothing) = nothing
_earth_properties_formula(identifier::Symbol) = EarthProps.FD.Formula(identifier)
function _earth_properties_formula(
        selection::FormulaSpec{ID, Order}
) where {ID, Order}
    _direct(selection, :earth_properties)
    return EarthProps.FD.Formula(Val(ID); selection.overrides...)
end
_earth_properties_formula(formula) = formula

_equivalent_earth(::Nothing) = nothing
_equivalent_earth(identifier::Symbol) = EHEM.AfterFD(identifier)
_equivalent_earth(sequence::EHEM.AbstractSequence) = sequence
_equivalent_earth(rule::EHEM.AbstractRule) = EHEM.AfterFD(rule)

function _direct(::FormulaSpec{ID, Order}, owner::Symbol) where {ID, Order}
    Order === :default || throw(ArgumentError(
        "formula order is only valid for equivalent_earth; got :$Order for $owner"
    ))
    return nothing
end

function _ehem_rule(selection::FormulaSpec{:Layer})
    required = (:layer,)
    unknown = setdiff(keys(selection.overrides), required)
    isempty(unknown) || throw(ArgumentError(
        "unknown assumptions for equivalent-earth policy :Layer: $(collect(unknown))"
    ))
    haskey(selection.overrides, :layer) || throw(ArgumentError(
        "equivalent-earth policy :Layer requires `layer`"
    ))
    layer = selection.overrides.layer
    layer isa Int || throw(ArgumentError(
        "equivalent-earth policy :Layer requires an integer `layer`"
    ))
    return EHEM.Layer(layer)
end

function _ehem_rule(selection::FormulaSpec{ID}) where {ID}
    return EHEM.Formula(Val(ID); selection.overrides...)
end

_ehem_order(::Val{:default}, rule) = EHEM.AfterFD(rule)
_ehem_order(::Val{:after}, rule) = EHEM.AfterFD(rule)
_ehem_order(::Val{:before}, rule) = EHEM.BeforeFD(rule)

function _equivalent_earth(
        selection::FormulaSpec{ID, Order}
) where {ID, Order}
    return _ehem_order(Val(Order), _ehem_rule(selection))
end

function Formulation(;
        internal_impedance = formula(:default),
        insulation_impedance = formula(:default),
        earth_impedance = formula(:default),
        insulation_admittance = formula(:default),
        semicon_admittance = formula(:default),
        earth_admittance = formula(:default),
        earth_properties = formula(:default),
        equivalent_earth = formula(:default),
        options::NamedTuple = (;)
)
    return LineParametersFormulation(;
        internal_impedance = _internal_impedance_formula(internal_impedance),
        insulation_impedance = _insulation_impedance_formula(insulation_impedance),
        earth_impedance = _earth_impedance_formula(earth_impedance),
        insulation_admittance = _insulation_admittance_formula(insulation_admittance),
        semicon_admittance = _semicon_admittance_formula(semicon_admittance),
        earth_admittance = _earth_admittance_formula(earth_admittance),
        earth_properties = _earth_properties_formula(earth_properties),
        equivalent_earth = _equivalent_earth(equivalent_earth),
        options = formulation_options(Val(LineParametersFormulation), options)
    )
end
