"""
$(TYPEDEF)

Define one materialised line-parameter calculation. Operating temperature and
analysis frequencies are fields of the problem.

$(TYPEDFIELDS)
"""
struct LineParametersProblem{T <: Real} <: AbstractProblemDefinition
    "Physical cable system."
    system::LineCableSystem{T}
    "Operating temperature \\[°C\\]."
    temperature::T
    "Static earth model."
    earth_props::EarthModel{T}
    "Strictly positive, sorted analysis frequencies \\[Hz\\]."
    frequencies::Vector{T}
end

Base.eltype(::LineParametersProblem{T}) where {T} = T
Base.eltype(::Type{LineParametersProblem{T}}) where {T} = T

function _check_line_parameters_problem(problem::LineParametersProblem)
    validate(problem.system)
    validate(problem.earth_props)
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
        # Temperature is formulation-independent, but the current linear
        # correction range uses a reference material only when the radial
        # analytical adapter is available. Unsupported physical geometry is
        # rejected later by the backend that requests `AnalyticalInput`.
        design.effective === nothing && continue
        reference = first(design.effective).conductor.material.T0
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
earth model, and frequencies to one real scalar type.
"""
function LineParametersProblem(
        system::LineCableSystem;
        temperature::Real = oftype(float(system.line_length), 20),
        earth_props::EarthModel,
        frequencies::AbstractVector{<:Real} = [oftype(float(system.line_length), 50)]
)
    isempty(frequencies) && throw(ArgumentError("frequencies cannot be empty"))
    T = promote_type(
        eltype(system), typeof(float(temperature)), eltype(earth_props),
        typeof(float(first(frequencies)))
    )
    problem = LineParametersProblem{T}(
        convert(LineCableSystem{T}, system),
        convert(T, float(temperature)),
        convert(EarthModel{T}, earth_props),
        T[convert(T, float(value)) for value in frequencies]
    )
    return validate(problem)
end

"""
$(TYPEDEF)

Define an explicit cable-constant calculation for one materialised design.

$(TYPEDFIELDS)
"""
struct CableConstantsProblem{D <: CableDesign, S, R <: Real} <: AbstractProblemDefinition
    design::D
    separation::S
    earth_resistivity::R

    function CableConstantsProblem(
            design::D,
            separation::S,
            earth_resistivity::R
    ) where {D <: CableDesign, S <: Union{Nothing, Real}, R <: Real}
        separation === nothing || separation > zero(separation) ||
            throw(DomainError(separation, "separation must be positive"))
        earth_resistivity > zero(earth_resistivity) || throw(DomainError(
            earth_resistivity, "earth_resistivity must be positive"
        ))
        return new{D, S, R}(
            design, separation, earth_resistivity
        )
    end
end

function CableConstantsProblem(
        design::CableDesign;
        separation = nothing,
        earth_resistivity = oftype(
            float(design.effective === nothing ? 100.0 :
                  first(design.effective).conductor.material.rho),
            100
        )
)
    cable_constants_problem(design, separation, earth_resistivity)
end

"""
$(TYPEDSIGNATURES)

Construct a cable-constant problem from resolved inputs.

ParametricBuilder extends this owner-visible staging boundary for explicit
finite sources.
"""
function cable_constants_problem(
        design::CableDesign,
        separation::Union{Nothing, Real},
        earth_resistivity::Real
)
    CableConstantsProblem(design, separation, earth_resistivity)
end

"""
$(TYPEDEF)

Store the concrete physical methods selected for an analytical calculation.

$(TYPEDFIELDS)
"""
struct AnalyticalFormulation{B, M <: NamedTuple, O <: NamedTuple} <: AbstractFormulation
    backend::B
    methods::M
    options::O
end

function AnalyticalFormulation(;
        internal_impedance::InternalImpedanceFormulation,
        insulation_impedance::InsulationImpedanceFormulation,
        earth_impedance::EarthImpedanceFormulation,
        insulation_admittance::InsulationAdmittanceFormulation,
        earth_admittance::EarthAdmittanceFormulation,
        earth_properties,
        modal_transform::Union{AbstractTransformFormulation, Nothing},
        equivalent_earth::Union{AbstractEHEMFormulation, Nothing},
        options::NamedTuple
)
    methods = (;
        internal_impedance, insulation_impedance, earth_impedance,
        insulation_admittance, earth_admittance, earth_properties,
        modal_transform, equivalent_earth
    )
    return AnalyticalFormulation(Val(:analytical), methods, options)
end

function Base.getproperty(formulation::AnalyticalFormulation, name::Symbol)
    name in fieldnames(typeof(formulation)) && return getfield(formulation, name)
    methods = getfield(formulation, :methods)
    haskey(methods, name) && return getproperty(methods, name)
    return getfield(formulation, name)
end

function Formulation(
        ::Val{:analytical};
        internal_impedance::InternalImpedanceFormulation = InternalImpedance.ScaledBessel(),
        insulation_impedance::InsulationImpedanceFormulation = InsulationImpedance.Lossless(),
        earth_impedance::EarthImpedanceFormulation = EarthImpedance.Papadopoulos(),
        insulation_admittance::InsulationAdmittanceFormulation = InsulationAdmittance.Lossless(),
        earth_admittance::EarthAdmittanceFormulation = EarthAdmittance.Papadopoulos(),
        earth_properties = EarthProperties.CPEarth(),
        modal_transform::Union{AbstractTransformFormulation, Nothing} = nothing,
        equivalent_earth::Union{AbstractEHEMFormulation, Nothing} = nothing,
        options::NamedTuple = (;)
)
    return AnalyticalFormulation(;
        internal_impedance, insulation_impedance, earth_impedance,
        insulation_admittance, earth_admittance, earth_properties,
        modal_transform, equivalent_earth,
        options = formulation_options(Val(AnalyticalFormulation), options)
    )
end
