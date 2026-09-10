"""
$(TYPEDEF)

Define one scalar line-parameter calculation over a completed cable system.
Operating temperature and analysis frequencies are fields of the problem.

$(TYPEDFIELDS)
"""
struct LineParametersProblem{
    T <: Real,
    S <: LineCableSystem{T},
    P <: Union{Nothing, Vector{Complex{T}}},
    E <: EarthModel{T}
} <: AbstractProblemDefinition
    "Physical cable system."
    system::S
    "Operating temperature \\[°C\\]."
    temperature::T
    "Static earth model."
    earth_props::E
    "Strictly positive, sorted analysis frequencies \\[Hz\\]."
    frequencies::Vector{T}
    "Optional longitudinal propagation constants aligned with frequency \\[1/m\\]."
    Γ::P

    function LineParametersProblem{T, S, P, E}(
            system::S,
            temperature::T,
            earth_props::E,
            frequencies::Vector{T},
            Γ::P
    ) where {
            T <: Real,
            S <: LineCableSystem{T},
            P <: Union{Nothing, Vector{Complex{T}}},
            E <: EarthModel{T}
    }
        return validate(new{T, S, P, E}(
            system,
            temperature,
            earth_props,
            frequencies,
            Γ
        ))
    end
end

Base.eltype(::LineParametersProblem{T}) where {T} = T
Base.eltype(::Type{LineParametersProblem{T}}) where {T} = T

function validate(problem::LineParametersProblem)
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
    isfinite(problem.temperature) || throw(DomainError(problem.temperature,
        "operating temperature must be finite"))
    isempty(problem.frequencies) && throw(ArgumentError("frequencies cannot be empty"))
    all(value -> isfinite(value) && value > zero(value), problem.frequencies) ||
        throw(DomainError(
            problem.frequencies, "frequencies must be positive and finite"
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
    return problem
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
    converted_earth = convert(EarthModel{T}, earth_props)
    propagation = Γ === nothing ? nothing :
                  Complex{T}[convert(Complex{T}, value) for value in Γ]
    return LineParametersProblem{
        T,
        typeof(converted_system),
        typeof(propagation),
        typeof(converted_earth)
    }(
        converted_system,
        convert(T, float(temperature)),
        converted_earth,
        T[convert(T, float(value)) for value in frequencies],
        propagation
    )
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
struct LineParametersFormulation{M <: NamedTuple, O <: NamedTuple, D <: NamedTuple} <:
       AbstractFormulation
    "Owner-resolved physical methods; context-dependent defaults remain deferred."
    methods::M
    "Shared physical computation options."
    options::O
    "Requested selections retained before owner and problem-context resolution."
    definitions::D
end

function LineParametersFormulation(methods::NamedTuple, options::NamedTuple)
    LineParametersFormulation(methods, options, methods)
end

function LineParametersFormulation(;
        internal_impedance::InternalImpedanceFormulation,
        insulation_impedance::InsulationImpedanceFormulation,
        earth_impedance::Union{EarthImpedanceFormulation, NamedTuple},
        insulation_admittance::InsulationAdmittanceFormulation,
        semicon_admittance::SemiconAdmittanceFormulation,
        earth_admittance::Union{EarthAdmittanceFormulation, NamedTuple},
        earth_properties,
        pipe_impedance::PipeImpedanceFormulation,
        temperature_dependence::Union{Nothing, TemperatureDependent.Formula} = TemperatureDependent.Formula(:default),
        options::NamedTuple
)
    methods = (;
        internal_impedance, insulation_impedance, earth_impedance,
        insulation_admittance, semicon_admittance, earth_admittance, earth_properties,
        pipe_impedance, temperature_dependence
    )
    return LineParametersFormulation(methods, options)
end

function _line_formulation(
        internal_impedance,
        insulation_impedance,
        earth_impedance,
        insulation_admittance,
        semicon_admittance,
        earth_admittance,
        earth_properties,
        pipe_impedance,
        temperature_dependence,
        options::NamedTuple
)
    selected = LineParametersFormulation(;
        internal_impedance = InternalImpedance.Formula(internal_impedance),
        insulation_impedance = InsulationImpedance.Formula(insulation_impedance),
        earth_impedance = Formulation(EarthImpedance.Formula, earth_impedance),
        insulation_admittance = InsulationAdmittance.Formula(insulation_admittance),
        semicon_admittance = SemiconAdmittance.Formula(semicon_admittance),
        earth_admittance = Formulation(EarthAdmittance.Formula, earth_admittance),
        earth_properties = earth_properties === nothing ? nothing :
                           Earth.FrequencyDependent.Formula(earth_properties),
        pipe_impedance = PipeImpedance.Formula(pipe_impedance),
        temperature_dependence = temperature_dependence === nothing ? nothing :
                                 TemperatureDependent.Formula(temperature_dependence),
        options = formulation_options(LineParametersFormulation, options)
    )
    definitions = (; internal_impedance, insulation_impedance,
        earth_impedance = earth_impedance isa NamedTuple ?
            NamedTuple{keys(selected.methods.earth_impedance)}(earth_impedance) : earth_impedance,
        insulation_admittance, semicon_admittance,
        earth_admittance = earth_admittance isa NamedTuple ?
            NamedTuple{keys(selected.methods.earth_admittance)}(earth_admittance) : earth_admittance,
        earth_properties, pipe_impedance, temperature_dependence)
    return LineParametersFormulation(selected.methods, selected.options, definitions)
end

"""
$(TYPEDSIGNATURES)

Select the complete physical-method bundle for a line-parameter calculation.

`earth_impedance` and `earth_admittance` each accept one formula or a NamedTuple
with exactly `air`, `earth`, and `mixed` selections. For a physical horizontal
air/soil two-half-space model, these select `(s,t)=(1,1)`, `(2,2)`, and the two
cross-layer mutual directions. Each required kind/layer case is validated against
the selected equation. Missing cases have no implicit fallback. The shorthand is
rejected for layered soil; scalar multilayer and explicit equivalent-earth
selections retain their own contracts. Formula hooks and numerical options remain
local to each selected entry.

`temperature_dependence=formula(:default)` selects the Materials-owned linear
resistivity law. `nothing` retains reference resistivity. Operating temperature
belongs to the problem; reference temperature and coefficients belong to each
material.

Each formula slot and the complete `options` tuple accepts either one scalar
selection or an explicit
[`Grid`](@ref LineCableModels.ParametricBuilder.Grid)/
[`Gridspace`](@ref LineCableModels.ParametricBuilder.Gridspace) source. Scalar
inputs return one [`LineParametersFormulation`](@ref). Varying inputs return a
`Gridspace{LineParametersFormulation}` whose points contain only completed,
owner-resolved formula values. Context-dependent `:default` selections remain
deferred until a scalar problem supplies its geometry and earth context.

`combine=:product` forms the Cartesian product among varying fields in this
formulation. `combine=:zip` aligns equally sized fields and broadcasts
singletons. This composition is independent of the Cartesian product between
problem points and formulation points performed by
[`Combinatorial`](@ref LineCableModels.ParametricBuilder.Combinatorial).
"""
function Formulation(;
        internal_impedance = formula(:default),
        insulation_impedance = formula(:default),
        earth_impedance = formula(:default),
        insulation_admittance = formula(:default),
        semicon_admittance = formula(:default),
        earth_admittance = formula(:default),
        earth_properties = formula(:default),
        pipe_impedance = formula(:default),
        temperature_dependence = formula(:default),
        options = (;),
        combine::Symbol = :product
)
    values = (
        internal_impedance,
        insulation_impedance,
        earth_impedance,
        insulation_admittance,
        semicon_admittance,
        earth_admittance,
        earth_properties,
        pipe_impedance,
        temperature_dependence,
        options
    )
    return parameterize(
        LineParametersFormulation,
        _line_formulation,
        values;
        combine
    )
end

"""
$(TYPEDSIGNATURES)

Expose complete requested and resolved formula selections and reduction options.
Formula hooks retain their concrete callables; serialization belongs to the writer.
"""
function Base.NamedTuple(value::LineParametersFormulation)
    record = function (selected)
        selected === nothing && return nothing
        selected isa Symbol && return NamedTuple(formula(selected))
        selected isa NamedTuple && return map(record,selected)
        return NamedTuple(selected)
    end
    # Retain complete records without specializing result containers on each
    # parameter tuple or hook type. The stored values remain unchanged.
    Record=NamedTuple{(:backend,:requested,:methods,:options),
        Tuple{Symbol,NamedTuple,NamedTuple,NamedTuple}}
    return Record((:coaxial,map(record,value.definitions),map(record,value.methods),value.options))
end
