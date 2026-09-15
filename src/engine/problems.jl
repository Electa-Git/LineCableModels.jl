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
    DataModel.clearance_geometry(problem.system.designs, problem.system.positions;
        required = problem.system.clearances, interface = true, adjust = false)
    phases = unique(problem.system.connection_order)
    active = filter(>(0), phases)
    isempty(active) && throw(ArgumentError(
        "at least one conductor must be assigned to an active phase",
    ))
    maximum(active) <= nphases(problem.system) || throw(DomainError(
        active,
        "an active-phase assignment exceeds the number of distinct active phases"
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
    converted_system = DataModel.interface_clearance(convert(LineCableSystem{T}, system))
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

- `connections`: Terminal-to-active-phase declarations; use one-based active
  phase IDs and `0` for grounded/eliminated conductors.
- `environment`: Optional physical environment declaration.
- `system_id`: Stable system identifier.
- `line_length`: Physical line length in meters.
- `temperature`: Operating temperature in °C.
- `earth_props`: Static earth model.
- `frequencies`: Positive sorted analysis frequencies in Hz.
- `Γ`: Optional longitudinal propagation constants aligned with `frequencies`
  in inverse meters.
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

"""Identify the owned coaxial calculation without report numbering."""
description(::Type{<:LineParametersFormulation}; compact::Bool=false) = description(LineCableModelsCoaxial;compact)
description(::LineParametersFormulation; compact::Bool=false) = description(LineParametersFormulation;compact)
formula_id(::Type{<:LineParametersFormulation}) = :coaxial
formula_id(::LineParametersFormulation) = :coaxial

"""
$(TYPEDSIGNATURES)

Iterate ordered formula-slot owners relevant to a physical quantity. With
`quantity=nothing`, retain every declared slot. This is a declaration read,
not problem-dependent formula resolution.
"""
function Base.pairs(::Type{LineParametersFormulation}; quantity=nothing)
    selected=(internal_impedance=InternalImpedance.Formula,
        insulation_impedance=InsulationImpedance.Formula,
        earth_impedance=EarthImpedance.Formula,
        shunt_model=ShuntModel.Formula,
        insulation_admittance=InsulationAdmittance.Formula,
        semicon_admittance=SemiconAdmittance.Formula,
        earth_admittance=EarthAdmittance.Formula,
        earth_properties=Earth.FrequencyDependent.Formula,
        pipe_impedance=PipeImpedance.Formula,
        temperature_dependence=TemperatureDependent.Formula)
    quantity===nothing && return pairs(selected)
    q=quantity isa Units.Quantity ? quantity : Grammar.request_quantity(quantity)
    series=q in (Units.quantity(Z),Units.quantity(R),Units.quantity(X),Units.quantity(L),
        Units.quantity(Z,abs),Units.quantity(Z,angle))
    shunt=q in (Units.quantity(Y),Units.quantity(G),Units.quantity(B),Units.quantity(C),
        Units.quantity(Y,abs),Units.quantity(Y,angle))
    series || shunt || throw(ArgumentError("no line-formulation selections for $q"))
    omitted=series ? (:shunt_model,:insulation_admittance,:semicon_admittance,:earth_admittance) :
        (:internal_impedance,:insulation_impedance,:earth_impedance,:pipe_impedance)
    return pairs((; (key=>value for (key,value) in pairs(selected) if key ∉ omitted)...))
end

description(::Type{LineParametersFormulation},::Val{:internal_impedance}) = "internal Z"
description(::Type{LineParametersFormulation},::Val{:insulation_impedance}) = "insulation Z"
description(::Type{LineParametersFormulation},::Val{:earth_impedance}) = "earth Z"
description(::Type{LineParametersFormulation},::Val{:shunt_model}) = "shunt geometry"
description(::Type{LineParametersFormulation},::Val{:insulation_admittance}) = "insulation Y"
description(::Type{LineParametersFormulation},::Val{:semicon_admittance}) = "semicon Y"
description(::Type{LineParametersFormulation},::Val{:earth_admittance}) = "earth Y"
description(::Type{LineParametersFormulation},::Val{:earth_properties}) = "soil law"
description(::Type{LineParametersFormulation},::Val{:pipe_impedance}) = "pipe Z"
description(::Type{LineParametersFormulation},::Val{:temperature_dependence}) = "temperature law"

"""Expose typed children and controls without serializing the formulation."""
Base.pairs(value::LineParametersFormulation; quantity=nothing) =
    pairs(LineParametersFormulation,(methods=value.methods,requested=map(formulation_options,value.definitions),options=value.options);quantity)
formulation_options(value::LineParametersFormulation) = value.options
formulation_options(::Type{LineParametersFormulation},retained::NamedTuple,::Val{:retained}) = retained.options

"""
$(TYPEDSIGNATURES)

Iterate owner-scoped selections from a formulation's declared child interface.
`retained.methods` contains typed children (including owner-bound saved leaves);
`retained.requested` contains their explicit controls. `pairs(owner; quantity)`
owns child order and relevance. The empty route identifies the owner itself.
"""
function Base.pairs(::Type{LineParametersFormulation}, retained::NamedTuple;
        quantity=nothing,owner=LineParametersFormulation)
    entries=Pair{Tuple,Any}[(owner,()) => (owner => (options=retained.options,))]
    for (slot,family) in pairs(owner;quantity)
        selected=retained.methods[slot]
        requested=retained.requested[slot]
        if selected isa NamedTuple
            children = pairs(family)
            Set(keys(selected)) == Set(key for (key,_) in children) || throw(ArgumentError(
                "retained $slot selections do not match the owning formula slots"))
            for (route,_) in children
                value=selected[route]
                controls=requested[route]
                push!(entries,(owner,(slot,route)) => (value===nothing || ismissing(value) ? value : value => controls))
            end
        else
            controls=requested
            push!(entries,(owner,(slot,)) => (selected===nothing || ismissing(selected) ? selected : selected => controls))
        end
    end
    return entries
end

function LineParametersFormulation(methods::NamedTuple, options::NamedTuple)
    LineParametersFormulation(methods, options, methods)
end

function LineParametersFormulation(;
        internal_impedance::Union{InternalImpedanceFormulation, NamedTuple},
        insulation_impedance::InsulationImpedanceFormulation,
        earth_impedance::Union{EarthImpedanceFormulation, NamedTuple},
        insulation_admittance::InsulationAdmittanceFormulation,
        semicon_admittance::SemiconAdmittanceFormulation,
        earth_admittance::Union{EarthAdmittanceFormulation, NamedTuple},
        shunt_model::ShuntModelFormulation = ShuntModel.Formula(:default),
        earth_properties,
        pipe_impedance::PipeImpedanceFormulation,
        temperature_dependence::Union{Nothing, TemperatureDependent.Formula} = TemperatureDependent.Formula(:default),
        options::NamedTuple
)
    methods = (;
        internal_impedance, insulation_impedance, earth_impedance, shunt_model,
        insulation_admittance, semicon_admittance, earth_admittance, earth_properties,
        pipe_impedance, temperature_dependence
    )
    return LineParametersFormulation(methods, options)
end

function _line_formulation(
        internal_impedance,
        insulation_impedance,
        earth_impedance,
        shunt_model,
        insulation_admittance,
        semicon_admittance,
        earth_admittance,
        earth_properties,
        pipe_impedance,
        temperature_dependence,
        options::NamedTuple
)
    selected = LineParametersFormulation(;
        internal_impedance = Formulation(InternalImpedance.Formula, internal_impedance),
        insulation_impedance = InsulationImpedance.Formula(insulation_impedance),
        earth_impedance = Formulation(EarthImpedance.Formula, earth_impedance),
        shunt_model = ShuntModel.Formula(shunt_model),
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
    definitions = (; internal_impedance = internal_impedance isa NamedTuple ?
            NamedTuple{keys(selected.methods.internal_impedance)}(internal_impedance) : internal_impedance,
        insulation_impedance,
        earth_impedance = earth_impedance isa NamedTuple ?
            NamedTuple{keys(selected.methods.earth_impedance)}(earth_impedance) : earth_impedance,
        shunt_model, insulation_admittance, semicon_admittance,
        earth_admittance = earth_admittance isa NamedTuple ?
            NamedTuple{keys(selected.methods.earth_admittance)}(earth_admittance) : earth_admittance,
        earth_properties, pipe_impedance, temperature_dependence)
    return LineParametersFormulation(selected.methods, selected.options, definitions)
end

"""
$(TYPEDSIGNATURES)

Select the complete physical-method bundle for a line-parameter calculation.

`shunt_model=:default` (or `:coaxial`) selects annular local shunt geometry.
`:boundary` explicitly prepares a lossless wire/tape boundary correction.
This choice is independent of `insulation_admittance` and `semicon_admittance`,
which select material constitutive laws. Boundary numerical controls and an
explicit fallback belong to `formula(:boundary; options, parameters)`.

`internal_impedance` accepts one formula or complete `inner`, `outer`, and
`transfer` selections. Each selected formula owns its callable overrides and
numerical controls. Only the surface kinds required by the geometry are
evaluated; unsupported or explicitly unused overrides fail during preflight.

`earth_impedance` and `earth_admittance` each accept one formula or a NamedTuple
with exactly `air`, `earth`, and `mixed` selections. For a physical horizontal
air/soil two-half-space model, these select `(s,t)=(1,1)`, `(2,2)`, and the two
cross-layer mutual directions. Each required kind/layer case is validated against
the selected equation. Missing cases have no implicit fallback. The shorthand is
rejected for layered soil; scalar multilayer and explicit equivalent-earth
selections retain their own requirements. Formula hooks and numerical options remain
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
        shunt_model = formula(:default),
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
        shunt_model,
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
