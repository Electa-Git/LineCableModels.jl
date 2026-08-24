"""
$(TYPEDEF)

Place one materialised cable design in a line cross-section.

$(TYPEDFIELDS)
"""
struct CablePosition{T <: Real}
    design_data::CableDesign{T}
    horz::T
    vert::T
    conn::Vector{Int}
end

Base.eltype(::CablePosition{T}) where {T} = T
Base.eltype(::Type{CablePosition{T}}) where {T} = T

function _check_cable_position(position::CablePosition)
    length(position.conn) == length(position.design_data.components) ||
        throw(DimensionMismatch("phase mapping must match the component count"))
    all(>=(0), position.conn) ||
        throw(DomainError(position.conn, "phase assignments must be nonnegative"))
    iszero(position.vert) && throw(DomainError(
        position.vert,
        "a cable centre cannot lie on the air/earth interface"
    ))
    radius = outer_radius(position.design_data)
    abs(position.vert) >= radius || throw(DomainError(
        position.vert,
        "the cable cross-section crosses the air/earth interface"
    ))
    return nothing
end

function Validation.rules(::Type{<:CablePosition})
    (Validation.OwnerRule(:cable_position_geometry, _check_cable_position),)
end

function CablePosition(
        design::CableDesign,
        horizontal::Real,
        vertical::Real,
        connections::Union{Nothing, AbstractDict} = nothing
)
    names = getproperty.(design.components, :id)
    mapping = connections === nothing ?
              [index == 1 ? 1 : 0 for index in eachindex(names)] :
    begin
        unknown = setdiff(String.(keys(connections)), names)
        isempty(unknown) || throw(KeyError(first(unknown)))
        [Int(get(connections, name, 0)) for name in names]
    end
    T = promote_type(
        eltype(design), typeof(float(horizontal)), typeof(float(vertical))
    )
    return validate(CablePosition{T}(
        convert(CableDesign{T}, design),
        convert(T, float(horizontal)),
        convert(T, float(vertical)),
        mapping
    ))
end

function Base.convert(::Type{CablePosition{T}}, position::CablePosition) where {T <: Real}
    return validate(CablePosition{T}(
        convert(CableDesign{T}, position.design_data),
        convert(T, position.horz),
        convert(T, position.vert),
        copy(position.conn)
    ))
end

function Base.convert(::Type{CablePosition{T}}, position::CablePosition{T}) where {T <:
                                                                                   Real}
    position
end

"""
$(TYPEDEF)

Store positioned cables and the materialised line length.

$(TYPEDFIELDS)
"""
mutable struct LineCableSystem{T <: Real}
    system_id::String
    line_length::T
    cables::Vector{CablePosition{T}}
end

Base.eltype(::LineCableSystem{T}) where {T} = T
Base.eltype(::Type{LineCableSystem{T}}) where {T} = T

ncables(system::LineCableSystem) = length(system.cables)
function nphases(system::LineCableSystem)
    assignments = unique(collect(Iterators.flatten(
        position.conn for position in system.cables
    )))
    return count(>(0), assignments)
end

function _check_line_cable_system(system::LineCableSystem)
    isempty(system.system_id) && throw(ArgumentError("system_id cannot be empty"))
    system.line_length > zero(system.line_length) ||
        throw(DomainError(system.line_length, "line length must be positive"))
    isempty(system.cables) && throw(ArgumentError("a line system requires one cable"))
    foreach(validate, system.cables)
    for left_index in eachindex(system.cables)
        for right_index in (left_index + 1):lastindex(system.cables)
            overlaps(system.cables[left_index], system.cables[right_index]) &&
                throw(DomainError(
                    (left_index, right_index),
                    "cable cross-sections overlap"
                ))
        end
    end
    return nothing
end

function Validation.rules(::Type{<:LineCableSystem})
    (Validation.OwnerRule(:line_cable_system_structure, _check_line_cable_system),)
end

function LineCableSystem(
        system_id::AbstractString,
        line_length::Real,
        position::CablePosition
)
    T = promote_type(typeof(float(line_length)), eltype(position))
    return validate(LineCableSystem{T}(
        String(system_id),
        convert(T, float(line_length)),
        CablePosition{T}[convert(CablePosition{T}, position)]
    ))
end

function LineCableSystem(
        system_id::AbstractString,
        line_length::Real,
        positions::AbstractVector{<:CablePosition}
)
    isempty(positions) && throw(ArgumentError("a line system requires one cable"))
    T = promote_type(typeof(float(line_length)), eltype.(positions)...)
    return validate(LineCableSystem{T}(
        String(system_id),
        convert(T, float(line_length)),
        CablePosition{T}[convert(CablePosition{T}, item) for item in positions]
    ))
end

function LineCableSystem(
        system_id::AbstractString,
        line_length::Real,
        design::CableDesign,
        horizontal::Real,
        vertical::Real,
        connections::Union{Nothing, AbstractDict} = nothing
)
    return LineCableSystem(
        system_id,
        line_length,
        CablePosition(design, horizontal, vertical, connections)
    )
end

function add!(system::LineCableSystem{T}, position::CablePosition{T}) where {T}
    foreach(system.cables) do existing
        overlaps(existing, position) && throw(DomainError(
            (existing, position),
            "cable cross-sections overlap"
        ))
    end
    candidate = validate(LineCableSystem{T}(
        system.system_id,
        system.line_length,
        CablePosition{T}[system.cables; position]
    ))
    system.cables = candidate.cables
    return system
end

function add!(system::LineCableSystem{T}, position::CablePosition{U}) where {T, U}
    throw(ArgumentError(
        "cannot add CablePosition{$U} to LineCableSystem{$T}; explicitly convert " *
        "the complete system before mutation",
    ))
end

function add!(
        system::LineCableSystem,
        design::CableDesign,
        horizontal::Real,
        vertical::Real,
        connections::Union{Nothing, AbstractDict} = nothing
)
    return add!(system, CablePosition(design, horizontal, vertical, connections))
end

function Base.convert(::Type{LineCableSystem{T}}, system::LineCableSystem) where {T <: Real}
    return validate(LineCableSystem{T}(
        system.system_id,
        convert(T, system.line_length),
        CablePosition{T}[convert(CablePosition{T}, position) for position in system.cables]
    ))
end

function Base.convert(::Type{LineCableSystem{T}}, system::LineCableSystem{T}) where {T <:
                                                                                     Real}
    system
end

include("dataframe.jl")

function Base.show(io::IO, ::MIME"text/plain", system::LineCableSystem)
    println(
        io,
        "LineCableSystem \"$(system.system_id)\": [line_length=$(system.line_length), " *
        "num_cables=$(ncables(system)), num_phases=$(nphases(system))]"
    )
    println(io, "└─ $(ncables(system))-element CablePosition:")
    for (index, position) in enumerate(system.cables)
        prefix = index == ncables(system) ? "   └─" : "   ├─"
        components = getproperty.(position.design_data.components, :id)
        mapping = join(
            ["$component→$phase" for (component, phase) in zip(components, position.conn)],
            ", "
        )
        println(
            io,
            "$prefix CableDesign \"$(position.design_data.cable_id)\": " *
            "[horz=$(round(position.horz, sigdigits=4)), " *
            "vert=$(round(position.vert, sigdigits=4)), conn=($mapping)]"
        )
    end
end
