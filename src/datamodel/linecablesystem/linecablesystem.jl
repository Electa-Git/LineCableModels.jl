"""
$(TYPEDEF)

Store authoritative cable placements and their eager global terminal state.

`designs`, `positions`, `connections`, and `environment` are declarations.
Global geometry, terminal order, terminal indices, and the flattened
connection order are derived by the constructor.

$(TYPEDFIELDS)
"""
struct LineCableSystem{
        T <: Real,
        D <: AbstractVector,
        P <: AbstractVector{<:Pose2{T}},
        C <: AbstractVector,
        E,
        G <: AbstractVector
}
    "Stable system identifier."
    system_id::String
    "Physical line length [m]."
    line_length::T
    "Eager cable designs in placement order."
    designs::D
    "Cable poses in the system frame."
    positions::P
    "Per-cable terminal connection declarations in terminal order."
    connections::C
    "Optional physical environment declaration."
    environment::E
    "Resolved cable regions placed in the system frame."
    geometry::G
    "Global `(cable, terminal)` records in deterministic order."
    terminal_order::Vector{NamedTuple{(:cable, :terminal), Tuple{Int, Symbol}}}
    "Global terminal index for every resolved system region."
    terminal_map::Vector{Int}
    "Connection assignments aligned with `terminal_order`."
    connection_order::Vector{Int}
end

Base.eltype(::LineCableSystem{T}) where {T} = T
Base.eltype(::Type{<:LineCableSystem{T}}) where {T} = T

ncables(system::LineCableSystem) = length(system.designs)
nphases(system::LineCableSystem) =
    length(unique(filter(>(0), system.connection_order)))

function LineCableSystem(
        designs::AbstractVector{<:CableDesign};
        positions,
        connections = nothing,
        orientations = nothing,
        environment = nothing,
        system_id::AbstractString = "line-cable-system",
        line_length::Real = 1
)
    # 1. Place every eager design in the system frame. Position tuples are a
    # constructor shorthand only; the stored declaration is always Pose2.
    isempty(designs) && throw(ArgumentError("a line system requires one cable"))
    identifier = String(system_id)
    isempty(identifier) && throw(ArgumentError("system_id cannot be empty"))
    isfinite(line_length) && line_length > zero(line_length) || throw(DomainError(
        line_length,
        "line length must be positive"
    ))
    position_values = if positions isa Pose2
        Pose2[positions]
    elseif positions isa Tuple && length(designs) == 1 &&
            length(positions) in (2, 3) && all(value -> value isa Real, positions)
        Pose2[Pose2(
            positions[1],
            positions[2],
            length(positions) == 3 ? positions[3] : 0
        )]
    elseif positions isa AbstractVector || positions isa Tuple
        Pose2[
            position isa Pose2 ? position :
            position isa Tuple && length(position) in (2, 3) &&
            all(value -> value isa Real, position) ?
            Pose2(
                position[1],
                position[2],
                length(position) == 3 ? position[3] : 0
            ) : throw(ArgumentError(
                "positions must contain Pose2 values or `(x, y[, φ])` tuples"
            ))
            for position in positions
        ]
    else
        throw(ArgumentError(
            "positions must be a Pose2, coordinate tuple, or position collection"
        ))
    end
    length(position_values) == length(designs) || throw(DimensionMismatch(
        "position count must match the cable-design count"
    ))
    if orientations !== nothing
        values = orientations isa Real && length(designs) == 1 ?
                 (orientations,) : orientations
        (values isa AbstractVector || values isa Tuple) || throw(ArgumentError(
            "orientations must be a real value or one real value per cable"
        ))
        length(values) == length(designs) || throw(DimensionMismatch(
            "orientation count must match the cable-design count"
        ))
        oriented = Pose2[]
        for (position, orientation) in zip(position_values, values)
            orientation isa Real && isfinite(orientation) || throw(DomainError(
                orientation,
                "cable orientations must be finite real values"
            ))
            push!(oriented, Pose2(position.x, position.y, orientation))
        end
        position_values = oriented
    end

    T = promote_type(
        typeof(float(line_length)),
        (eltype(design) for design in designs)...,
        (eltype(position) for position in position_values)...
    )
    poses = Pose2{T}[convert(Pose2{T}, position) for position in position_values]
    declared_designs = CableDesign[designs...]

    # 2. Resolve global primitive and terminal order while placing the already
    # resolved cable geometry. Cable order and each design's terminal order are
    # retained verbatim.
    terminal_order = NamedTuple{(:cable, :terminal), Tuple{Int, Symbol}}[]
    terminal_offsets = Int[]
    offset = 0
    for (cable_index, design) in enumerate(declared_designs)
        push!(terminal_offsets, offset)
        for terminal in design.terminal_order
            push!(terminal_order, (cable = cable_index, terminal = terminal))
            offset += 1
        end
    end

    # 3. Resolve connection definitions independently for every design.
    declarations = if connections === nothing
        fill(nothing, length(designs))
    elseif length(designs) == 1 &&
            (connections isa AbstractDict || connections isa NamedTuple ||
             connections isa AbstractVector{<:Integer})
        Any[connections]
    elseif connections isa AbstractVector || connections isa Tuple
        length(connections) == length(designs) || throw(DimensionMismatch(
            "connection declaration count must match the cable-design count"
        ))
        collect(connections)
    else
        throw(ArgumentError(
            "connections must be one declaration per cable design"
        ))
    end
    normalized_connections = Vector{Int}[]
    for (cable_index, (design, declaration)) in
        enumerate(zip(declared_designs, declarations))
        names = design.terminal_order
        values = if declaration === nothing
            [index == 1 ? cable_index : 0 for index in eachindex(names)]
        elseif declaration isa AbstractVector{<:Integer}
            length(declaration) == length(names) || throw(DimensionMismatch(
                "connection vector for cable $cable_index must match its terminal count"
            ))
            Int[declaration...]
        elseif declaration isa NamedTuple
            unknown = setdiff(collect(keys(declaration)), names)
            isempty(unknown) || throw(KeyError(first(unknown)))
            Int[Int(get(declaration, name, 0)) for name in names]
        elseif declaration isa AbstractDict
            normalized = Dict{Symbol, Int}()
            for (name, value) in declaration
                key = Symbol(name)
                key in names || throw(KeyError(key))
                value isa Integer || throw(ArgumentError(
                    "connection assignments must be integers"
                ))
                normalized[key] = Int(value)
            end
            Int[get(normalized, name, 0) for name in names]
        else
            throw(ArgumentError(
                "each connection declaration must be a mapping or integer vector"
            ))
        end
        all(>=(0), values) || throw(DomainError(
            values,
            "phase assignments must be nonnegative"
        ))
        push!(normalized_connections, values)
    end
    connection_order = collect(Iterators.flatten(normalized_connections))

    # 4. Resolve global geometry and environment relations. The air/earth
    # interface remains `y == 0`; a richer environment may add formulation
    # checks later without changing the physical system declaration.
    global_geometry = ResolvedRegion[]
    terminal_map = Int[]
    for (cable_index, (design, pose)) in enumerate(zip(declared_designs, poses))
        iszero(pose.y) && throw(DomainError(
            pose.y,
            "a cable centre cannot lie on the air/earth interface"
        ))
        radius = outer_radius(design)
        abs(pose.y) >= radius || throw(DomainError(
            pose.y,
            "the cable cross-section crosses the air/earth interface"
        ))
        for (source, local_terminal) in zip(
                design.geometry.regions,
                design.terminal_map
        )
            push!(global_geometry, ResolvedRegion(
                source.region,
                PlacedShape(source.shape, pose),
                source.terminal,
                source.overlength,
                source.turns,
                source.pattern_depth,
                source.path_depth
            ))
            push!(
                terminal_map,
                iszero(local_terminal) ? 0 :
                terminal_offsets[cable_index] + local_terminal
            )
        end
    end
    for left in eachindex(declared_designs)
        for right in (left + 1):length(declared_designs)
            distance = hypot(
                poses[left].x - poses[right].x,
                poses[left].y - poses[right].y
            )
            limit = outer_radius(declared_designs[left]) +
                    outer_radius(declared_designs[right])
            tolerance = oftype(limit, 1e-8) * max(limit, one(limit))
            distance + tolerance < limit && throw(DomainError(
                (left, right),
                "cable cross-sections overlap"
            ))
        end
    end

    # 5–6. Freeze declarations and all current downstream ordering together.
    return LineCableSystem{
        T,
        typeof(declared_designs),
        typeof(poses),
        typeof(normalized_connections),
        typeof(environment),
        typeof(global_geometry)
    }(
        identifier,
        convert(T, float(line_length)),
        declared_designs,
        poses,
        normalized_connections,
        environment,
        global_geometry,
        terminal_order,
        terminal_map,
        connection_order
    )
end

function LineCableSystem(
        design::CableDesign;
        position = Pose2(0, -1, 0),
        kwargs...
)
    return LineCableSystem(CableDesign[design]; positions = Pose2[position], kwargs...)
end

function Base.show(io::IO, ::MIME"text/plain", system::LineCableSystem)
    println(
        io,
        "LineCableSystem \"$(system.system_id)\": [line_length=$(system.line_length), " *
        "cables=$(ncables(system)), terminals=$(length(system.terminal_order))]"
    )
    for (index, (design, position, connections)) in enumerate(zip(
            system.designs,
            system.positions,
            system.connections
    ))
        prefix = index == ncables(system) ? "└─" : "├─"
        mapping = join(
            ["$terminal→$phase" for (terminal, phase) in
             zip(design.terminal_order, connections)],
            ", "
        )
        println(
            io,
            "$prefix CableDesign \"$(design.cable_id)\": " *
            "[x=$(round(position.x, sigdigits=4)), " *
            "y=$(round(position.y, sigdigits=4)), φ=$(round(position.φ, sigdigits=4)), " *
            "connections=($mapping)]"
        )
    end
end

include("dataframe.jl")
