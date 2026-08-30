"""
$(TYPEDEF)

Store completed cable placements and their global terminal state.

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
    "Completed cable designs in placement order."
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

    function LineCableSystem{T, D, P, C, E, G}(
            system_id::String,
            line_length::T,
            designs::D,
            positions::P,
            connections::C,
            environment::E,
            geometry::G,
            terminal_order::Vector{
                NamedTuple{(:cable, :terminal), Tuple{Int, Symbol}}
            },
            terminal_map::Vector{Int},
            connection_order::Vector{Int}
    ) where {
            T <: Real,
            D <: AbstractVector,
            P <: AbstractVector{<:Pose2{T}},
            C <: AbstractVector,
            E,
            G <: AbstractVector
    }
        return new{T, D, P, C, E, G}(
            system_id,
            line_length,
            designs,
            positions,
            connections,
            environment,
            geometry,
            terminal_order,
            terminal_map,
            connection_order
        )
    end
end

Base.eltype(::LineCableSystem{T}) where {T} = T
Base.eltype(::Type{<:LineCableSystem{T}}) where {T} = T

ncables(system::LineCableSystem) = length(system.designs)
nphases(system::LineCableSystem) =
    length(unique(filter(>(0), system.connection_order)))

function build(
        ::Type{LineCableSystem},
        designs,
        placements,
        connections,
        environment,
        system_id::AbstractString,
        line_length::Real;
        combine::Symbol = :product
)
    combine in (:product, :zip) || throw(ArgumentError(
        "combine must be :product or :zip"
    ))
    declared_input = if designs isa CableDesign
        typeof(designs)[designs]
    elseif designs isa AbstractVector || designs isa Tuple
        all(design -> design isa CableDesign, designs) || throw(ArgumentError(
            "designs must contain completed CableDesign objects"
        ))
        collect(designs)
    else
        throw(ArgumentError("designs must be a CableDesign or a design collection"))
    end

    # 1. Place every completed design in the system frame. Coordinate tuples
    # are construction shorthand; the stored declaration is always Pose2.
    isempty(declared_input) && throw(ArgumentError("a line system requires one cable"))
    identifier = String(system_id)
    isempty(identifier) && throw(ArgumentError("system_id cannot be empty"))
    isfinite(line_length) && line_length > zero(line_length) || throw(DomainError(
        line_length,
        "line length must be positive"
    ))
    position_values = if placements isa Pose2
        Pose2[placements]
    elseif placements isa Tuple && length(placements) in (2, 3) &&
            all(value -> value isa Real, placements)
        Pose2[Pose2(
            placements[1],
            placements[2],
            length(placements) == 3 ? placements[3] : 0
        )]
    elseif placements isa AbstractVector || placements isa Tuple
        Pose2[
            position isa Pose2 ? position :
            position isa Tuple && length(position) in (2, 3) &&
            all(value -> value isa Real, position) ?
            Pose2(
                position[1],
                position[2],
                length(position) == 3 ? position[3] : 0
            ) : throw(ArgumentError(
                "placements must contain Pose2 values or `(x, y[, φ])` tuples"
            ))
            for position in placements
        ]
    else
        throw(ArgumentError(
            "placements must be a Pose2, coordinate tuple, or placement collection"
        ))
    end
    declared_designs = if length(declared_input) == 1 && length(position_values) > 1
        fill(only(declared_input), length(position_values))
    else
        declared_input
    end
    length(position_values) == length(declared_designs) || throw(DimensionMismatch(
        "placement count must match the cable-design count"
    ))

    T = promote_type(
        typeof(float(line_length)),
        (eltype(design) for design in declared_designs)...,
        (eltype(position) for position in position_values)...
    )
    poses = Pose2{T}[convert(Pose2{T}, position) for position in position_values]

    # 2. Establish global primitive and terminal order while retaining cable
    # and local terminal order verbatim.
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
        fill(nothing, length(declared_designs))
    elseif length(declared_designs) == 1 &&
            (connections isa AbstractDict || connections isa NamedTuple ||
             connections isa AbstractVector{<:Integer})
        Any[connections]
    elseif connections isa AbstractVector || connections isa Tuple
        length(connections) == length(declared_designs) || throw(DimensionMismatch(
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

    # 4. Place the resolved cable geometry and validate its relation to the
    # air/earth interface and the other cable cross-sections.
    global_geometry = PlacedRegion[]
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
            primitive = resolve(pose, source.primitive)
            push!(global_geometry, PlacedRegion(
                source.source,
                primitive,
                source.terminal,
                (patterns = source.placement.patterns,),
                source.paths
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

    # 5. Freeze declarations and their completed ordering together.
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
