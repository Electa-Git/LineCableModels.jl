"""
$(TYPEDEF)

Store completed cable placements and their global terminal state.

`designs`, `declared_positions`, `connections`, and `environment` are declarations.
`positions` contains the resolved poses after automatic exterior-clearance
adjustment. Touching cables are separated by at least 1 μm plus the propagated
uncertainty reserve; genuinely overlapping nominal designs are rejected.
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
    "Cable poses before automatic clearance adjustment \\[m, m, rad\\]."
    declared_positions::P
    "Retained pairwise clearance requirements; diagonal entries concern the interface \\[m\\]."
    clearances::Matrix{T}
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
            connection_order::Vector{Int},
            declared_positions::P,
            clearances::Matrix{T}
    ) where {
            T <: Real,
            D <: AbstractVector,
            P <: AbstractVector{<:Pose2{T}},
            C <: AbstractVector,
            E,
            G <: AbstractVector
    }
        return validate(new{T, D, P, C, E, G}(
            system_id,
            line_length,
            designs,
            positions,
            declared_positions,
            clearances,
            connections,
            environment,
            geometry,
            terminal_order,
            terminal_map,
            connection_order
        ))
    end
end

"""
$(TYPEDSIGNATURES)

Return a system whose exterior circles satisfy their retained clearance from
the air-earth interface. Reuse the input when no adjustment is necessary.
This construction step is called when a line problem supplies the earth model.
"""
function interface_clearance(system::LineCableSystem)
    reference = system.declared_positions
    context = _CLEARANCE_CONTEXT[]
    sampling = context !== nothing && context.sampling[]
    reference_centres = sampling ? get(context.references, system.clearances, nothing) : nothing
    positions, _, displacement = clearance_geometry(system.designs, system.positions;
        required = system.clearances, reference, reference_centres, interface = true, sampling)
    iszero(displacement) && return system
    result = build(LineCableSystem, system.designs, positions, system.connections,
        system.environment, system.system_id, system.line_length;
        _clearances = system.clearances, _declared_positions = system.declared_positions,
        _interface = true, _rebuild = true)
    sampling && reference_centres !== nothing &&
        (context.references[result.clearances] = reference_centres)
    _record_clearance_adjustment(system.system_id, displacement)
    return result
end

Base.eltype(::LineCableSystem{T}) where {T} = T
Base.eltype(::Type{<:LineCableSystem{T}}) where {T} = T

function realize(rng::Random.AbstractRNG, point::Gridpoint{LineCableSystem}, distribution)
    return realize_clearance(rng, point, distribution)
end

ncables(system::LineCableSystem) = length(system.designs)
nphases(system::LineCableSystem) = length(unique(filter(>(0), system.connection_order)))

function validate(system::LineCableSystem)
    isempty(system.system_id) && throw(ArgumentError(
        "LineCableSystem.system_id cannot be empty"
    ))
    isfinite(system.line_length) && system.line_length > zero(system.line_length) ||
        throw(DomainError(
            system.line_length,
            "LineCableSystem.line_length must be positive and finite"
        ))
    isempty(system.designs) && throw(ArgumentError(
        "LineCableSystem.designs must contain at least one CableDesign"
    ))
    length(system.positions) == length(system.designs) || throw(DimensionMismatch(
        "LineCableSystem.positions must contain one Pose2 per design; received " *
        "$(length(system.positions)) positions for $(length(system.designs)) designs"
    ))
    length(system.declared_positions) == length(system.designs) || throw(DimensionMismatch(
        "LineCableSystem.declared_positions must contain one pose per cable"))
    length(system.connections) == length(system.designs) || throw(DimensionMismatch(
        "LineCableSystem.connections must contain one declaration per design; " *
        "received $(length(system.connections)) for $(length(system.designs)) designs"
    ))
    expected_terminals = sum(length(design.terminal_order) for design in system.designs)
    length(system.terminal_order) == expected_terminals || throw(DimensionMismatch(
        "LineCableSystem.terminal_order must contain $expected_terminals entries; " *
        "received $(length(system.terminal_order))"
    ))
    length(system.connection_order) == expected_terminals || throw(DimensionMismatch(
        "LineCableSystem.connection_order must contain $expected_terminals entries; " *
        "received $(length(system.connection_order))"
    ))
    expected_regions = sum(length(design.geometry.regions) for design in system.designs)
    length(system.geometry) == expected_regions || throw(DimensionMismatch(
        "LineCableSystem.geometry must contain $expected_regions resolved regions; " *
        "received $(length(system.geometry))"
    ))
    length(system.terminal_map) == expected_regions || throw(DimensionMismatch(
        "LineCableSystem.terminal_map must contain $expected_regions entries; " *
        "received $(length(system.terminal_map))"
    ))

    terminal_index = 0
    region_index = 0
    terminal_offset = 0
    for (cable_index, design) in pairs(system.designs)
        validate(design)
        pose = system.positions[cable_index]
        all(isfinite, (pose.x, pose.y, pose.φ)) || throw(DomainError(
            (pose.x, pose.y, pose.φ),
            "LineCableSystem.positions[$cable_index] must contain finite coordinates"
        ))
        connection = system.connections[cable_index]
        connection isa AbstractVector{<:Integer} || throw(ArgumentError(
            "LineCableSystem.connections[$cable_index] must be an integer vector; " *
            "received $(typeof(connection))"
        ))
        length(connection) == length(design.terminal_order) ||
            throw(DimensionMismatch(
                "LineCableSystem.connections[$cable_index] must contain " *
                "$(length(design.terminal_order)) terminal assignments; received " *
                "$(length(connection))"
            ))
        all(>=(0), connection) || throw(DomainError(
            connection,
            "LineCableSystem.connections[$cable_index] must be nonnegative"
        ))
        for (local_terminal, name) in pairs(design.terminal_order)
            terminal_index += 1
            expected = (cable = cable_index, terminal = name)
            system.terminal_order[terminal_index] == expected ||
                throw(DimensionMismatch(
                    "LineCableSystem.terminal_order[$terminal_index] must be " *
                    "$(repr(expected)); received " *
                    "$(repr(system.terminal_order[terminal_index]))"
                ))
            system.connection_order[terminal_index] == connection[local_terminal] ||
                throw(DimensionMismatch(
                    "LineCableSystem.connection_order[$terminal_index] must match " *
                    "connections[$cable_index][$local_terminal]"
                ))
        end
        for (local_region, placed) in pairs(design.geometry.regions)
            region_index += 1
            resolved = system.geometry[region_index]
            resolved.source == placed.source || throw(DimensionMismatch(
                "LineCableSystem.geometry[$region_index].source does not match " *
                "designs[$cable_index].geometry.regions[$local_region].source"
            ))
            resolved.terminal === placed.terminal || throw(DimensionMismatch(
                "LineCableSystem.geometry[$region_index].terminal does not match " *
                "designs[$cable_index].geometry.regions[$local_region].terminal"
            ))
            local_index = design.terminal_map[local_region]
            expected = iszero(local_index) ? 0 : terminal_offset + local_index
            system.terminal_map[region_index] == expected || throw(DimensionMismatch(
                "LineCableSystem.terminal_map[$region_index] must be $expected; " *
                "received $(system.terminal_map[region_index])"
            ))
        end
        terminal_offset += length(design.terminal_order)
    end

    clearance_geometry(system.designs, system.positions;
        required = system.clearances, interface = system.environment !== nothing,
        adjust = false)
    return system
end

function build(
        ::Type{LineCableSystem},
        designs,
        placements,
        connections,
        environment,
        system_id::AbstractString,
        line_length::Real;
        combine::Symbol = :product,
        _clearances = nothing,
        _declared_positions = nothing,
        _interface::Bool = environment !== nothing,
        _rebuild::Bool = false
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
        Pose2[position isa Pose2 ? position :
              position isa Tuple && length(position) in (2, 3) &&
              all(value -> value isa Real, position) ?
              Pose2(
                  position[1],
                  position[2],
                  length(position) == 3 ? position[3] : 0
              ) :
              throw(ArgumentError(
                  "placements must contain Pose2 values or `(x, y[, φ])` tuples"
              ))
              for position in placements]
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
    original_poses = _declared_positions === nothing ? copy(poses) :
                     Pose2{T}[convert(Pose2{T}, position) for position in _declared_positions]
    context = _CLEARANCE_CONTEXT[]
    sampling = context !== nothing && context.sampling[]
    reference = _rebuild ? poses : original_poses
    reference_centres = nothing
    if context !== nothing && !_rebuild && sampling
        context.cursor[] += 1
        context.cursor[] <= length(context.records) || throw(ArgumentError(
            "sampled construction produced more cable systems than its declaration"))
        record = context.records[context.cursor[]]
        record.system_id == identifier && record.cables == getproperty.(declared_designs, :cable_id) ||
            throw(ArgumentError("sampled cable-system layout differs from its declaration"))
        _clearances = record.clearances
        reference = record.positions
        reference_centres = record.centres
    end
    poses, clearances, displacement = clearance_geometry(declared_designs, poses;
        required = _clearances, reference, reference_centres, interface = _interface, sampling)
    if context !== nothing && !_rebuild && !sampling
        push!(context.records, (system_id = identifier,
            cables = getproperty.(declared_designs, :cable_id),
            clearances = nominal.(clearances), positions = original_poses,
            centres = getproperty.(_clearance_exterior.(declared_designs, original_poses), :centre)))
    end

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
    for (cable_index, (design, declaration)) in enumerate(zip(declared_designs, declarations))
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

    # 4. Place the resolved cable geometry. An explicit environment owns any
    # interface constraint; formulation-specific media checks occur in the
    # problem that consumes the system.
    global_geometry = PlacedRegion[]
    terminal_map = Int[]
    for (cable_index, (design, pose)) in enumerate(zip(declared_designs, poses))
        for (source, local_terminal) in zip(
            design.geometry.regions,
            design.terminal_map
        )
            push!(global_geometry, resolve(pose, source))
            push!(
                terminal_map,
                iszero(local_terminal) ? 0 :
                terminal_offsets[cable_index] + local_terminal
            )
        end
    end
    # 5. Freeze declarations and their completed ordering together.
    system = LineCableSystem{
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
        connection_order,
        original_poses,
        clearances
    )
    if sampling && !_rebuild
        context.references[system.clearances] = reference_centres
    end
    _record_clearance_adjustment(identifier, displacement)
    return system
end
