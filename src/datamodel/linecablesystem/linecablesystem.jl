"""
$(TYPEDEF)

Store completed cable placements and their global terminal state.

`designs`, `connections`, and `environment` are declarations. `positions` holds the
declared poses, with every exterior gap below 1 μm raised to 1 μm: touching cables are
1 μm apart. Overlapping exteriors throw [`OverlapError`](@ref). Global geometry, terminal
order, terminal indices, and the flattened connection order are derived by the
constructor.

$(TYPEDFIELDS)
"""
struct LineCableSystem{
    T <: Real,
    D <: AbstractVector,
    P <: AbstractVector{<:Pose2{T}},
    C <: AbstractVector{<:AbstractVector{<:Integer}},
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
    "Active phase IDs aligned with `terminal_order`. Zero selects a conductor for elimination."
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
            C <: AbstractVector{<:AbstractVector{<:Integer}},
            E,
            G <: AbstractVector
    }
        return validate(new{T, D, P, C, E, G}(
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
        ))
    end
end

line_length(system::LineCableSystem) = system.line_length

"""
$(TYPEDSIGNATURES)

Return a system whose cable exteriors lie at least 1 μm from the air-earth interface. An
interface gap from zero up to 1 μm is raised to 1 μm, and a cable that crosses the
interface throws [`OverlapError`](@ref). Reuse the input when nothing moves. A line
problem calls it when it supplies the earth model.
"""
function interface_clearance(system::LineCableSystem)
    positions, displacement = clearance_geometry(system.designs, system.positions;
        pairs = false, interface = true)
    iszero(displacement) && return system
    return build(LineCableSystem, system.designs, positions, system.connections,
        system.environment, system.system_id, system.line_length)
end

Base.eltype(::LineCableSystem{T}) where {T} = T
Base.eltype(::Type{<:LineCableSystem{T}}) where {T} = T

Base.convert(::Type{LineCableSystem{T}}, system::LineCableSystem{T}) where {T <: Real} = system

# A wider scalar type rebuilds the system through its own `build`, from its built
# positions.
function Base.convert(::Type{LineCableSystem{T}}, system::LineCableSystem) where {T <: Real}
    return build(LineCableSystem, system.designs, system.positions, system.connections,
        system.environment, system.system_id, convert(T, system.line_length))::LineCableSystem{T}
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
        pose = system.positions[cable_index]
        all(isfinite, (pose.x, pose.y, pose.φ)) || throw(DomainError(
            (pose.x, pose.y, pose.φ),
            "LineCableSystem.positions[$cable_index] must contain finite coordinates"
        ))
        connection = system.connections[cable_index]
        length(connection) == length(design.terminal_order) ||
            throw(DimensionMismatch(
                "LineCableSystem.connections[$cable_index] must contain " *
                "$(length(design.terminal_order)) terminal assignments; received " *
                "$(length(connection))"
            ))
        all(>=(0), connection) || throw(DomainError(
            connection,
            "LineCableSystem.connections[$cable_index] must contain nonnegative " *
            "active phase IDs or zero"
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
        interface = system.environment !== nothing, adjust = false)
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
        combine::Symbol = :product
)
    combine in (:product, :zip) || throw(ArgumentError(
        "combine must be :product or :zip"
    ))
    declared_input = _declared_designs(designs)

    # 1. Place every completed design in the system frame. Coordinate tuples
    # are construction shorthand. The stored declaration is always Pose2.
    isempty(declared_input) && throw(ArgumentError("a line system requires one cable"))
    identifier = String(system_id)
    isempty(identifier) && throw(ArgumentError("system_id cannot be empty"))
    isfinite(line_length) && line_length > zero(line_length) || throw(DomainError(
        line_length,
        "line length must be positive"
    ))
    position_values = _declared_poses(placements)
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
    positions, _ = clearance_geometry(declared_designs, poses;
        interface = environment !== nothing)

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
    declarations = _declared_connections(connections, length(declared_designs))
    normalized_connections = Vector{Int}[]
    for (cable_index, (design, declaration)) in enumerate(zip(declared_designs, declarations))
        values = _connection_values(declaration, design.terminal_order, cable_index)
        all(>=(0), values) || throw(DomainError(
            values,
            "connection assignments must contain nonnegative active phase IDs or zero"
        ))
        push!(normalized_connections, values)
    end
    connection_order = collect(Iterators.flatten(normalized_connections))

    # 4. Place the resolved cable geometry. An explicit environment owns any
    # interface constraint. Formulation-specific media checks occur in the
    # problem that consumes the system.
    global_geometry = PlacedRegion[]
    terminal_map = Int[]
    for (cable_index, (design, pose)) in enumerate(zip(declared_designs, positions))
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
    return LineCableSystem{
        T,
        typeof(declared_designs),
        typeof(positions),
        typeof(normalized_connections),
        typeof(environment),
        typeof(global_geometry)
    }(
        identifier,
        convert(T, float(line_length)),
        declared_designs,
        positions,
        normalized_connections,
        environment,
        global_geometry,
        terminal_order,
        terminal_map,
        connection_order
    )
end

# The designs of a system in placement order: one completed design, or a collection of them.
_declared_designs(design::CableDesign) = typeof(design)[design]
function _declared_designs(designs::Union{AbstractVector, Tuple})
    return [_declared_design(design) for design in designs]
end
function _declared_designs(designs)
    throw(ArgumentError(
        "designs must be a CableDesign or a design collection; " *
        "got a $(nameof(typeof(designs)))"
    ))
end
_declared_design(design::CableDesign) = design
function _declared_design(design)
    throw(ArgumentError(
        "designs must contain completed CableDesign objects; " *
        "got a $(nameof(typeof(design)))"
    ))
end

# The placements of a system are one pose or coordinate tuple `(x, y[, φ])`, or a
# collection of them. Coordinate tuples are construction shorthand, and the stored
# declaration is always Pose2.
function _declared_poses(
        placement::Union{Pose2, Tuple{Real, Real}, Tuple{Real, Real, Real}}
)
    return Pose2[_declared_pose(placement)]
end
function _declared_poses(placements::Union{AbstractVector, Tuple})
    return Pose2[_declared_pose(position) for position in placements]
end
function _declared_poses(placements)
    throw(ArgumentError(
        "placements must be a Pose2, coordinate tuple, or placement collection; " *
        "got a $(nameof(typeof(placements)))"
    ))
end
_declared_pose(position::Pose2) = position
_declared_pose(position::Tuple{Real, Real}) = Pose2(position[1], position[2], 0)
_declared_pose(position::Tuple{Real, Real, Real}) = Pose2(position[1], position[2], position[3])
function _declared_pose(position)
    throw(ArgumentError(
        "placements must contain Pose2 values or `(x, y[, φ])` tuples; " *
        "got a $(nameof(typeof(position)))"
    ))
end

# Each design has one connection declaration. With one design, the declaration may come
# alone.
_declared_connections(::Nothing, count::Int) = fill(nothing, count)
function _declared_connections(connections::Union{AbstractVector, Tuple}, count::Int)
    length(connections) == count || throw(DimensionMismatch(
        "connection declaration count must match the cable-design count"
    ))
    return collect(connections)
end
function _declared_connections(declaration::AbstractVector{<:Integer}, count::Int)
    count == 1 && return Any[declaration]
    length(declaration) == count || throw(DimensionMismatch(
        "connection declaration count must match the cable-design count"
    ))
    return collect(declaration)
end
function _declared_connections(declaration::Union{AbstractDict, NamedTuple}, count::Int)
    count == 1 || throw(ArgumentError(
        "connections must be one declaration per cable design"
    ))
    return Any[declaration]
end
function _declared_connections(connections, count::Int)
    throw(ArgumentError(
        "connections must be one declaration per cable design"
    ))
end

# The active phase ID of each terminal of one design, zero for elimination. Without a
# declaration, the design's first terminal takes the cable index.
function _connection_values(::Nothing, names, cable_index::Int)
    return [index == 1 ? cable_index : 0 for index in eachindex(names)]
end
function _connection_values(declaration::AbstractVector, names, cable_index::Int)
    length(declaration) == length(names) || throw(DimensionMismatch(
        "connection vector for cable $cable_index must match its terminal count"
    ))
    return Int[_connection_assignment(value) for value in declaration]
end
function _connection_values(declaration::NamedTuple, names, cable_index::Int)
    unknown = setdiff(collect(keys(declaration)), names)
    isempty(unknown) || throw(KeyError(first(unknown)))
    return Int[Int(get(declaration, name, 0)) for name in names]
end
function _connection_values(declaration::AbstractDict, names, cable_index::Int)
    normalized = Dict{Symbol, Int}()
    for (name, value) in declaration
        key = Symbol(name)
        key in names || throw(KeyError(key))
        normalized[key] = _connection_assignment(value)
    end
    return Int[get(normalized, name, 0) for name in names]
end
function _connection_values(declaration, names, cable_index::Int)
    throw(ArgumentError(
        "each connection declaration must be a mapping or integer vector; " *
        "got a $(nameof(typeof(declaration)))"
    ))
end
_connection_assignment(value::Integer) = Int(value)
function _connection_assignment(value)
    throw(ArgumentError(
        "connection assignments must be integers; got a $(nameof(typeof(value)))"
    ))
end

"""
$(TYPEDSIGNATURES)

Place completed cable designs in one system. Without keywords, the system takes the
identifier `"line-cable-system"` and a length of 1 m, and it has no connections or
environment.
Returns a `LineCableSystem`, or a `Gridspace{LineCableSystem}` when an input is a `Grid` or
`Gridspace`.
"""
function build(
        ::Type{LineCableSystem},
        designs,
        placements;
        connections = nothing,
        environment = nothing,
        system_id = "line-cable-system",
        line_length = 1,
        combine::Symbol = :product
)
    values = (designs, placements, connections, environment, system_id, line_length)
    return parameterize(
        LineCableSystem, (values...) -> build(LineCableSystem, values...), values; combine)
end

"""
$(TYPEDSIGNATURES)

Build one system from placed cables. Each placed cable is a `(design, pose, connections)`
record, and a formation is a collection of them that contributes its records in order.
Returns a `LineCableSystem`, or a `Gridspace{LineCableSystem}` when an input is a `Grid` or
`Gridspace`.
"""
function build(
        ::Type{LineCableSystem},
        placements::Union{Tuple, AbstractVector, AbstractGrid, Gridspace};
        environment = nothing,
        system_id = "line-cable-system",
        line_length = 1,
        combine::Symbol = :product
)
    entries(placement) = (placement,)
    entries(formation::Union{Tuple, AbstractVector}) = formation
    caller = function (placements, environment, system_id, line_length)
        records = [validate(record, LineCableSystem)
                   for placement in placements for record in entries(placement)]
        isempty(records) && throw(ArgumentError(
            "placed cable declarations must be a nonempty collection"))
        return build(LineCableSystem, getproperty.(records, :design),
            getproperty.(records, :pose), getproperty.(records, :connections),
            environment, system_id, line_length)
    end
    values = (placements, environment, system_id, line_length)
    return parameterize(LineCableSystem, caller, values; combine)
end

# A placed cable checked for the system build.
function validate(
        record::NamedTuple{(:design, :pose, :connections), <:Tuple{CableDesign, Pose2, Any}},
        ::Type{LineCableSystem}
)
    return record
end

function validate(record, ::Type{LineCableSystem})
    throw(ArgumentError("each placed cable requires design, pose, and connections; " *
                        "got a $(nameof(typeof(record)))"))
end

Commons.input_fields(::Type{<:LineCableSystem}) = (line_length=(name="line length",unit="m"),
    positions=(name="position",unit="m"))
