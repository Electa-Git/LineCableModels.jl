"""
$(TYPEDEF)

Arrange physical members while retaining their independent terminal identities.

$(TYPEDFIELDS)
"""
struct Assembly{A, E <: AbstractCablePart, P, H, C, N} <: AbstractCablePart
    "Pose relative to the containing frame."
    at::A
    "Repeated physical member."
    item::E
    "Cross-sectional placement definition."
    pattern::P
    "Longitudinal path definition or `nothing`."
    path::H
    "Compaction definition or `nothing`."
    compact::C
    "Member-name prefix or exact member names."
    names::N

    function Assembly(
            at::A,
            item::E,
            pattern::P,
            path::H,
            compact::C,
            names::N
    ) where {A, E <: AbstractCablePart, P, H, C, N}
        at isa Pose2 || throw(ArgumentError("assembly pose must resolve to Pose2"))
        names isa Union{Symbol, AbstractVector{Symbol}, Tuple{Vararg{Symbol}}} ||
            throw(ArgumentError("assembly names must be a Symbol or symbol collection"))
        names isa Symbol && isempty(String(names)) &&
            throw(ArgumentError("assembly name prefix cannot be empty"))
        names isa Symbol || all(name -> !isempty(String(name)), names) ||
            throw(ArgumentError("assembly names cannot be empty"))
        return new{A, E, P, H, C, N}(at, item, pattern, path, compact, names)
    end
end

function resolve(context::EmptyBoundary, assembly::Assembly)
    child = resolve(EmptyBoundary(), assembly.item)
    poses = placements(assembly.pattern, child, assembly.compact)
    count = length(poses)
    count > 0 || throw(ArgumentError("assembly placement cannot be empty"))
    member_names = assembly.names isa Symbol ?
                   [Symbol(assembly.names, :_, index) for index in 1:count] :
                   collect(assembly.names)
    length(member_names) == count || throw(DimensionMismatch(
        "assembly requires $count names; got $(length(member_names))"
    ))
    allunique(member_names) || throw(ArgumentError("assembly names must be unique"))

    regions = ResolvedRegion[]
    terminals = Symbol[]
    for (member_name, pose) in zip(member_names, poses)
        placed_at = assembly.at * pose
        radius = hypot(pose.x, pose.y)
        path_factor = assembly.path === nothing ? one(eltype(placed_at)) :
                      overlength(assembly.path, radius)
        local_turns = assembly.path === nothing ? zero(eltype(placed_at)) :
                      inv(pitch(assembly.path, radius))
        pattern_depth = assembly.pattern === nothing ? 0 : 1
        path_depth = assembly.path === nothing ? 0 : 1
        terminal_map = Dict{Symbol, Symbol}()
        if length(child.terminals) == 1
            terminal_map[only(child.terminals)] = member_name
            push!(terminals, member_name)
        else
            for terminal in child.terminals
                qualified = Symbol(member_name, :__, terminal)
                terminal_map[terminal] = qualified
                push!(terminals, qualified)
            end
        end
        for source in child.regions
            terminal = source.terminal
            if terminal !== nothing
                terminal = terminal_map[terminal]
            elseif source.region.material.kind === :conductor
                terminal = member_name
                member_name in terminals || push!(terminals, member_name)
            end
            push!(regions, ResolvedRegion(
                source.region,
                PlacedShape(source.shape, placed_at),
                terminal,
                source.overlength * path_factor,
                source.turns + local_turns,
                source.pattern_depth + pattern_depth,
                source.path_depth + path_depth
            ))
        end
    end
    allunique(terminals) || throw(ArgumentError("assembly terminal names must be unique"))

    child_extent = support(boundary(child))
    outer_radius = maximum(hypot(pose.x, pose.y) + child_extent for pose in poses)
    local_boundary = DiskShape(outer_radius)
    return ResolvedPart(regions, PlacedShape(local_boundary, assembly.at), terminals)
end

function resolve(context::AbstractShape, assembly::Assembly)
    throw(ArgumentError(
        "an Assembly requires explicit placement inside a Stack or Enclosure"
    ))
end
