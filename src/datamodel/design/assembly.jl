"Store one explicitly placed member of an `Assembly`."
struct AssemblyMember{E <: AbstractCablePart, P <: Pose2}
    item::E
    at::P
end

"""
Store the resolved, potentially disconnected boundary of an `Assembly`.

The member boundaries remain exact. In particular, an assembly of sector
cores is not replaced by a circular envelope. Consumers that require one
containing domain must request an explicit `Enclosure`.
"""
struct AssemblyShape{T <: Real, S <: Tuple} <: AbstractShape{T}
    members::S

    function AssemblyShape(members::S) where {S <: Tuple}
        isempty(members) && throw(ArgumentError(
            "an assembly boundary requires at least one member"
        ))
        all(member -> member isa AbstractShape, members) || throw(ArgumentError(
            "assembly boundary members must be resolved shapes"
        ))
        T = promote_type(map(eltype, members)...)
        return new{T, S}(members)
    end
end

boundary(shape::AssemblyShape) = shape
area(shape::AssemblyShape) = sum(area, shape.members)
perimeter(shape::AssemblyShape) = sum(perimeter, shape.members)
function support(shape::AssemblyShape, angle::Real)
    maximum(member -> support(member, angle), shape.members)
end
support(shape::AssemblyShape) = maximum(support, shape.members)

function centroid(shape::AssemblyShape)
    areas = map(area, shape.members)
    total = sum(areas)
    centres = map(centroid, shape.members)
    return (
        sum(index -> areas[index] * centres[index][1], eachindex(areas)) / total,
        sum(index -> areas[index] * centres[index][2], eachindex(areas)) / total
    )
end

function resolve(at::Pose2, shape::AssemblyShape)
    AssemblyShape(map(member -> resolve(at, member), shape.members))
end

AssemblyMember(item::AbstractCablePart) = AssemblyMember(item, Pose2(0, 0, 0))

function Base.:(==)(left::AssemblyMember, right::AssemblyMember)
    left.item == right.item && left.at == right.at
end

"""
$(TYPEDEF)

Arrange physical members while retaining their independent terminal
identities.

An assembly stores either one prototype with a repetition pattern or an
explicit tuple of independently placed members. Both forms resolve through
the same terminal-preserving operation.

$(TYPEDFIELDS)
"""
struct Assembly{A, E, P, H, C, N} <: AbstractCablePart
    "Pose relative to the containing frame."
    at::A
    "Repeated prototype or tuple of explicit members."
    item::E
    "Cross-sectional placement definition for a repeated prototype."
    pattern::P
    "Longitudinal path definition or `nothing`."
    path::H
    "Compaction definition or `nothing`."
    compact::C
    "Exact repeated-member terminal names, or `nothing` for explicit members."
    names::N

    function Assembly(
            at::A,
            item::E,
            pattern::P,
            path::H,
            compact::C,
            names::N
    ) where {A, E, P, H, C, N}
        at isa Pose2 || throw(ArgumentError("assembly pose must resolve to Pose2"))
        repeated = item isa AbstractCablePart
        explicit = item isa Tuple && !isempty(item) &&
                   all(member -> member isa AssemblyMember, item)
        repeated || explicit ||
            throw(ArgumentError(
                "assembly item must be one prototype or explicit placed members"
            ))
        if repeated
            pattern === nothing && throw(ArgumentError(
                "a repeated assembly requires a placement pattern"
            ))
            names isa Union{Nothing, AbstractVector{Symbol}, Tuple{Vararg{Symbol}}} ||
                throw(ArgumentError(
                    "assembly names must be exact symbols or nothing"
                ))
            names === nothing ||
                all(name -> !isempty(String(name)), names) ||
                throw(ArgumentError(
                    "assembly names cannot be empty"
                ))
        else
            pattern === nothing && path === nothing && compact === nothing &&
            names === nothing || throw(ArgumentError(
                "explicit assembly members own their poses and terminal identities"
            ))
        end
        return new{A, E, P, H, C, N}(at, item, pattern, path, compact, names)
    end
end

function Base.:(==)(left::Assembly, right::Assembly)
    left.at == right.at && left.item == right.item &&
        left.pattern == right.pattern && left.path == right.path &&
        left.compact == right.compact && left.names == right.names
end

function _append_assembly_member!(
        regions::Vector{PlacedRegion},
        child::CableGeometry,
        assembly_at::Pose2,
        member_at::Pose2;
        terminal_map = identity,
        pattern = nothing,
        member::Int = 1,
        path = nothing
)
    radius = hypot(member_at.x, member_at.y)
    extent = zero(support(boundary(child)))
    for source in child.regions
        placed = resolve(assembly_at * member_at, source)
        terminal = source.terminal === nothing ? nothing : terminal_map(source.terminal)
        local_primitive = resolve(member_at, source.primitive)
        extent = max(extent, support(local_primitive))
        patterns = pattern === nothing ? placed.placement.patterns :
                   (placed.placement.patterns...,
            (pattern = pattern, member = member, pose = member_at))
        paths = path === nothing ? source.paths :
                (source.paths..., (path = path, radius = radius))
        push!(regions,
            PlacedRegion(
                source.source,
                placed.primitive,
                terminal,
                (patterns = patterns,),
                paths
            ))
    end
    return extent
end

function _resolve_repeated(assembly::Assembly)
    child = resolve(EmptyBoundary(), assembly.item)
    child_terminals = unique(Symbol[source.terminal
                                    for source in child.regions
                                    if source.terminal !== nothing])
    length(child_terminals) <= 1 || throw(ArgumentError(
        "a repeated assembly prototype may resolve at most one terminal"
    ))
    poses = placements(assembly.pattern, child, assembly.compact)
    count = length(poses)
    count > 0 || throw(ArgumentError("assembly placement cannot be empty"))
    member_names = if isempty(child_terminals)
        assembly.names === nothing || throw(ArgumentError(
            "a zero-terminal repeated assembly does not accept terminal names"
        ))
        fill(nothing, count)
    else
        assembly.names === nothing && throw(ArgumentError(
            "a terminal-bearing repeated assembly requires exact names"
        ))
        values = collect(assembly.names)
        length(values) == count || throw(DimensionMismatch(
            "assembly requires $count names; got $(length(values))"
        ))
        allunique(values) || throw(ArgumentError("assembly names must be unique"))
        values
    end

    regions = PlacedRegion[]
    boundaries = AbstractShape[]
    for (member, (name, placement)) in enumerate(zip(member_names, poses))
        pose = _placement_pose(placement)
        mapping = isempty(child_terminals) ? identity :
                  terminal -> terminal === only(child_terminals) ? name : terminal
        _append_assembly_member!(
            regions,
            child,
            assembly.at,
            pose;
            terminal_map = mapping,
            pattern = assembly.pattern,
            member,
            path = assembly.path
        )
        push!(boundaries, resolve(assembly.at * pose, boundary(child)))
    end
    return CableGeometry(regions, AssemblyShape(Tuple(boundaries)))
end

function _resolve_explicit(assembly::Assembly)
    regions = PlacedRegion[]
    terminals = Symbol[]
    boundaries = AbstractShape[]
    for member in assembly.item
        child = resolve(EmptyBoundary(), member.item)
        child_terminals = unique(Symbol[source.terminal
                                        for source in child.regions
                                        if source.terminal !== nothing])
        for terminal in child_terminals
            terminal in terminals && throw(ArgumentError(
                "explicit assembly terminal :$terminal is not unique"
            ))
            push!(terminals, terminal)
        end
        _append_assembly_member!(
            regions, child, assembly.at, member.at
        )
        push!(boundaries, resolve(assembly.at * member.at, boundary(child)))
    end
    return CableGeometry(regions, AssemblyShape(Tuple(boundaries)))
end

function resolve(
        ::EmptyBoundary,
        assembly::Assembly{<:Any, <:AbstractCablePart}
)
    _resolve_repeated(assembly)
end

resolve(
    ::EmptyBoundary,
    assembly::Assembly{<:Any, <:Tuple}
) = _resolve_explicit(assembly)

function resolve(::AbstractShape, ::Assembly)
    throw(ArgumentError(
        "an Assembly requires explicit placement inside a Stack or Enclosure"
    ))
end
