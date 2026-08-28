"""
$(TYPEDEF)

Contain a physical cable object in an explicit cross-section, fill, and wall.

$(TYPEDFIELDS)
"""
struct Enclosure{
        A,
        S <: AbstractPrimitive,
        E <: AbstractCablePart,
        F,
        W
} <: AbstractCablePart
    "Physical enclosure identity."
    tag::Symbol
    "Pose relative to the containing frame."
    at::A
    "Containing cross-sectional primitive."
    shape::S
    "Enclosed physical object."
    item::E
    "Filling material or filling Region."
    fill::F
    "Optional physical wall."
    wall::W

    function Enclosure(
            tag::Symbol,
            at::A,
            shape::S,
            item::E,
            fill::F,
            wall::W
    ) where {A, S <: AbstractPrimitive, E <: AbstractCablePart, F, W}
        isempty(String(tag)) && throw(ArgumentError("enclosure tag cannot be empty"))
        at isa Pose2 || throw(ArgumentError("enclosure pose must resolve to Pose2"))
        fill isa Union{Material, Region} ||
            throw(ArgumentError("enclosure fill must be Material or Region"))
        wall isa Union{Nothing, AbstractCablePart} ||
            throw(ArgumentError("enclosure wall must be a cable part or nothing"))
        return new{A, S, E, F, W}(tag, at, shape, item, fill, wall)
    end
end

function resolve(context::EmptyBoundary, enclosure::Enclosure)
    container = resolve(EmptyBoundary(), enclosure.shape)
    contents = resolve(EmptyBoundary(), enclosure.item)
    inner_extent = support(boundary(contents))
    outer_extent = support(container)
    inner_extent < outer_extent || throw(DomainError(
        inner_extent,
        "enclosed geometry must fit strictly inside the enclosure boundary"
    ))

    regions = ResolvedRegion[]
    for source in contents.regions
        push!(regions, ResolvedRegion(
            source.region,
            PlacedShape(source.shape, enclosure.at),
            source.terminal,
            source.overlength,
            source.turns,
            source.pattern_depth,
            source.path_depth
        ))
    end

    fill_result = if enclosure.fill isa Material
        container isa DiskShape || throw(ArgumentError(
            "material fill currently requires a circular enclosure"
        ))
        fill_region = Region(
            Symbol(enclosure.tag, :_fill),
            Annulus(inner_extent, container.r),
            enclosure.fill
        )
        resolve(EmptyBoundary(), fill_region)
    else
        resolve(boundary(contents), enclosure.fill)
    end
    isapprox(support(boundary(fill_result)), outer_extent) || throw(DomainError(
        support(boundary(fill_result)),
        "enclosure fill must reach the containing boundary"
    ))
    for source in fill_result.regions
        push!(regions, ResolvedRegion(
            source.region,
            PlacedShape(source.shape, enclosure.at),
            source.terminal,
            source.overlength,
            source.turns,
            source.pattern_depth,
            source.path_depth
        ))
    end

    outer = container
    if enclosure.wall !== nothing
        wall_result = resolve(container, enclosure.wall)
        for source in wall_result.regions
            push!(regions, ResolvedRegion(
                source.region,
                PlacedShape(source.shape, enclosure.at),
                source.terminal,
                source.overlength,
                source.turns,
                source.pattern_depth,
                source.path_depth
            ))
        end
        outer = boundary(wall_result)
    end
    return ResolvedPart(
        regions,
        PlacedShape(outer, enclosure.at),
        copy(contents.terminals)
    )
end

function resolve(context::AbstractShape, enclosure::Enclosure)
    throw(ArgumentError(
        "an Enclosure requires explicit placement inside a Stack or Enclosure"
    ))
end
