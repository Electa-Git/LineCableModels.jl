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

function resolve(
        container::DiskShape,
        contents::AbstractShape,
        material::Material,
        tag::Symbol
)
    fill_region = Region(
        Symbol(tag, :_fill),
        Annulus(support(contents), container.r),
        material
    )
    return resolve(EmptyBoundary(), fill_region)
end

function resolve(
        ::AbstractShape,
        ::AbstractShape,
        ::Material,
        ::Symbol
)
    throw(ArgumentError("material fill currently requires a circular enclosure"))
end

resolve(
    ::AbstractShape,
    contents::AbstractShape,
    fill::Region,
    ::Symbol
) = resolve(contents, fill)

function resolve(context::EmptyBoundary, enclosure::Enclosure)
    container = resolve(EmptyBoundary(), enclosure.shape)
    contents = resolve(EmptyBoundary(), enclosure.item)
    inner_extent = support(boundary(contents))
    outer_extent = support(container)
    inner_extent < outer_extent || throw(DomainError(
        inner_extent,
        "enclosed geometry must fit strictly inside the enclosure boundary"
    ))

    regions = PlacedRegion[]
    for source in contents.regions
        shape = PlacedShape(source.shape, enclosure.at)
        push!(regions, PlacedRegion(
            source.source,
            shape,
            source.terminal,
            (pose = shape.at, patterns = source.placement.patterns),
            source.paths
        ))
    end

    fill_result = resolve(
        container,
        boundary(contents),
        enclosure.fill,
        enclosure.tag
    )
    isapprox(support(boundary(fill_result)), outer_extent) || throw(DomainError(
        support(boundary(fill_result)),
        "enclosure fill must reach the containing boundary"
    ))
    for source in fill_result.regions
        shape = PlacedShape(source.shape, enclosure.at)
        push!(regions, PlacedRegion(
            source.source,
            shape,
            source.terminal,
            (pose = shape.at, patterns = source.placement.patterns),
            source.paths
        ))
    end

    outer = container
    if enclosure.wall !== nothing
        wall_result = resolve(container, enclosure.wall)
        for source in wall_result.regions
            shape = PlacedShape(source.shape, enclosure.at)
            push!(regions, PlacedRegion(
                source.source,
                shape,
                source.terminal,
                (pose = shape.at, patterns = source.placement.patterns),
                source.paths
            ))
        end
        outer = boundary(wall_result)
    end
    return CableGeometry(regions, PlacedShape(outer, enclosure.at))
end

function resolve(context::AbstractShape, enclosure::Enclosure)
    throw(ArgumentError(
        "an Enclosure requires explicit placement inside a Stack or Enclosure"
    ))
end
