"""
$(TYPEDEF)

Contain a physical cable object in an explicit cross-section, fill, and wall.

$(TYPEDFIELDS)
"""
struct Enclosure{
        A,
        S <: AbstractPrimitiveDefinition,
        E <: AbstractCablePart,
        F,
        W
} <: AbstractCablePart
    "Physical enclosure identity."
    tag::Symbol
    "Pose relative to the containing frame."
    at::A
    "Containing cross-sectional primitive."
    primitive::S
    "Enclosed physical object."
    item::E
    "Filling material or filling Region."
    fill::F
    "Optional physical wall."
    wall::W

    function Enclosure(
            tag::Symbol,
            at::A,
            primitive::S,
            item::E,
            fill::F,
            wall::W
    ) where {A, S <: AbstractPrimitiveDefinition, E <: AbstractCablePart, F, W}
        isempty(String(tag)) && throw(ArgumentError("enclosure tag cannot be empty"))
        at isa Pose2 || throw(ArgumentError("enclosure pose must resolve to Pose2"))
        fill isa Union{Material, Region} ||
            throw(ArgumentError("enclosure fill must be Material or Region"))
        wall isa Union{Nothing, AbstractCablePart} ||
            throw(ArgumentError("enclosure wall must be a cable part or nothing"))
        return new{A, S, E, F, W}(tag, at, primitive, item, fill, wall)
    end
end

function resolve(
        container::Disk,
        contents::AbstractPrimitive,
        material::Material,
        tag::Symbol
)
    fill_region = Region(
        Symbol(tag, :_fill),
        AnnulusDefinition(support(contents), container.r),
        material
    )
    return resolve(EmptyBoundary(), fill_region)
end

function resolve(
        ::AbstractPrimitive,
        ::AbstractPrimitive,
        ::Material,
        ::Symbol
)
    throw(ArgumentError("material fill currently requires a circular enclosure"))
end

resolve(
    ::AbstractPrimitive,
    contents::AbstractPrimitive,
    fill::Region,
    ::Symbol
) = resolve(contents, fill)

function resolve(context::EmptyBoundary, enclosure::Enclosure)
    container = resolve(EmptyBoundary(), enclosure.primitive)
    contents = resolve(EmptyBoundary(), enclosure.item)
    inner_extent = support(boundary(contents))
    outer_extent = support(container)
    inner_extent < outer_extent || throw(DomainError(
        inner_extent,
        "enclosed geometry must fit strictly inside the enclosure boundary"
    ))

    regions = PlacedRegion[]
    for source in contents.regions
        primitive = resolve(enclosure.at, source.primitive)
        push!(regions, PlacedRegion(
            source.source,
            primitive,
            source.terminal,
            (patterns = source.placement.patterns,),
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
        primitive = resolve(enclosure.at, source.primitive)
        push!(regions, PlacedRegion(
            source.source,
            primitive,
            source.terminal,
            (patterns = source.placement.patterns,),
            source.paths
        ))
    end

    outer = container
    if enclosure.wall !== nothing
        wall_result = resolve(container, enclosure.wall)
        for source in wall_result.regions
            primitive = resolve(enclosure.at, source.primitive)
            push!(regions, PlacedRegion(
                source.source,
                primitive,
                source.terminal,
                (patterns = source.placement.patterns,),
                source.paths
            ))
        end
        outer = boundary(wall_result)
    end
    return CableGeometry(regions, resolve(enclosure.at, outer))
end

function resolve(context::AbstractPrimitive, enclosure::Enclosure)
    throw(ArgumentError(
        "an Enclosure requires explicit placement inside a Stack or Enclosure"
    ))
end
