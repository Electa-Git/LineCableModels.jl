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
    ) where {A, S <: AbstractPrimitive, E <: AbstractCablePart, F, W}
        isempty(String(tag)) && throw(ArgumentError("enclosure tag cannot be empty"))
        at isa Pose2 || throw(ArgumentError("enclosure pose must resolve to Pose2"))
        fill isa Union{Material, Region} ||
            throw(ArgumentError("enclosure fill must be Material or Region"))
        wall isa Union{Nothing, AbstractCablePart} ||
            throw(ArgumentError("enclosure wall must be a cable part or nothing"))
        return new{A, S, E, F, W}(tag, at, primitive, item, fill, wall)
    end
end

Base.:(==)(left::Enclosure, right::Enclosure) =
    left.tag == right.tag && left.at == right.at &&
    left.primitive == right.primitive && left.item == right.item &&
    left.fill == right.fill && left.wall == right.wall

function resolve(
        container::Disk,
        holes::Tuple,
        material::Material,
        tag::Symbol
)
    length(holes) == 1 || return _difference_fill(container, holes, material, tag)
    return _disk_fill(container, only(holes), material, tag)
end

function _disk_fill(
        container::Disk,
        contents::Union{Disk, Annulus},
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


_disk_fill(container::Disk, contents, material::Material, tag::Symbol) =
    _difference_fill(container, (contents,), material, tag)

function _difference_fill(container, holes, material, tag)
    source = Region(Symbol(tag, :_fill), _definition(container), material)
    primitive = _DifferencePrimitive(container, holes)
    area(primitive) > zero(eltype(primitive)) || throw(DomainError(
        area(primitive), "enclosure fill must have positive area"
    ))
    return CableGeometry(PlacedRegion[PlacedRegion(source, primitive)], boundary(container))
end

resolve(container::AbstractShape, holes::Tuple, material::Material, tag::Symbol) =
    _difference_fill(container, holes, material, tag)

_definition(primitive::Disk) = Disk(primitive.r)
_definition(primitive::Rectangle) = Rectangle(primitive.w, primitive.h)
_definition(primitive::Ellipse) = Ellipse(primitive.a, primitive.b)
_definition(primitive::Sector) = Sector(
    primitive.ri, primitive.ro, primitive.φ0, primitive.span
)
_definition(primitive::Annulus) = Annulus(primitive.ri, primitive.ro)
_definition(primitive::Polygon) = Polygon(primitive.points)

resolve(
    ::AbstractPrimitive,
    holes::Tuple,
    fill::Region,
    ::Symbol
) = length(holes) == 1 ? resolve(only(holes), fill) : throw(ArgumentError(
    "an explicit fill Region requires one occupied boundary"
))

function _occupied_boundaries(item::AbstractCablePart)
    return (boundary(resolve(EmptyBoundary(), item)),)
end

function _occupied_boundaries(
        assembly::Assembly{<:Any, <:AbstractCablePart}
)
    child = resolve(EmptyBoundary(), assembly.item)
    return Tuple(
        resolve(assembly.at * _placement_pose(pose), boundary(child))
        for pose in placements(assembly.pattern, child, assembly.compact)
    )
end

function _occupied_boundaries(assembly::Assembly{<:Any, <:Tuple})
    return Tuple(
        resolve(
            assembly.at * member.at,
            boundary(resolve(EmptyBoundary(), member.item))
        )
        for member in assembly.item
    )
end

function _contained(container::Disk, child::AbstractShape)
    return support(child) <= container.r
end

function _contained(container::Rectangle, child::AbstractShape)
    return support(child, 0) <= container.w / 2 &&
           support(child, pi) <= container.w / 2 &&
           support(child, pi / 2) <= container.h / 2 &&
           support(child, -pi / 2) <= container.h / 2
end

function _contained(::AbstractShape, ::AbstractShape)
    throw(ArgumentError(
        "enclosure containment is not implemented for this resolved shape pair"
    ))
end

function resolve(context::EmptyBoundary, enclosure::Enclosure)
    container = resolve(EmptyBoundary(), enclosure.primitive)
    contents = resolve(EmptyBoundary(), enclosure.item)
    holes = _occupied_boundaries(enclosure.item)
    all(hole -> _contained(container, hole), holes) || throw(DomainError(
        holes,
        "enclosed geometry must fit inside the enclosure boundary"
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
        holes,
        enclosure.fill,
        enclosure.tag
    )
    outer_extent = support(container)
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

function resolve(context::AbstractShape, enclosure::Enclosure)
    throw(ArgumentError(
        "an Enclosure requires explicit placement inside a Stack or Enclosure"
    ))
end
