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

function Base.:(==)(left::Enclosure, right::Enclosure)
    left.tag == right.tag && left.at == right.at &&
        left.primitive == right.primitive && left.item == right.item &&
        left.fill == right.fill && left.wall == right.wall
end

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
        contents::Disk,
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

function _disk_fill(container::Disk, contents, material::Material, tag::Symbol)
    _difference_fill(container, (contents,), material, tag)
end

function _difference_fill(container, holes, material, tag)
    source = Region(Symbol(tag, :_fill), _definition(container), material)
    primitive = DifferenceShape(container, holes)
    area(primitive) > zero(eltype(primitive)) || throw(DomainError(
        area(primitive), "enclosure fill must have positive area"
    ))
    return CableGeometry(PlacedRegion[PlacedRegion(source, primitive)], boundary(container))
end

function resolve(container::AbstractShape, holes::Tuple, material::Material, tag::Symbol)
    _difference_fill(container, holes, material, tag)
end

_definition(primitive::Disk) = Disk(primitive.r)
_definition(primitive::Rectangle) = Rectangle(primitive.w, primitive.h)
_definition(primitive::Ellipse) = Ellipse(primitive.a, primitive.b)
_definition(shape::SectorShape) = shape.primitive
_definition(primitive::Annulus) = Annulus(primitive.ri, primitive.ro)
_definition(primitive::Polygon) = Polygon(primitive.points)

function resolve(
        ::AbstractPrimitive,
        holes::Tuple,
        fill::Region,
        ::Symbol
)
    length(holes) == 1 ? resolve(only(holes), fill) :
    throw(ArgumentError(
        "an explicit fill Region requires one occupied boundary"
    ))
end

function _contained(container::Disk, child::AbstractShape)
    extent = support(child)
    return extent < container.r || isapprox(extent, container.r)
end

function _contained(container::Rectangle, child::AbstractShape)
    return support(child, 0) <= container.w / 2 &&
           support(child, pi) <= container.w / 2 &&
           support(child, pi / 2) <= container.h / 2 &&
           support(child, -pi / 2) <= container.h / 2
end

function _contained(container::Ellipse, child::AbstractShape)
    tolerance = sqrt(eps(float(max(container.a, container.b))))
    angles = range(zero(container.a), oftype(container.a, 2pi); length = 4097)
    return all(Iterators.take(angles, 4096)) do angle
        support(child, angle) <= support(container, angle) + tolerance
    end
end

function _contained(container::Annulus, child::Disk)
    offset = hypot(
        child.at.x - container.at.x,
        child.at.y - container.at.y
    )
    return (offset + child.r < container.ro ||
            isapprox(offset + child.r, container.ro)) &&
           (offset - child.r > container.ri ||
            isapprox(offset - child.r, container.ri))
end

function _contained(container::Annulus, child::Annulus)
    concentric = isapprox(child.at.x, container.at.x) &&
                 isapprox(child.at.y, container.at.y)
    return concentric &&
           (child.ro < container.ro || isapprox(child.ro, container.ro)) &&
           (child.ri > container.ri || isapprox(child.ri, container.ri))
end

function _contained(container::Annulus, child::BentStrip)
    concentric = isapprox(child.at.x, container.at.x) &&
                 isapprox(child.at.y, container.at.y)
    return concentric &&
           (child.ro < container.ro || isapprox(child.ro, container.ro)) &&
           (child.ri > container.ri || isapprox(child.ri, container.ri))
end

function _contained(::AbstractShape, ::AbstractShape)
    throw(ArgumentError(
        "enclosure containment is not implemented for this resolved shape pair"
    ))
end

function fill_holes(::AbstractCablePart, contents::CableGeometry)
    return Tuple(source.primitive for source in contents.regions)
end

fill_holes(::Enclosure, contents::CableGeometry) = (boundary(contents),)

function fill_holes(assembly::Assembly, contents::CableGeometry)
    outer = boundary(contents)
    outer isa AssemblyShape || throw(ArgumentError(
        "an enclosed assembly must retain its member boundaries"
    ))
    return outer.members
end

function resolve(context::EmptyBoundary, enclosure::Enclosure)
    container = resolve(EmptyBoundary(), enclosure.primitive)
    contents = resolve(EmptyBoundary(), enclosure.item)
    holes = enclosure.fill isa Material ?
            fill_holes(enclosure.item, contents) : (boundary(contents),)
    all(hole -> _contained(container, hole), holes) || throw(DomainError(
        holes,
        "enclosed geometry must fit inside the enclosure boundary"
    ))

    regions = PlacedRegion[]
    for source in contents.regions
        placed = resolve(enclosure.at, source)
        push!(regions,
            PlacedRegion(
                source.source,
                placed.primitive,
                source.terminal,
                placed.placement,
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
        placed = resolve(enclosure.at, source)
        push!(regions,
            PlacedRegion(
                source.source,
                placed.primitive,
                source.terminal,
                placed.placement,
                source.paths
            ))
    end

    outer = container
    if enclosure.wall !== nothing
        wall_result = resolve(container, enclosure.wall)
        for source in wall_result.regions
            placed = resolve(enclosure.at, source)
            push!(regions,
                PlacedRegion(
                    source.source,
                    placed.primitive,
                    source.terminal,
                    placed.placement,
                    source.paths
                ))
        end
        outer = boundary(wall_result)
    end
    return CableGeometry(regions, resolve(enclosure.at, boundary(outer)))
end

function resolve(context::AbstractShape, enclosure::Enclosure)
    container = resolve(EmptyBoundary(), enclosure.primitive)
    container isa Annulus || throw(ArgumentError(
        "a contextual Enclosure requires an annular containing primitive"
    ))
    placed = resolve(enclosure.at, container)
    context isa Disk &&
    isapprox(context.r, placed.ri) &&
    isapprox(context.at.x, placed.at.x) &&
    isapprox(context.at.y, placed.at.y) || throw(DomainError(
        context,
        "an annular Enclosure must continue the preceding circular boundary"
    ))
    return resolve(EmptyBoundary(), enclosure)
end
