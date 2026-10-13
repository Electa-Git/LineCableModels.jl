"""
$(TYPEDEF)

Mark a resolved fill or wall region that establishes an enclosure boundary.

Geometry methods use this marker to identify an enclosure boundary.
"""
struct EnclosureBoundary end

# The marker has no geometry, so it is unchanged when its region moves.
resolve(::Pose2, marker::EnclosureBoundary) = marker

"""
$(TYPEDEF)

Contain a physical cable object in an explicit cross-section, fill, and wall.

$(TYPEDFIELDS)
"""
struct Enclosure{
    A <: Pose2,
    S <: AbstractPrimitive,
    E <: AbstractCablePart,
    F <: Union{Material, Region},
    W <: Union{Nothing, AbstractCablePart}
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
    ) where {
            A <: Pose2,
            S <: AbstractPrimitive,
            E <: AbstractCablePart,
            F <: Union{Material, Region},
            W <: Union{Nothing, AbstractCablePart}
    }
        tag === Symbol("") && throw(ArgumentError("enclosure tag cannot be empty"))
        return new{A, S, E, F, W}(tag, at, primitive, item, fill, wall)
    end
end

"""
Contain a cable part in `primitive` with `fill` and an optional `wall`, or declare a
finite space of such enclosures.
"""
function Enclosure(
        tag::Union{Symbol, AbstractGrid, Gridspace},
        item::Union{AbstractCablePart, AbstractGrid, Gridspace};
        at = Pose2(0, 0, 0),
        primitive,
        fill,
        wall = nothing,
        combine::Symbol = :product
)
    values = (tag, at, primitive, item, fill, wall)
    return parameterize(Enclosure, Enclosure, values; combine)
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
    length(holes) == 1 && return _disk_fill(container, only(holes), material, tag)
    # A coaxial terminal may supply a disk plus several annuli as its occupied
    # shapes. Their union is one disk only when the radial intervals are complete.
    # Preserve that exact annular fill instead of manufacturing a Boolean region
    # that computational flattening would classify as non-radial interstitial fill.
    if all(holes) do shape
            shape isa Union{Disk, Annulus} &&
                iszero(shape.at.x) && iszero(shape.at.y) &&
                iszero(container.at.x) && iszero(container.at.y)
        end
        intervals = sort!([(r_in(shape), r_ex(shape)) for shape in holes]; by = first)
        outer = zero(container.r)
        tolerance = 100 * eps(typeof(float(_geometry_scalar(container.r)))) * container.r
        for (inner, radius) in intervals
            inner <= outer + tolerance || return _difference_fill(container, holes, material, tag)
            outer = max(outer, radius)
        end
        return _disk_fill(container, Disk(outer), material, tag)
    end
    return _difference_fill(container, holes, material, tag)
end

function _disk_fill(
        container::Disk,
        contents::Disk,
        material::Material,
        tag::Symbol
)
    # Only contents at the center leave an annular gap. Off-center contents leave a
    # crescent.
    contents.at.x == container.at.x && contents.at.y == container.at.y ||
        return _difference_fill(container, (contents,), material, tag)
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
    tolerance = 5.0e-6 * container.r
    return extent <= container.r + tolerance
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

function _contained(container::SectorShape, child::Disk)
    accommodates(container, centroid(child), child.r)
end

function _contained(container::SectorShape, child::Polygon)
    cosine = cos(child.at.φ)
    sine = sin(child.at.φ)
    tolerance = 5.0e-6 * container.primitive.r_back
    return all(child.points) do point
        resolved = (
            child.at.x + cosine * point[1] - sine * point[2],
            child.at.y + sine * point[1] + cosine * point[2]
        )
        _geometry_scalar(clearance(container, resolved) + tolerance) >= 0
    end
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

# The holes that a material fill leaves are those of the enclosed item. A fill Region has
# a declared geometry, and it surrounds the outline of the contents.
function fill_holes(::Material, item::AbstractCablePart, contents::CableGeometry)
    return fill_holes(item, contents)
end
fill_holes(::Region, ::AbstractCablePart, contents::CableGeometry) = (boundary(contents),)

function fill_holes(group::Group, contents::CableGeometry)
    if group.boundary isa Disk && last(contents.regions).source.primitive isa Rectangle
        outer = boundary(contents)
        center = first(contents.regions).primitive
        if center isa Disk && center.at.x == outer.at.x && center.at.y == outer.at.y
            # The formation already resolved its exact occupied disk.
            return (outer,)
        end
    end
    return Tuple(source.primitive for source in contents.regions)
end

# An assembly, grouped or enclosed, keeps the outline of each member.
function fill_holes(
        ::Group{<:Pose2, <:Assembly},
        contents::CableGeometry{<:AbstractVector{<:PlacedRegion}, <:AssemblyShape}
)
    return boundary(contents).members
end
function fill_holes(::Group{<:Pose2, <:Assembly}, ::CableGeometry)
    throw(ArgumentError("a grouped assembly must retain its member boundaries"))
end
function fill_holes(
        ::Assembly,
        contents::CableGeometry{<:AbstractVector{<:PlacedRegion}, <:AssemblyShape}
)
    return boundary(contents).members
end
function fill_holes(::Assembly, ::CableGeometry)
    throw(ArgumentError("an enclosed assembly must retain its member boundaries"))
end

# An explicit fill Region always places its geometry, and a material fills only the area
# that the contents leave.
_has_fill(::Region, remaining_area, tolerance) = true
_has_fill(::Material, remaining_area, tolerance) = remaining_area > tolerance

function resolve(context::EmptyBoundary, enclosure::Enclosure)
    container = resolve(EmptyBoundary(), enclosure.primitive)
    contents = resolve(EmptyBoundary(), enclosure.item)
    return _resolve_enclosure(enclosure, container, contents)
end

function _resolve_enclosure(enclosure::Enclosure, container, contents::CableGeometry)
    holes = fill_holes(enclosure.fill, enclosure.item, contents)
    all(hole -> _contained(container, hole), holes) || throw(DomainError(
        holes,
        "enclosed geometry must fit inside the enclosure boundary"
    ))
    # A complete circular partition is one occupied disk even through terminal
    # and stack wrappers. Incomplete formations retain their individual holes.
    occupied = boundary(contents)
    if container isa Disk && occupied isa Disk &&
       occupied.at.x == container.at.x && occupied.at.y == container.at.y &&
       isapprox(sum(area, holes), area(occupied);
            rtol=0, atol=geometry_tolerance(area(occupied)))
        holes = (occupied,)
    end

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

    remaining_area = area(container) - sum(area, holes; init = zero(area(container)))
    tolerance = geometry_tolerance(area(container))
    remaining_area >= -tolerance || throw(DomainError(
        remaining_area, "enclosure contents exceed the containing boundary area"
    ))
    if _has_fill(enclosure.fill, remaining_area, tolerance)
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
                    (
                        patterns = (
                            placed.placement.patterns...,
                            (
                                owner = Enclosure,
                                pattern = EnclosureBoundary(),
                                member = 1,
                                pose = enclosure.at
                            )
                        ),
                    ),
                    source.paths
                ))
        end
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
                    (
                        patterns = (
                            placed.placement.patterns...,
                            (
                                owner = Enclosure,
                                pattern = EnclosureBoundary(),
                                member = 1,
                                pose = enclosure.at
                            )
                        ),
                    ),
                    source.paths
                ))
        end
        outer = boundary(wall_result)
    end
    return CableGeometry(regions, resolve(enclosure.at, boundary(outer)))
end

function resolve(context::Disk, enclosure::Enclosure{<:Pose2, <:Annulus})
    container = resolve(EmptyBoundary(), enclosure.primitive)
    placed = resolve(enclosure.at, container)
    validate(context, enclosure, placed)
    inner = Disk(context.r, container.at)
    contents = resolve(inner, enclosure.item)
    effective = Annulus(context.r, container.ro, container.at)
    return _resolve_enclosure(enclosure, effective, contents)
end

function resolve(context::AbstractShape, ::Enclosure{<:Pose2, <:Annulus})
    throw(DomainError(
        context,
        "an annular Enclosure must be concentric with and outside the preceding circular boundary"
    ))
end

function resolve(::AbstractShape, ::Enclosure)
    throw(ArgumentError(
        "a contextual Enclosure requires an annular containing primitive"
    ))
end

# The preceding circular boundary of an annular Enclosure lies inside its container and
# shares its center.
function validate(context::Disk, enclosure::Enclosure, container::Annulus)
    isapprox(context.at.x, container.at.x) &&
    isapprox(context.at.y, container.at.y) &&
    context.r <= container.ri + geometry_tolerance(container.ri) || throw(DomainError(
        context,
        "an annular Enclosure must be concentric with and outside the preceding circular boundary"
    ))
    validate(context, enclosure.fill, container)
    return context
end

# A filled course defines its complement down to the preceding physical boundary.
# Explicit wire radii remain fixed. Contextual rings use this disk. An explicit fill
# Region has its own geometry and is not extended, so it must continue that disk.
validate(context::Disk, ::Material, ::Annulus) = context
function validate(context::Disk, ::Region, container::Annulus)
    isapprox(context.r, container.ri;
        rtol=0, atol=geometry_tolerance(container.ri)) || throw(DomainError(
        context, "an explicit fill Region must continue the preceding boundary"))
    return context
end
