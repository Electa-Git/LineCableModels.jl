"""
$(TYPEDEF)

Supertype of the declared parts of a cable: a region, or a node that layers, names, repeats,
lays or encloses parts.
"""
abstract type AbstractCablePart end

# A part in its own frame resolves from the center, where nothing lies inside it yet.
function resolve(::EmptyBoundary, part::AbstractCablePart)
    return resolve(
        (regions = PlacedRegion[], outer = EmptyBoundary(), terminals = NamedTuple[],
            repeats = NamedTuple[], lays = NamedTuple[], fills = NamedTuple[]),
        part)
end

"""
$(TYPEDEF)

Declare a part of one material in one cross-section.

$(TYPEDFIELDS)
"""
struct Region{P <: Union{AbstractPrimitive, Offset}, M <: Material} <: AbstractCablePart
    "Name of the region within the cable."
    tag::Symbol
    "Declared cross-section: a primitive, or an offset from the outline of the part inside."
    primitive::P
    "Material of the region. Its class tells whether the region conducts."
    material::M
    function Region(tag::Symbol, primitive::P,
            material::M) where {P <: Union{AbstractPrimitive, Offset}, M <: Material}
        return new{P, M}(tag, primitive, material)
    end
end

# A region resolves its primitive against the outline of the parts inside it.
function resolve(described::NamedTuple, region::Region)
    shape = resolve(described.outer, region.primitive)
    return merge(described,
        (regions = [described.regions; PlacedRegion(region, shape)],
            outer = boundary(shape)))
end

"""
$(TYPEDEF)

Layer parts outward. Each part resolves against the outline of the parts inside it.

$(TYPEDFIELDS)
"""
struct Stack{V <: AbstractVector{<:AbstractCablePart}} <: AbstractCablePart
    "Parts from the inside out."
    items::V
    function Stack(items::V) where {V <: AbstractVector{<:AbstractCablePart}}
        return validate(new{V}(items))
    end
end

"""
$(TYPEDSIGNATURES)

Layer `items` outward, from the inside out.

# Errors

- Throws `ArgumentError` when no part is given.
"""
Stack(items::AbstractCablePart...) = Stack(AbstractCablePart[items...])

function validate(stack::Stack)
    isempty(stack.items) && throw(ArgumentError("Stack.items must hold at least one part"))
    return stack
end

resolve(described::NamedTuple, stack::Stack) = foldl(resolve, stack.items; init = described)

"""
$(TYPEDEF)

Assign the conductive regions of a part to one terminal, which its name identifies within
its scope. The regions of one terminal form one node of the circuit. A terminal around a part
takes the conductive regions of every terminal inside it, the copies of a repeat included,
and the parts of one name in one scope share one terminal.

$(TYPEDFIELDS)
"""
struct Terminal{P <: AbstractCablePart} <: AbstractCablePart
    "Name that identifies the terminal within its scope."
    name::Symbol
    "Part whose conductive regions belong to the terminal."
    part::P
    function Terminal(name::Symbol, part::P) where {P <: AbstractCablePart}
        return new{P}(name, part)
    end
end

# The conductive regions of the part leave every terminal inside it and join the terminal of
# its name in the scope of the part, the one with an empty scope. A terminal left without
# regions has no row.
function resolve(described::NamedTuple, terminal::Terminal)
    inner = resolve(described, terminal.part)
    conductive = [index
                  for index in (length(described.regions) + 1):length(inner.regions)
                  if inner.regions[index].source.material isa Material{<:Real, :conductor}]
    named = [index
             for row in inner.terminals if row.name == terminal.name && isempty(row.scope)
             for index in row.regions]
    others = [(name = row.name, scope = row.scope,
                  regions = setdiff(row.regions, conductive))
              for row in inner.terminals
              if !(row.name == terminal.name && isempty(row.scope))]
    rows = [others;
            (name = terminal.name, scope = NamedTuple[],
                regions = sort!(union(named, conductive)))]
    terminals = sort!([row for row in rows if !isempty(row.regions)];
        by = row -> first(row.regions))
    return merge(inner, (terminals = terminals,))
end

"""
$(TYPEDEF)

Repeat a part at the placements of a pattern. The part is declared in its own frame, and the
pattern places each copy about the center of the enclosing part. A repeated region is a
strand, and the pattern bends a rectangular strand into an annular sector of the same area
on its course. Each copy is a scope: the terminals declared in the part are distinct
terminals in every copy. The outline of the copies is the circle about the center that
contains them.

$(TYPEDFIELDS)
"""
struct Repeat{P <: AbstractLayoutPattern, R <: AbstractCablePart} <: AbstractCablePart
    "Pattern that places the copies."
    pattern::P
    "Part that the pattern copies."
    part::R
    function Repeat(pattern::P, part::R) where {
            P <: AbstractLayoutPattern, R <: AbstractCablePart}
        return new{P, R}(pattern, part)
    end
end

# A copy at the center of the cable has no part inside it.
validate(shape::AbstractShape, ::Repeat, ::EmptyBoundary) = shape

# A copy lies outside the outline of the parts inside the repeat.
function validate(shape::AbstractShape, repeat::Repeat, outline::AbstractShape)
    radius = support(outline)
    tolerance = sqrt(eps(float(radius))) * max(one(radius), radius)
    distance(shape) + tolerance >= radius || throw(DomainError(distance(shape),
        "a copy in a repeat must not enter the part inside it, of outer radius $radius m"))
    return shape
end

# The pattern places a strand, and it bends a rectangular strand into an annular sector on
# its course.
function resolve(described::NamedTuple, repeat::Repeat{<:AbstractLayoutPattern, <:Region})
    course = placements(resolve(described.outer, repeat.pattern, repeat.part.primitive),
        repeat.part.primitive)
    shapes = [validate(resolve(at, course.shape), repeat, described.outer)
              for at in course.placements]
    start = length(described.regions)
    row = (placements = course.placements,
        copies = [index:index for index in (start + 1):(start + length(shapes))])
    return merge(described,
        (
            regions = [described.regions;
                       [PlacedRegion(repeat.part, shape) for shape in shapes]],
            outer = Circle(maximum(support, shapes)), repeats = [described.repeats; row]))
end

# The pattern places a composite part by its outline. Every copy repeats the rows of the part,
# with its regions placed and its indices shifted. The terminals of a copy take the copy as
# their outermost scope, and its lays take the placed axis. The row of the repeat precedes
# those of its copies.
function resolve(described::NamedTuple, repeat::Repeat)
    member = resolve(EmptyBoundary(), repeat.part)
    course = placements(resolve(described.outer, repeat.pattern, member.outer), member.outer)
    index = length(described.repeats) + 1
    offsets = [(at = at,
                   region = length(described.regions) + (k - 1) * length(member.regions),
                   repeat = index + (k - 1) * length(member.repeats))
               for (k, at) in enumerate(course.placements)]
    regions = [PlacedRegion(placed.source,
                   validate(resolve(offset.at, placed.shape), repeat, described.outer))
               for offset in offsets for placed in member.regions]
    terminals = [(name = row.name,
                     scope = [(repeat = index, copy = k);
                              [(repeat = entry.repeat + offset.repeat, copy = entry.copy)
                               for entry in row.scope]],
                     regions = row.regions .+ offset.region)
                 for (k, offset) in enumerate(offsets) for row in member.terminals]
    repeats = [(placements = course.placements,
                   copies = [(offset.region + 1):(offset.region + length(member.regions))
                             for offset in offsets]);
               [(placements = [offset.at * at for at in row.placements],
                    copies = [indices .+ offset.region for indices in row.copies])
                for offset in offsets for row in member.repeats]]
    lays = [(repeat = row.repeat + offset.repeat,
                radii = row.radii, pitches = row.pitches,
                overlengths = row.overlengths, axis = offset.at(row.axis))
            for offset in offsets for row in member.lays]
    fills = [(region = row.region + offset.region,
                 contents = row.contents .+ offset.region)
             for offset in offsets for row in member.fills]
    return merge(described,
        (regions = [described.regions; regions],
            outer = Circle(maximum(placed -> support(placed.shape), regions)),
            terminals = [described.terminals; terminals],
            repeats = [described.repeats; repeats], lays = [described.lays; lays],
            fills = [described.fills; fills]))
end

"""
$(TYPEDEF)

Lay the copies of a repeat along the cable. Each copy follows the helix about the center of
the repeat, at the distance of its centroid.

$(TYPEDFIELDS)
"""
struct Lay{H <: Helix, R <: Repeat} <: AbstractCablePart
    "Helix that the copies follow."
    helix::H
    "Repeat whose copies the helix lays."
    part::R
    function Lay(helix::H, part::R) where {H <: Helix, R <: Repeat}
        return new{H, R}(helix, part)
    end
end

# The row of the laid repeat is the first that the repeat adds.
function resolve(described::NamedTuple, lay::Lay)
    inner = resolve(described, lay.part)
    index = length(described.repeats) + 1
    radii = map(inner.repeats[index].copies) do indices
        shapes = [inner.regions[region].shape for region in indices]
        x, y = sum(shape -> area(shape) .* centroid(shape), shapes) ./ sum(area, shapes)
        return hypot(x, y)
    end
    row = (repeat = index, radii = radii,
        pitches = [lay.helix.dir * pitch(lay.helix, radius) for radius in radii],
        overlengths = [overlength(lay.helix, radius) for radius in radii],
        axis = (zero(eltype(radii)), zero(eltype(radii))))
    return merge(inner, (lays = [inner.lays; row],))
end

"""
$(TYPEDEF)

Fill a container around parts with a material. The fill region declares the container as its
primitive, either explicitly or as an offset from the part inside the enclosure, and takes
it without the regions of the contents. The outline of the enclosure is the outline of its
container. The container is a circle or an annulus, and the regions inside it lie between its
inner and outer radii. `resolve` throws `ArgumentError` for any other container.

$(TYPEDFIELDS)
"""
struct Enclosure{F <: Region, C <: AbstractCablePart} <: AbstractCablePart
    "Region of the fill material, whose primitive declares the container."
    fill::F
    "Parts inside the container."
    contents::C
    function Enclosure(fill::F, contents::C) where {F <: Region, C <: AbstractCablePart}
        return new{F, C}(fill, contents)
    end
end

# A circle or an annulus contains the regions between its inner and outer radii.
validate(container::Union{Circle, Annulus}, ::Enclosure) = container

function validate(container::AbstractShape, enclosure::Enclosure)
    throw(ArgumentError(
        "the container of :$(enclosure.fill.tag) is a $(nameof(typeof(container))); " *
        "an enclosure takes a circle or an annulus as its container"))
end

# A region inside an enclosure lies inside its container.
function validate(shape::AbstractShape, enclosure::Enclosure, container::AbstractShape)
    radius = support(container)
    tolerance = sqrt(eps(float(radius))) * max(one(radius), radius)
    inside = support(shape) <= radius + tolerance &&
             distance(shape) + tolerance >= distance(container)
    inside || throw(DomainError(support(shape),
        "a region inside :$(enclosure.fill.tag) must lie between " *
        "$(distance(container)) m and $radius m from its center"))
    return shape
end

function resolve(described::NamedTuple, enclosure::Enclosure)
    container = validate(resolve(described.outer, enclosure.fill.primitive), enclosure)
    inner = resolve(described, enclosure.contents)
    contents = (length(described.regions) + 1):length(inner.regions)
    holes = [validate(inner.regions[index].shape, enclosure, container)
             for index in contents]
    placed = PlacedRegion(enclosure.fill, DifferenceShape(container, holes))
    row = (region = length(inner.regions) + 1, contents = collect(contents))
    return merge(inner,
        (regions = [inner.regions; placed], outer = boundary(container),
            fills = [inner.fills; row]))
end
