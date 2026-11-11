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
            assemblies = NamedTuple[], repeats = NamedTuple[], lays = NamedTuple[], fills = NamedTuple[]),
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
    conductive = Int[index
                     for index in (length(described.regions) + 1):length(inner.regions)
                     if inner.regions[index].source.material isa
                        Material{<:Real, :conductor}]
    named = Int[index
                for row in inner.terminals
                if row.name == terminal.name && isempty(row.scope)
                for index in row.regions]
    others = NamedTuple[(name = row.name, scope = row.scope,
                            regions = setdiff(row.regions, conductive))
                        for row in inner.terminals
                        if !(row.name == terminal.name && isempty(row.scope))]
    rows = NamedTuple[others;
                      (name = terminal.name, scope = Int[],
                          regions = sort!(union(named, conductive)))]
    terminals = sort!(NamedTuple[row for row in rows if !isempty(row.regions)];
        by = row -> first(row.regions))
    return merge(inner, (terminals = terminals,))
end

"""
$(TYPEDEF)

Repeat a part at the placements of a pattern. The part is declared in its own frame, and the
pattern places each copy about the center of the enclosing part. A repeated region is a
strand, and the pattern bends a rectangular strand into an annular sector of the same area
on its course. Each copy is a scope: the terminals declared in the part are distinct
terminals in every copy. A copy of a part other than a region is a sub-cable of the design.
The outline of the copies is the circle about the center that contains them.

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
    distance(shape) + tolerance(radius) >= radius || throw(DomainError(distance(shape),
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

# The pattern places a composite part by its outline. Every copy is a sub-cable and repeats
# the rows of the part, with its regions placed and its indices shifted. The terminals of a
# copy take the row of the copy as their outermost scope, and its lays take the placed axis.
# The row of the repeat precedes those of its copies, and the row of a copy precedes those of
# the sub-cables inside it.
function resolve(described::NamedTuple, repeat::Repeat)
    member = resolve(EmptyBoundary(), repeat.part)
    course = placements(resolve(described.outer, repeat.pattern, member.outer), member.outer)
    index = length(described.repeats) + 1
    m = length(member.assemblies)
    base = length(described.assemblies)
    offsets = [(at = at,
                   region = length(described.regions) +
                            (k - 1) * length(member.regions),
                   repeat = index + (k - 1) * length(member.repeats),
                   assembly = base + (k - 1) * (m + 1) + 1)
               for (k, at) in enumerate(course.placements)]
    regions = [PlacedRegion(placed.source,
                   validate(resolve(offset.at, placed.shape), repeat, described.outer))
               for offset in offsets for placed in member.regions]
    rows = NamedTuple[]
    for offset in offsets
        push!(rows,
            (at = offset.at,
                regions = Int[index + offset.region
                              for index in eachindex(member.regions)
                              if all(row -> index ∉ row.regions, member.assemblies)],
                members = Int[index + offset.assembly
                              for index in eachindex(member.assemblies)
                              if all(row -> index ∉ row.members, member.assemblies)]))
        append!(rows,
            NamedTuple[(at = offset.at * row.at,
                           regions = row.regions .+ offset.region, members = row.members .+
                                                                             offset.assembly)
                       for row in member.assemblies])
    end
    terminals = NamedTuple[(name = row.name,
                               scope = Int[offset.assembly; row.scope .+ offset.assembly],
                               regions = row.regions .+ offset.region)
                           for offset in offsets
                           for row in member.terminals]
    repeats = NamedTuple[(placements = course.placements,
                             copies = UnitRange{Int}[(offset.region + 1):(offset.region + length(member.regions))
                                                     for offset in offsets]);
                         NamedTuple[(
                                        placements = [offset.at * at
                                                      for at in row.placements],
                                        copies = UnitRange{Int}[indices .+ offset.region
                                                                for indices in row.copies])
                                    for offset in offsets for row in member.repeats]]
    lays = NamedTuple[(
                          regions = UnitRange{Int}[indices .+ offset.region
                                                   for indices in row.regions],
                          radii = row.radii, pitches = row.pitches,
                          overlengths = row.overlengths, axis = offset.at(row.axis))
                      for offset in offsets for row in member.lays]
    fills = NamedTuple[(region = row.region + offset.region,
                           contents = row.contents .+ offset.region)
                       for offset in offsets for row in member.fills]
    return merge(described,
        (regions = [described.regions; regions],
            outer = Circle(maximum(placed -> support(placed.shape), regions)),
            terminals = NamedTuple[described.terminals; terminals],
            assemblies = NamedTuple[described.assemblies; rows],
            repeats = NamedTuple[described.repeats; repeats], lays = NamedTuple[described.lays;
                                                                                lays],
            fills = NamedTuple[described.fills; fills]))
end

"""
$(TYPEDEF)

Place sub-cables of different parts in the frame of the part that encloses them. Each member
is a part in its own frame, and its placement puts it in the frame of the assembly. Each
member is a scope: the terminals declared in it are distinct terminals of the design. The
members lie apart from each other and outside the part inside the assembly. The outline of
the assembly is the circle about the center that contains its members.

$(TYPEDFIELDS)
"""
struct Assembly{M <: Tuple{Vararg{Pair{<:Placement, <:AbstractCablePart}}}} <:
       AbstractCablePart
    "Members of the assembly, each a placement and the part placed there."
    members::M
    function Assembly(members::M) where {M <: Tuple{Vararg{Pair{
            <:Placement, <:AbstractCablePart}}}}
        return new{M}(members)
    end
end

# A member lies apart from the other members when there is no part inside the assembly.
function validate(shape::AbstractShape, ::Assembly, circle::Circle,
        others::AbstractVector, ::EmptyBoundary)
    shape ⊆ PerforatedShape(circle, others) || throw(DomainError(shape,
        "a member of an assembly must not overlap another member or the part inside the assembly"))
    return shape
end

# A member lies apart from the other members and the part inside the assembly.
function validate(shape::AbstractShape, ::Assembly, circle::Circle,
        others::AbstractVector, inner::AbstractShape)
    T = eltype(circle)
    shape ⊆ PerforatedShape(circle, AbstractShape{T}[inner; others]) ||
        throw(DomainError(shape,
            "a member of an assembly must not overlap another member or the part inside the assembly"))
    return shape
end

"""
$(TYPEDSIGNATURES)

Resolve each member in its own frame and place its regions and rows in the enclosing frame.
"""
function resolve(described::NamedTuple, assembly::Assembly)
    resolved = [resolve(EmptyBoundary(), part) for (at, part) in assembly.members]
    outlines = [resolve(at, member.outer)
                for ((at, part), member) in zip(assembly.members, resolved)]
    regions = [PlacedRegion(placed.source, resolve(at, placed.shape))
               for ((at, part), member) in zip(assembly.members, resolved)
               for placed in member.regions]
    circle = Circle(maximum(placed -> support(placed.shape), regions))
    T = eltype(circle)
    for (k, outline) in enumerate(outlines)
        others = AbstractShape{T}[shape for (j, shape) in enumerate(outlines) if j != k]
        validate(outline, assembly, circle, others, described.outer)
    end
    rows = NamedTuple[]
    terminals = NamedTuple[]
    repeats = NamedTuple[]
    lays = NamedTuple[]
    fills = NamedTuple[]
    for (k, ((at, part), member)) in enumerate(zip(assembly.members, resolved))
        shift = length(described.regions) +
                sum(j -> length(resolved[j].regions), 1:(k - 1); init = 0)
        base = length(described.assemblies) + length(rows) + 1
        push!(rows,
            (at = at,
                regions = Int[index + shift
                              for index in eachindex(member.regions)
                              if all(row -> index ∉ row.regions, member.assemblies)],
                members = Int[index + base
                              for index in eachindex(member.assemblies)
                              if all(row -> index ∉ row.members, member.assemblies)]))
        append!(rows,
            NamedTuple[(at = at * row.at, regions = row.regions .+ shift,
                           members = row.members .+ base) for row in member.assemblies])
        append!(terminals,
            NamedTuple[(name = row.name, scope = Int[base; row.scope .+ base],
                           regions = row.regions .+ shift) for row in member.terminals])
        append!(repeats,
            NamedTuple[(placements = [at * placement for placement in row.placements],
                           copies = UnitRange{Int}[indices .+ shift
                                                   for indices in row.copies])
                       for row in member.repeats])
        append!(lays,
            NamedTuple[(
                           regions = UnitRange{Int}[indices .+ shift
                                                    for indices in row.regions],
                           radii = row.radii, pitches = row.pitches, overlengths = row.overlengths,
                           axis = at(row.axis)) for row in member.lays])
        append!(fills,
            NamedTuple[(region = row.region + shift, contents = row.contents .+ shift)
                       for row in member.fills])
    end
    return merge(described,
        (regions = [described.regions; regions], outer = circle,
            terminals = NamedTuple[described.terminals; terminals],
            assemblies = NamedTuple[described.assemblies; rows],
            repeats = NamedTuple[described.repeats; repeats],
            lays = NamedTuple[described.lays; lays], fills = NamedTuple[described.fills;
                                                                        fills]))
end

"""
$(TYPEDEF)

Lay the copies of a repeat along the cable. Each copy follows the helix about the center of
the repeat, at the lay radius of its centroid:

```math
r = \\left| \\frac{\\sum_i A_i \\mathbf c_i}{\\sum_i A_i} \\right|,
```

where ``A_i`` \\[m²\\] and ``\\mathbf c_i`` \\[m\\] are the area and centroid of each region of the
copy. The helix has the [`pitch`](@ref) and [`overlength`](@ref) at that radius, the pitch
signed by the direction of the helix. A copy at the radius 0 lies straight, with an infinite
pitch and an overlength of 1.

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

# The laid regions are the copies of the first row of repeats that the repeat adds.
function resolve(described::NamedTuple, lay::Lay)
    inner = resolve(described, lay.part)
    index = length(described.repeats) + 1
    radii = map(inner.repeats[index].copies) do indices
        shapes = [inner.regions[region].shape for region in indices]
        return hypot(centroid(shapes)...)
    end
    row = (regions = inner.repeats[index].copies, radii = radii,
        pitches = [lay.helix.dir * pitch(lay.helix, radius) for radius in radii],
        overlengths = [overlength(lay.helix, radius) for radius in radii],
        axis = (zero(eltype(radii)), zero(eltype(radii))))
    return merge(inner, (lays = NamedTuple[inner.lays; row],))
end

"""
$(TYPEDEF)

Fill a container around parts with a material. The fill region declares the container as its
primitive, either explicitly or as an offset from the part inside the enclosure. It takes the
container without the outlines of the outermost regions of its contents, those that lie inside
no other content's outline. The design refuses an empty hollow inside such an outline as a
gap.

The container is a primitive or a resolved cross-section. Its outline is the outline of the
enclosure, and `⊆` tests whether each content lies within it.

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

# A region inside an enclosure lies within its container.
function validate(shape::AbstractShape, enclosure::Enclosure, container::AbstractShape)
    shape ⊆ container || throw(DomainError(support(shape),
        "a region inside :$(enclosure.fill.tag) must lie within its container"))
    return shape
end

"""
$(TYPEDSIGNATURES)

Resolve the fill as its container minus the contents' outlines.
The holes are the outlines of the contents that lie within no other content's outline, by `⊆`.
"""
function resolve(described::NamedTuple, enclosure::Enclosure)
    container = resolve(described.outer, enclosure.fill.primitive)
    inner = resolve(described, enclosure.contents)
    contents = (length(described.regions) + 1):length(inner.regions)
    shapes = [validate(inner.regions[index].shape, enclosure, container)
              for index in contents]
    outlines = map(boundary, shapes)
    holes = [outlines[i]
             for i in eachindex(shapes)
             if !any(shapes[i] ⊆ outlines[j] for j in eachindex(shapes) if j != i)]
    placed = PlacedRegion(enclosure.fill, PerforatedShape(container, holes))
    row = (region = length(inner.regions) + 1, contents = collect(contents))
    return merge(inner,
        (regions = [inner.regions; placed], outer = boundary(container),
            fills = [inner.fills; row]))
end

function Base.map(f, region::Region)
    return Region(region.tag, map(f, region.primitive), map(f, region.material))
end
Base.map(f, stack::Stack) = Stack(map.(f, stack.items))
Base.map(f, terminal::Terminal) = Terminal(terminal.name, map(f, terminal.part))
Base.map(f, repeat::Repeat) = Repeat(map(f, repeat.pattern), map(f, repeat.part))
"""
$(TYPEDSIGNATURES)

Map `f` over each member's placement and part.
"""
function Base.map(f, assembly::Assembly)
    Assembly(Tuple(map(f, at) => map(f, part)
    for (at, part) in assembly.members))
end
Base.map(f, lay::Lay) = Lay(map(f, lay.helix), map(f, lay.part))
function Base.map(f, enclosure::Enclosure)
    return Enclosure(map(f, enclosure.fill), map(f, enclosure.contents))
end
