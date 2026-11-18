"""
$(TYPEDEF)

Supertype of the declared parts of a cable: a region, or a node that layers, names, repeats,
lays, groups or encloses parts.
"""
abstract type AbstractCablePart end

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

"""
$(TYPEDEF)

Represent one region as a cable or a system places it: its declaration and its resolved
cross-section. It is the leaf of a resolved part.

$(TYPEDFIELDS)
"""
struct CablePart{R <: Region, S <: AbstractShape} <: AbstractCablePart
    "Declared region."
    source::R
    "Resolved cross-section in the frame of the part, the cable or the system that places the region."
    shape::S
    function CablePart(source::R, shape::S) where {R <: Region, S <: AbstractShape}
        return new{R, S}(source, shape)
    end
end

# A region resolves its primitive against the outline of the parts inside it. Each part adds
# the regions that it places to `regions`, in the frame that it resolves in. It returns its
# resolved part and its outline as a pair.
function resolve(outline::Union{EmptyBoundary, AbstractShape}, region::Region,
        regions::Vector{CablePart})
    shape = resolve(outline, region.primitive)
    push!(regions, CablePart(region, shape))
    return Pair{AbstractCablePart, AbstractShape}(last(regions), boundary(shape))
end

# The placed region with its cross-section in `T`: itself when the cross-section is in `T`,
# and a conversion otherwise. A perforated cross-section always converts, so that its holes are
# a vector of `AbstractShape{T}` whatever they are.
function CablePart(::Type{T},
        placed::CablePart{<:Region, <:AbstractShape{T}}) where {T <: Real}
    return placed
end
function CablePart(::Type{T}, placed::CablePart) where {T <: Real}
    return CablePart(placed.source, map(x -> convert(T, x), placed.shape))
end
function CablePart(::Type{T},
        placed::CablePart{<:Region, <:PerforatedShape{T}}) where {T <: Real}
    return CablePart(placed.source, map(x -> convert(T, x), placed.shape))
end

# A leaf of the resolved part takes the next placed region of the cable. It is a region of the
# sub-cable around it, the last row of its scope.
function resolve(rows::NamedTuple, scope::Vector{Int}, ::CablePart)
    push!(rows.regions, rows.placed[length(rows.regions) + 1])
    push!(rows.assemblies[last(scope)].regions, length(rows.regions))
    return rows
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

# Each part of a stack resolves against the outline of the parts inside it.
function resolve(outline::Union{EmptyBoundary, AbstractShape}, stack::Stack,
        regions::Vector{CablePart})
    items = AbstractCablePart[]
    for item in stack.items
        resolved = resolve(outline, item, regions)
        push!(items, resolved.first)
        outline = resolved.second
    end
    return Pair{AbstractCablePart, AbstractShape}(Stack(items), outline)
end

function resolve(rows::NamedTuple, scope::Vector{Int}, stack::Stack)
    foreach(item -> resolve(rows, scope, item), stack.items)
    return rows
end

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

function resolve(outline::Union{EmptyBoundary, AbstractShape}, terminal::Terminal,
        regions::Vector{CablePart})
    resolved = resolve(outline, terminal.part, regions)
    return Pair{AbstractCablePart, AbstractShape}(Terminal(terminal.name, resolved.first), resolved.second)
end

# The conductive regions of the part leave every terminal inside it and join the terminal of
# its name in the scope of the part. A terminal left without regions has no row. The scope of
# a row lists the sub-cables that contain its terminal, below the cable itself.
function resolve(rows::NamedTuple, scope::Vector{Int}, terminal::Terminal)
    start = length(rows.regions)
    resolve(rows, scope, terminal.part)
    conductive = Int[index
                     for index in (start + 1):length(rows.regions)
                     if rows.regions[index].source.material isa
                        Material{<:Real, :conductor}]
    own = scope[2:end]
    named = Int[index
                for row in rows.terminals
                if row.name == terminal.name && row.scope == own
                for index in row.regions]
    others = NamedTuple[(name = row.name, scope = row.scope,
                            regions = setdiff(row.regions, conductive))
                        for row in rows.terminals
                        if !(row.name == terminal.name && row.scope == own)]
    kept = NamedTuple[others;
                      (name = terminal.name, scope = own,
                          regions = sort!(union(named, conductive)))]
    append!(empty!(rows.terminals),
        sort!(NamedTuple[row for row in kept if !isempty(row.regions)];
            by = row -> first(row.regions)))
    return rows
end

"""
$(TYPEDEF)

Repeat a part at the placements of a pattern. The part is declared in its own frame, and the
pattern places each copy about the center of the enclosing part. A repeated region is a
strand, and the pattern bends a rectangular strand into an annular sector of the same area
on its course. Each copy is a scope: the terminals declared in the part are distinct
terminals in every copy. A copy of a part other than a region is a sub-cable of the design.
The outline of the copies is the circle about the center that contains them.

A resolved repeat keeps one copy in its own frame. Its pattern then lists the placements of
the copies.

$(TYPEDFIELDS)
"""
struct Repeat{P <: Union{AbstractLayoutPattern, AbstractVector{<:Placement}},
    R <: AbstractCablePart} <: AbstractCablePart
    "Pattern that places the copies, or the placements themselves in a resolved repeat."
    pattern::P
    "Part that the pattern copies."
    part::R
    function Repeat(pattern::P,
            part::R) where {
            P <: Union{AbstractLayoutPattern, AbstractVector{<:Placement}},
            R <: AbstractCablePart}
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
function resolve(outline::Union{EmptyBoundary, AbstractShape},
        repeat::Repeat{<:AbstractLayoutPattern, <:Region}, regions::Vector{CablePart})
    course = placements(resolve(outline, repeat.pattern, repeat.part.primitive),
        repeat.part.primitive)
    shapes = [validate(resolve(at, course.shape), repeat, outline)
              for at in course.placements]
    append!(regions, [CablePart(repeat.part, shape) for shape in shapes])
    return Pair{AbstractCablePart, AbstractShape}(
        Repeat(course.placements, CablePart(repeat.part, course.shape)),
        Circle(maximum(support, shapes)))
end

# The pattern places a composite part by its outline. The part resolves once, in its own
# frame, and every region of every copy lies outside the part inside the repeat.
function resolve(outline::Union{EmptyBoundary, AbstractShape},
        repeat::Repeat{<:AbstractLayoutPattern}, regions::Vector{CablePart})
    inner = CablePart[]
    resolved = resolve(EmptyBoundary(), repeat.part, inner)
    course = placements(resolve(outline, repeat.pattern, resolved.second), resolved.second)
    copies = [CablePart(placed.source,
                  validate(resolve(at, placed.shape), repeat, outline))
              for at in course.placements for placed in inner]
    append!(regions, copies)
    return Pair{AbstractCablePart, AbstractShape}(
        Repeat(course.placements, resolved.first),
        Circle(maximum(placed -> support(placed.shape), copies)))
end

# Each copy of a strand is a region of the sub-cable around the repeat.
function resolve(rows::NamedTuple, scope::Vector{Int},
        repeat::Repeat{<:AbstractVector{<:Placement}, <:CablePart})
    start = length(rows.regions)
    foreach(at -> resolve(rows, scope, repeat.part), repeat.pattern)
    push!(rows.repeats,
        (placements = repeat.pattern,
            copies = [index:index for index in (start + 1):length(rows.regions)]))
    return rows
end

# Each copy of a composite part is a sub-cable. The row of the repeat precedes those of its
# copies.
function resolve(rows::NamedTuple, scope::Vector{Int},
        repeat::Repeat{<:AbstractVector{<:Placement}})
    copies = UnitRange{Int}[]
    push!(rows.repeats, (placements = repeat.pattern, copies = copies))
    for at in repeat.pattern
        start = length(rows.regions)
        resolve(rows, scope, at => repeat.part)
        push!(copies, (start + 1):length(rows.regions))
    end
    return rows
end

# A member at its placement is a sub-cable. Its row precedes those of the sub-cables inside
# it. The terminals declared in it take that row as their innermost scope. The member lists
# its rows in its own frame. The placement then puts the sub-cables inside it in the enclosing
# frame, with the placements of the repeats and the axes of the lays.
function resolve(rows::NamedTuple, scope::Vector{Int},
        member::Pair{<:Placement, <:AbstractCablePart})
    at, part = member
    row = length(rows.assemblies) + 1
    push!(rows.assemblies[last(scope)].members, row)
    push!(rows.assemblies, (at = at, regions = Int[], members = Int[]))
    repeats, lays = length(rows.repeats), length(rows.lays)
    resolve(rows, [scope; row], part)
    for index in (row + 1):length(rows.assemblies)
        nested = rows.assemblies[index]
        rows.assemblies[index] = (
            at = at * nested.at, regions = nested.regions, members = nested.members)
    end
    for index in (repeats + 1):length(rows.repeats)
        nested = rows.repeats[index]
        rows.repeats[index] = (
            placements = [at * placement for placement in nested.placements],
            copies = nested.copies)
    end
    for index in (lays + 1):length(rows.lays)
        nested = rows.lays[index]
        rows.lays[index] = (regions = nested.regions, radii = nested.radii,
            pitches = nested.pitches, overlengths = nested.overlengths,
            axis = at(nested.axis))
    end
    return rows
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

# Each member resolves in its own frame. The members lie apart from each other and outside
# the part inside the assembly. The outline is the circle about the center that contains
# their regions.
function resolve(outline::Union{EmptyBoundary, AbstractShape}, assembly::Assembly,
        regions::Vector{CablePart})
    inner = [CablePart[] for member in assembly.members]
    resolved = [resolve(EmptyBoundary(), part, own)
                for ((at, part), own) in zip(assembly.members, inner)]
    outlines = [resolve(member.first, part.second)
                for (member, part) in zip(assembly.members, resolved)]
    copies = [CablePart(placed.source, resolve(at, placed.shape))
              for ((at, part), own) in zip(assembly.members, inner) for placed in own]
    circle = Circle(maximum(placed -> support(placed.shape), copies))
    T = eltype(circle)
    for (k, shape) in enumerate(outlines)
        others = AbstractShape{T}[other for (j, other) in enumerate(outlines) if j != k]
        validate(shape, assembly, circle, others, outline)
    end
    append!(regions, copies)
    return Pair{AbstractCablePart, AbstractShape}(
        Assembly(Tuple(member.first => part.first
        for (member, part) in zip(assembly.members, resolved))),
        circle)
end

function resolve(rows::NamedTuple, scope::Vector{Int}, assembly::Assembly)
    foreach(member -> resolve(rows, scope, member), assembly.members)
    return rows
end

"""
$(TYPEDEF)

Lay the copies of a repeat along the cable. Each copy follows the helix about the center of
the repeat, at the lay radius of its centroid:

```math
r = \\lvert \\mathbf c \\rvert,
```

where ``\\mathbf c`` is the centroid of the regions of the copy, in \\[m\\] from that center.
The helix has the [`pitch`](@ref) and [`overlength`](@ref) at that radius, the pitch signed
by the direction of the helix. A copy at the radius 0 lies straight, with an infinite pitch
and an overlength of 1.

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

function resolve(outline::Union{EmptyBoundary, AbstractShape}, lay::Lay,
        regions::Vector{CablePart})
    resolved = resolve(outline, lay.part, regions)
    return Pair{AbstractCablePart, AbstractShape}(Lay(lay.helix, resolved.first), resolved.second)
end

# The laid regions are the copies of the first row of repeats that the repeat adds. The lay
# radius of a copy is the distance of its centroid from the center of the repeat. That center
# is the origin of the sub-cable around the lay, placed in the frame of the cable.
function resolve(rows::NamedTuple, scope::Vector{Int}, lay::Lay)
    index = length(rows.repeats) + 1
    resolve(rows, scope, lay.part)
    copies = rows.repeats[index].copies
    centers = [centroid([rows.regions[region].shape for region in indices])
               for indices in copies]
    center = foldr((row, point) -> rows.assemblies[row].at(point), scope[2:end];
        init = map(zero, first(centers)))
    radii = [distance(point .- center) for point in centers]
    push!(rows.lays,
        (regions = copies, radii = radii,
            pitches = [lay.helix.dir * pitch(lay.helix, radius) for radius in radii],
            overlengths = [overlength(lay.helix, radius) for radius in radii],
            axis = (zero(eltype(radii)), zero(eltype(radii)))))
    return rows
end

"""
$(TYPEDEF)

Gather the primitives of one conductor, or of one dielectric layer, in a container and fill
the space between them with a material. The fill region declares the container as its
primitive, either explicitly or as an offset from the part inside the group. It takes the
container without the outlines of the outermost regions of its contents, those that lie
inside no other content's outline. The design refuses an empty hollow inside such an outline
as a gap.

The container is a primitive or a resolved cross-section. Its outline is the outline of the
group, and `⊆` tests whether each content lies within it.

$(TYPEDFIELDS)
"""
struct Group{F <: Union{Region, CablePart}, C <: AbstractCablePart} <:
       AbstractCablePart
    "Region of the fill material, whose primitive declares the container."
    fill::F
    "Primitives inside the container."
    contents::C
    function Group(fill::F,
            contents::C) where {F <: Union{Region, CablePart}, C <: AbstractCablePart}
        return new{F, C}(fill, contents)
    end
end

"""
$(TYPEDEF)

Contain sub-cables in a filler, such as the cores of a pipe or the cable of a duct. The fill
region declares the container as its primitive, either explicitly or as an offset from the
part inside the enclosure. It takes the container without the outlines of the outermost
regions of its contents, those that lie inside no other content's outline. The design refuses
an empty hollow inside such an outline as a gap.

The container is a primitive or a resolved cross-section. Its outline is the outline of the
enclosure, and `⊆` tests whether each content lies within it.

$(TYPEDFIELDS)
"""
struct Enclosure{F <: Union{Region, CablePart}, C <: AbstractCablePart} <:
       AbstractCablePart
    "Region of the fill material, whose primitive declares the container."
    fill::F
    "Sub-cables inside the container."
    contents::C
    function Enclosure(fill::F,
            contents::C) where {F <: Union{Region, CablePart}, C <: AbstractCablePart}
        return new{F, C}(fill, contents)
    end
end

# A region inside a group or an enclosure lies within the container of its fill.
function validate(shape::AbstractShape, fill::Region, container::AbstractShape)
    shape ⊆ container || throw(DomainError(support(shape),
        "a region inside :$(fill.tag) must lie within its container"))
    return shape
end

"""
$(TYPEDSIGNATURES)

Resolve the fill of a group or an enclosure around its contents. The fill takes its container
minus the contents' outlines. The holes are the outlines of the contents that lie within no
other content's outline, by `⊆`. The placed fill joins `regions` after its contents.

# Returns

- The pair of the resolved contents and the outline of the fill.
"""
function resolve(outline::Union{EmptyBoundary, AbstractShape}, fill::Region,
        contents::AbstractCablePart, regions::Vector{CablePart})
    container = resolve(outline, fill.primitive)
    start = length(regions)
    resolved = resolve(outline, contents, regions)
    shapes = [validate(regions[index].shape, fill, container)
              for index in (start + 1):length(regions)]
    outlines = map(boundary, shapes)
    holes = [outlines[i]
             for i in eachindex(shapes)
             if !any(shapes[i] ⊆ outlines[j] for j in eachindex(shapes) if j != i)]
    push!(regions, CablePart(fill, PerforatedShape(container, holes)))
    return Pair{AbstractCablePart, AbstractShape}(resolved.first, boundary(container))
end

function resolve(outline::Union{EmptyBoundary, AbstractShape}, group::Group,
        regions::Vector{CablePart})
    resolved = resolve(outline, group.fill, group.contents, regions)
    return Pair{AbstractCablePart, AbstractShape}(Group(last(regions), resolved.first), resolved.second)
end

function resolve(outline::Union{EmptyBoundary, AbstractShape}, enclosure::Enclosure,
        regions::Vector{CablePart})
    resolved = resolve(outline, enclosure.fill, enclosure.contents, regions)
    return Pair{AbstractCablePart, AbstractShape}(Enclosure(last(regions), resolved.first), resolved.second)
end

# The fill is the region after its contents.
function resolve(rows::NamedTuple, scope::Vector{Int}, part::Union{Group, Enclosure})
    start = length(rows.regions)
    resolve(rows, scope, part.contents)
    contents = collect((start + 1):length(rows.regions))
    resolve(rows, scope, part.fill)
    push!(rows.fills, (region = length(rows.regions), contents = contents))
    return rows
end

# The outline of a resolved part is that of its outermost region. The copies of a strand lie
# in the circle about the center that contains them.
boundary(placed::CablePart) = boundary(placed.shape)
boundary(region::Region{<:AbstractPrimitive}) = boundary(region.primitive)
boundary(stack::Stack) = boundary(last(stack.items))
boundary(terminal::Terminal) = boundary(terminal.part)
boundary(lay::Lay) = boundary(lay.part)
boundary(group::Group{<:CablePart}) = boundary(group.fill.shape)
function boundary(repeat::Repeat{<:AbstractVector{<:Placement}, <:CablePart})
    return Circle(maximum(at -> support(resolve(at, repeat.part.shape)), repeat.pattern))
end

# A resolved part begins at its nearest region, and a group where its container begins.
distance(placed::CablePart) = distance(placed.shape)
distance(stack::Stack) = minimum(distance, stack.items)
distance(terminal::Terminal) = distance(terminal.part)
distance(lay::Lay) = distance(lay.part)
distance(group::Group{<:CablePart}) = distance(group.fill.shape.outer)
function distance(repeat::Repeat{<:AbstractVector{<:Placement}, <:CablePart})
    return minimum(at -> distance(resolve(at, repeat.part.shape)), repeat.pattern)
end

# `foreach` calls `f` on each region of a resolved part, from the inside out. A repeat gives
# the region of its part once for each copy, and a group its contents before its fill.
Base.foreach(f, placed::CablePart) = f(placed)
Base.foreach(f, stack::Stack) = foreach(item -> foreach(f, item), stack.items)
Base.foreach(f, terminal::Terminal) = foreach(f, terminal.part)
Base.foreach(f, lay::Lay) = foreach(f, lay.part)
function Base.foreach(f, repeat::Repeat{<:AbstractVector{<:Placement}})
    return foreach(at -> foreach(f, repeat.part), repeat.pattern)
end
function Base.foreach(f, group::Group{<:CablePart})
    foreach(f, group.contents)
    return foreach(f, group.fill)
end

# The strands of a resolved part are its conductive regions. Each lies in the frame `at`, with
# the overlength `ℓ` and the turns per unit length `N` of the lays above it, and carries its
# geometric mean radius `g`.
function resolve(strands::Vector, at::Placement, ℓ::Real, N::Real,
        placed::CablePart{<:Region{<:Union{AbstractPrimitive, Offset},
            <:Material{<:Real, :conductor}}})
    push!(strands,
        (S = resolve(at, placed.shape), material = placed.source.material, ℓ = ℓ, N = N,
            g = gmr(placed.shape, placed.source.material.mu_r)))
    return strands
end
resolve(strands::Vector, ::Placement, ::Real, ::Real, ::CablePart) = strands

function resolve(strands::Vector, at::Placement, ℓ::Real, N::Real, stack::Stack)
    foreach(item -> resolve(strands, at, ℓ, N, item), stack.items)
    return strands
end
function resolve(strands::Vector, at::Placement, ℓ::Real, N::Real, terminal::Terminal)
    return resolve(strands, at, ℓ, N, terminal.part)
end
function resolve(strands::Vector, at::Placement, ℓ::Real, N::Real, group::Group)
    return resolve(strands, at, ℓ, N, group.contents)
end

# Each part of the stack of a terminal lists a conductive region. The dielectric layers of a
# terminal follow it in the stack around it.
function resolve(strands::Vector, at::Placement, ℓ::Real, N::Real,
        terminal::Terminal{<:Stack})
    for item in terminal.part.items
        listed = length(strands)
        resolve(strands, at, ℓ, N, item)
        validate(strands, terminal, listed)
    end
    return strands
end
function validate(strands::Vector, terminal::Terminal, listed::Int)
    length(strands) > listed || throw(ArgumentError(
        "each part of the terminal :$(terminal.name) must hold a conductive region; " *
        "declare its dielectric layers after the terminal"))
    return strands
end

# The part of a repeat lists its conductive regions once, in its own frame, and each copy
# places them.
function resolve(strands::Vector, at::Placement, ℓ::Real, N::Real,
        repeat::Repeat{<:AbstractVector{<:Placement}})
    own = resolve(NamedTuple[], Placement(0, 0), one(ℓ), zero(N), repeat.part)
    for copy in repeat.pattern, strand in own

        push!(strands,
            (S = resolve(at * copy, strand.S), material = strand.material,
                ℓ = ℓ * strand.ℓ, N = N + strand.N, g = strand.g))
    end
    return strands
end

# A laid strand follows the helix at the lay radius `r` of its centroid. The overlengths of
# nested lays multiply, and their turns add.
function resolve(strands::Vector, at::Placement, ℓ::Real, N::Real,
        lay::Lay{<:Helix, <:Repeat{<:AbstractVector{<:Placement}, <:CablePart}})
    own = resolve(NamedTuple[], Placement(0, 0), one(ℓ), zero(N), lay.part.part)
    c = centroid(lay.part.part.shape)
    for copy in lay.part.pattern, strand in own

        r = distance(copy(c))
        push!(strands,
            (S = resolve(at * copy, strand.S), material = strand.material,
                ℓ = ℓ * overlength(lay.helix, r) * strand.ℓ,
                N = N + inv(lay.helix.dir * pitch(lay.helix, r)) + strand.N, g = strand.g))
    end
    return strands
end

function Base.map(f, region::Region)
    return Region(region.tag, map(f, region.primitive), map(f, region.material))
end
function Base.map(f, placed::CablePart)
    return CablePart(map(f, placed.source), map(f, placed.shape))
end
Base.map(f, stack::Stack) = Stack(map.(f, stack.items))
Base.map(f, terminal::Terminal) = Terminal(terminal.name, map(f, terminal.part))
Base.map(f, repeat::Repeat) = Repeat(map(f, repeat.pattern), map(f, repeat.part))
function Base.map(f, repeat::Repeat{<:AbstractVector{<:Placement}})
    return Repeat([map(f, at) for at in repeat.pattern], map(f, repeat.part))
end
"""
$(TYPEDSIGNATURES)

Map `f` over each member's placement and part.
"""
function Base.map(f, assembly::Assembly)
    Assembly(Tuple(map(f, at) => map(f, part)
    for (at, part) in assembly.members))
end
Base.map(f, lay::Lay) = Lay(map(f, lay.helix), map(f, lay.part))
Base.map(f, group::Group) = Group(map(f, group.fill), map(f, group.contents))
function Base.map(f, enclosure::Enclosure)
    return Enclosure(map(f, enclosure.fill), map(f, enclosure.contents))
end
