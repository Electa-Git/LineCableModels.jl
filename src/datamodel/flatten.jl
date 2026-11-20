"""
$(TYPEDSIGNATURES)

Return the cross-section of the one region that stands for a conductor. The conductor has
the outline `O`, and `inner` is the outline of the part inside it in its stack. Its metal has
the cross-sections `S`.

Inside a circle, the cross-section is a tube about the center of `O`. It reaches the
farthest extent of the metal from that center. It begins where the part inside it ends, and
it is a circle when there is none.

Any other outline moves inward until it meets the metal. The cross-section is the region
inside its parallel at the nearest distance from the metal to the outline. That distance is
zero for metal that touches the outline.

# Errors

- Throws `DomainError` where the inward parallel of the outline fails, as
  `parallel(contour, d)` refuses it. It fails for a nearest distance at or above the radius
  of a rounded corner of the outline.
"""
function flatten(O::Circle, ::EmptyBoundary, S::AbstractVector{<:AbstractShape})
    return Circle(maximum(shape -> support(resolve(inv(O.at), shape)), S), O.at)
end
function flatten(O::Circle, inner::Circle, S::AbstractVector{<:AbstractShape})
    return Annulus(inner.r, maximum(shape -> support(resolve(inv(O.at), shape)), S), O.at)
end
function flatten(O::AbstractPrimitive, ::EmptyBoundary, S::AbstractVector{<:AbstractShape})
    outline = Contour(O)
    s = minimum(clearance(a, b)[3]
    for shape in S for a in outline.outer for b in Contour(shape).outer)
    return ArcPolygon(parallel(outline, -s).outer)
end

"""
$(TYPEDSIGNATURES)

Reduce a terminal to one region with the material that stands for its metal, under the name
of the terminal. `inner` is the outline of the part inside the terminal in its stack, and an
`EmptyBoundary` when there is none. `O` is the outline of the conductor of the terminal, and
`strands` lists its conductive regions.

The conductive regions ``i`` of the terminal carry direct current in parallel. A region of
resistivity ``\\rho_i``, cross-section ``S_i`` and overlength ``\\ell_i`` along its helices
has the conductance ``G_i = A(S_i) / (\\rho_i \\ell_i)`` per unit length. Its share of the
current is ``w_i = G_i / \\sum_j G_j``. The terminal has the resistance
``R = 1 / \\sum_i G_i`` and the temperature coefficient ``\\alpha = \\sum_i w_i \\alpha_i``.
Its geometric mean radius ``g`` follows from

```math
\\log g = \\sum_i w_i^2 \\log g_i + 2 \\sum_{i < j} w_i w_j \\log \\operatorname{gmd}(S_i,S_j),
```

where ``g_i`` is the [`gmr`](@ref) of region ``i``. The cross-section ``S`` of the one
region is what [`flatten`](@ref) gives the outline `O` outside `inner`. It is a
circle or an annulus inside a round outline, and any other outline moved inward until it
meets the metal. The region is
the [`Region`](@ref) of the cross-section ``S`` with the resistance ``R`` and the geometric
mean radius ``g``. The method sums the differences
``\\log \\operatorname{gmr}(S, 0) - \\log g_i`` and
``\\log \\operatorname{gmr}(S, 0) - \\log \\operatorname{gmd}(S_i,S_j)``, because ``g`` lies
close to ``\\operatorname{gmr}(S, 0)`` for a thin tube.

The current of the terminal makes ``N = \\sum_i w_i N_i`` turns per unit length. ``N_i``
sums the inverse pitches of the helices that region ``i`` follows about the axis of the
conductor, signed by their direction. A helix about another axis only lengthens the region.
The strand lay of a rope member off the axis of the rope is such a helix.

A terminal with turns has a [`Lay`](@ref) over its region. The helix has the lay length
``1 / \\lvert N \\rvert`` and the direction of ``N``.

# Returns

- A [`Terminal`](@ref) of the same name over one [`Region`](@ref), or over a [`Lay`](@ref)
  of that region when its current has turns.

# Notes

- The one region stands for metals of one relative permeability. The method warns when the
  conductive regions of the terminal have several. The region then keeps their shares of the
  direct current. It leaves out the flux of the current of one metal in another, and the
  exchange of current between the metals above DC. Magnetic metal with a dielectric between
  it and the other metal is a terminal of its own, on the same connection.

# Errors

- Throws `ArgumentError` when the terminal does not contain a conductive region.
- Throws `ArgumentError` when the conductive regions of the terminal have different
  reference temperatures.
"""
function flatten(inner::Union{EmptyBoundary, AbstractShape}, terminal::Terminal,
        @nospecialize(O::AbstractShape), strands::Vector)
    validate(strands, flatten, terminal)
    S = [strand.S for strand in strands]
    media = validate([strand.material for strand in strands], flatten, terminal)
    return flatten(terminal.name, flatten(O, inner, S), S, media,
        [strand.ℓ for strand in strands], [strand.g for strand in strands],
        [strand.N for strand in strands])
end

# The sums of the association run on the columns of the conductive regions, each a vector
# of one type. `S` holds their cross-sections and `media` their materials. `ℓ` holds their
# overlengths and `g` their geometric mean radii. `N` holds their turns.
function flatten(name::Symbol, tube::AbstractPrimitive, S::Vector{<:AbstractShape},
        media::Vector{<:Material}, ℓ::Vector{<:Real}, g::Vector{<:Real},
        N::Vector{<:Real})
    G = [area(S[i]) / (media[i].rho * ℓ[i]) for i in eachindex(S)]
    w = G ./ sum(G)
    b = log(gmr(tube, 0))
    gap = sum(w[i]^2 * (b - log(g[i])) for i in eachindex(S)) +
          2sum(
        (w[i] * w[j] * (b - log(gmd(S[i], S[j])))
        for i in eachindex(S) for j in (i + 1):length(S));
        init = zero(eltype(w)))
    region = Region(name, tube, 1 / sum(G), exp(b - gap), first(media).T0,
        sum(w .* [medium.alpha for medium in media]))
    turns = sum(w .* N)
    iszero(nominal(turns)) && return Terminal(name, region)
    return Terminal(name,
        Lay(Helix(Pitch(1 / abs(turns)); dir = nominal(turns) > 0 ? 1 : -1), region))
end

# A terminal contains a conductive region.
function validate(strands::Vector, ::typeof(flatten), terminal::Terminal)
    isempty(strands) && throw(ArgumentError(
        "the terminal :$(terminal.name) does not contain a conductive region"))
    return strands
end

# The materials of the conductive regions of a terminal share one reference temperature.
# One region stands for metals of one relative permeability. With several, it keeps their
# shares of the direct current. It leaves out the flux of the current of one metal in
# another, and the exchange of current between them above DC.
function validate(media::Vector{<:Material}, ::typeof(flatten), terminal::Terminal)
    all(medium -> nominal(medium.T0) == nominal(first(media).T0), media) ||
        throw(ArgumentError(
            "the materials of the terminal :$(terminal.name) must share one reference temperature"))
    all(medium -> nominal(medium.mu_r) == nominal(first(media).mu_r), media) ||
        @warn "the terminal :$(terminal.name) merges metals of different relative permeability into one region; the region keeps their shares of the direct current and leaves out the flux of the current of one metal in another"
    return media
end

"""
$(TYPEDSIGNATURES)

Return a terminal over one region as its own flattened form. A circle or an annulus is the
region of its resolved cross-section. A layer that an offset declares around another
outline than a circle stays the resolved part. An offset declares a thickness and not a
cross-section. Any other region is its declaration. A lay over the region stays around it.

# Returns

- The [`Terminal`](@ref) over one [`Region`](@ref), over the [`Lay`](@ref) of one, or over
  the resolved [`CablePart`](@ref) of such a layer.
"""
function flatten(terminal::Terminal{<:CablePart{<:Region, <:Union{Circle, Annulus}}})
    placed = terminal.part
    return Terminal(terminal.name,
        Region(placed.source.tag, placed.shape, placed.source.material))
end
flatten(terminal::Terminal{<:CablePart{<:Region{<:Offset}, <:PerforatedShape}}) = terminal
flatten(terminal::Terminal{<:CablePart}) = Terminal(terminal.name, terminal.part.source)
flatten(terminal::Terminal{<:Region}) = terminal
function flatten(terminal::Terminal{<:Lay{
        <:Helix, <:CablePart{<:Region, <:Union{Circle, Annulus}}}})
    placed = terminal.part.part
    return Terminal(terminal.name,
        Lay(terminal.part.helix,
            Region(placed.source.tag, placed.shape, placed.source.material)))
end
flatten(terminal::Terminal{<:Lay{<:Helix, <:Region}}) = terminal

"""
$(TYPEDSIGNATURES)

Reduce a stack to one region per terminal and per dielectric layer, from the inside out.
Each terminal becomes one region, as [`flatten`](@ref) of a terminal gives it. A round
conductor is a tube, and the outline of any other conductor moves inward to its metal. The
fill of its group beyond it is a layer of its own. The layer has no width where the metal
ends at the container of the group. A group around a terminal is the conductor of that
terminal.

The dielectric parts of a terminal, before its conductive parts or after them, are layers
of the stack. Any other group is one layer of the composite material of its regions, by
their fractions of the area. A dielectric circle or annulus is its own layer, and any other
region is its declaration. An enclosure and an assembly keep their form around their
flattened members. A flattened stack is its own flattened form.

# Returns

- A [`Stack`](@ref) of regions, of terminals over one region each or over its lay, and of
  the enclosures and assemblies of the stack.

# Errors

- Throws `ArgumentError` when a part of the stack is neither a terminal, a group, an
  enclosure, an assembly nor a region.
- Throws `ArgumentError` when two parts of the stack declare one terminal.
- Throws `ArgumentError` when a terminal does not contain a conductive region, or when a
  dielectric part lies between its conductive parts.
- Throws `ArgumentError` when two terminals of the stack are in contact, with no dielectric
  of some width between them.
"""
function flatten(stack::Stack)
    items = AbstractCablePart[]
    foreach(item -> flatten(items, item), stack.items)
    return validate(Stack(items), flatten)
end

# A dielectric of some width lies between each terminal of a flattened stack and the next
# one. `gap` counts the parts from the terminal `i` to the next one.
function validate(stack::Stack, ::typeof(flatten))
    items = stack.items
    i = findfirst(item -> item isa Terminal, items)
    i === nothing && return stack
    while true
        gap = findfirst(item -> item isa Terminal, view(items, (i + 1):lastindex(items)))
        gap === nothing && break
        all(isempty, view(items, (i + 1):(i + gap - 1))) && throw(ArgumentError(
            "the terminals :$(items[i].name) and :$(items[i + gap].name) are in contact with no dielectric between them and cannot be two terminals. A dielectric between two terminals is a region of the design between them. Metal of one phase in contact is one terminal."))
        i += gap
    end
    return stack
end

"""
$(TYPEDSIGNATURES)

Reduce a cable design to the flattened stack of its resolved parts.

# Returns

- A [`Stack`](@ref), as [`flatten`](@ref) of a stack gives it.

# Errors

- Throws the errors of [`flatten`](@ref) of a stack.
"""
flatten(design::CableDesign) = flatten(design.tree)

# A terminal adds its one region to the flattened parts `items`, outside the last of them,
# then the fill beyond its metal. A terminal over one region does not read the outline
# inside it.
function flatten(items::Vector{AbstractCablePart}, terminal::Terminal)
    validate(terminal, flatten, items)
    flat = flatten(isempty(items) ? EmptyBoundary() : boundary(last(items)), terminal,
        boundary(terminal.part), resolve(NamedTuple[], Placement(0, 0), 1.0, 0.0, terminal))
    push!(items, flat)
    return flatten(items, flat, terminal.part)
end

# The conductive parts of the stack of a terminal are its conductor. Its dielectric parts
# are layers of the flattened stack: `before` of them lie inside the conductor and `after`
# of them around it. A stack inside the stack gives its parts in their order.
function flatten(items::Vector{AbstractCablePart}, terminal::Terminal{<:Stack})
    validate(terminal, flatten, items)
    parts = terminal.part.items
    k = findfirst(part -> part isa Stack, parts)
    k === nothing || return flatten(items,
        Terminal(terminal.name,
            Stack(AbstractCablePart[view(parts, 1:(k - 1)); parts[k].items;
                                    view(parts, (k + 1):lastindex(parts))])))
    strands = NamedTuple[]
    before = after = 0
    for part in parts
        listed = length(strands)
        resolve(strands, Placement(0, 0), 1.0, 0.0, part)
        if length(strands) > listed
            validate(terminal, flatten, after)
        elseif isempty(strands)
            before += 1
        else
            after += 1
        end
    end
    foreach(part -> flatten(items, part), Iterators.take(parts, before))
    flatten(items, terminal, parts[end - after], strands)
    foreach(part -> flatten(items, part), Iterators.drop(parts, length(parts) - after))
    return items
end

# The conductive parts of a terminal follow one another. `after` counts the dielectric parts
# that follow the conductive parts listed so far.
function validate(terminal::Terminal, ::typeof(flatten), after::Int)
    after == 0 || throw(ArgumentError(
        "a dielectric part lies between the conductive parts of the terminal " *
        ":$(terminal.name); two conductors with a dielectric between them are two terminals"))
    return terminal
end

# The conductor of the stack of a terminal adds its one region to `items`, then the fill
# beyond its metal. `conductor` is the outermost conductive part of the stack, and `strands`
# lists the conductive regions of the terminal.
function flatten(items::Vector{AbstractCablePart}, @nospecialize(terminal::Terminal),
        @nospecialize(conductor::AbstractCablePart), strands::Vector)
    flat = flatten(isempty(items) ? EmptyBoundary() : boundary(last(items)), terminal,
        boundary(conductor), strands)
    push!(items, flat)
    return flatten(items, flat, conductor)
end

# A layer that an offset declares around another outline than a circle is the whole
# conductor of a terminal with no other metal. So is a region under the lay of its turns.
# The terminal over such a part alone is its own flattened form, with the layers of the
# stack around it. Other metal of the terminal inside the part merges into its outline.
function flatten(items::Vector{AbstractCablePart}, terminal::Terminal,
        conductor::Union{CablePart{<:Region{<:Offset}, <:PerforatedShape},
            Lay{<:Helix, <:CablePart}},
        strands::Vector)
    length(strands) == 1 &&
        return push!(items, flatten(Terminal(terminal.name, conductor)))
    push!(items,
        flatten(isempty(items) ? EmptyBoundary() : boundary(last(items)), terminal,
            boundary(conductor), strands))
    return items
end

function flatten(items::Vector{AbstractCablePart},
        terminal::Terminal{<:Union{
            Region, CablePart, Lay{<:Helix, <:Union{Region, CablePart}}}})
    push!(items, flatten(validate(terminal, flatten, items)))
    return items
end

# A stack declares the conductor of a terminal in one part.
function validate(
        terminal::Terminal, ::typeof(flatten), items::Vector{AbstractCablePart})
    terminal.name in [item.name for item in items if item isa Terminal] &&
        throw(ArgumentError(
            "the terminal :$(terminal.name) is declared in two parts of its stack; " *
            "flatten reads the conductor of a terminal from one Terminal"))
    return terminal
end

# A stack inside a stack adds its parts in their order.
function flatten(items::Vector{AbstractCablePart}, stack::Stack)
    foreach(item -> flatten(items, item), stack.items)
    return items
end

# A group around a terminal is the conductor of that terminal, which begins where the
# container of the group does. Any other group is a dielectric layer.
function flatten(
        items::Vector{AbstractCablePart}, group::Group{<:CablePart, <:Terminal})
    return flatten(items, Terminal(group.contents.name, group))
end
function flatten(items::Vector{AbstractCablePart}, group::Group{<:CablePart, <:Stack})
    k = findfirst(item -> item isa Terminal, group.contents.items)
    k === nothing && return flatten(items, group.fill.shape.outer, group)
    return flatten(items, Terminal(group.contents.items[k].name, group))
end
function flatten(items::Vector{AbstractCablePart}, group::Group{<:CablePart})
    return flatten(items, group.fill.shape.outer, group)
end

# The regions of a dielectric group tile its container, a circle or an annulus. They are
# one layer of their composite material, each with its fraction of the area.
function flatten(items::Vector{AbstractCablePart}, container::Union{Circle, Annulus},
        group::Group)
    pieces = CablePart[]
    foreach(placed -> push!(pieces, placed), group)
    A = [area(placed.shape) for placed in pieces]
    push!(items,
        Region(group.fill.source.tag,
            Annulus(r_in(container), r_ex(container), container.at),
            Material([placed.source.material for placed in pieces], A ./ sum(A))))
    return items
end

# A dielectric circle or annulus is its own layer. Any other region is its declaration, and
# so is a flattened region.
function flatten(items::Vector{AbstractCablePart},
        placed::CablePart{<:Region, <:Union{Circle, Annulus}})
    push!(items, Region(placed.source.tag, placed.shape, placed.source.material))
    return items
end
function flatten(items::Vector{AbstractCablePart}, placed::CablePart)
    push!(items, placed.source)
    return items
end
function flatten(items::Vector{AbstractCablePart}, region::Region)
    push!(items, region)
    return items
end

# An enclosure and an assembly keep their form around their flattened members.
function flatten(items::Vector{AbstractCablePart}, part::Union{Enclosure, Assembly})
    push!(items, flatten(part))
    return items
end

function flatten(::Vector{AbstractCablePart}, part::AbstractCablePart)
    throw(ArgumentError(
        "flatten forms a layer of a stack from a group or a region; " *
        "a $(nameof(typeof(part))) of the stack is neither"))
end

# The fill of a round group reaches from the tube `core` of its terminal to its container.
# It is a layer of the stack. Where the metal ends at the container the layer has no width,
# and its two radii have their own derivatives. A stack ends at its last part, and any other
# part ends at its metal.
function flatten(items::Vector{AbstractCablePart}, core::Terminal,
        group::Group{<:CablePart{
            <:Region, <:PerforatedShape{<:Real, <:Union{Circle, Annulus}}}})
    container = group.fill.shape.outer
    push!(items,
        Region(group.fill.source.tag,
            Annulus(r_ex(boundary(core)), r_ex(container), container.at),
            group.fill.source.material))
    return items
end

# The fill of any other group reaches from the outline of its conductor to its container. It
# is the layer of the nearest distance between the two, which an offset declares.
function flatten(items::Vector{AbstractCablePart}, core::Terminal,
        group::Group{<:CablePart{<:Region, <:PerforatedShape}})
    container = Contour(group.fill.shape.outer)
    conductor = Contour(boundary(core))
    s = minimum(clearance(a, b)[3] for a in container.outer for b in conductor.outer)
    push!(items, Region(group.fill.source.tag, Offset(s), group.fill.source.material))
    return items
end
function flatten(items::Vector{AbstractCablePart}, core::Terminal, stack::Stack)
    return flatten(items, core, last(stack.items))
end
flatten(items::Vector{AbstractCablePart}, ::Terminal, ::AbstractCablePart) = items

"""
$(TYPEDSIGNATURES)

Flatten the contents of an enclosure. The enclosure keeps its form and its fill around
them. The fill of a resolved enclosure is its resolved part, whose cross-section is the
container without the outlines of the contents.

# Returns

- An [`Enclosure`](@ref) of the same fill around the flattened contents.
"""
flatten(enclosure::Enclosure) = Enclosure(enclosure.fill, flatten(enclosure.contents))

"""
$(TYPEDSIGNATURES)

Flatten each member of an assembly once, in its own frame, as a stack. The members keep
their placements.

# Returns

- An [`Assembly`](@ref) of the placements and the flattened stacks.
"""
function flatten(assembly::Assembly)
    return Assembly(Tuple(at => flatten(Stack(part)) for (at, part) in assembly.members))
end

"""
$(TYPEDSIGNATURES)

Flatten the part of a repeat once, in its own frame, as a stack. The copies keep their
placements.

# Returns

- A [`Repeat`](@ref) of the placements and the flattened stack.
"""
function flatten(repeat::Repeat{<:AbstractVector{<:Placement}})
    return Repeat(repeat.pattern, flatten(Stack(repeat.part)))
end
