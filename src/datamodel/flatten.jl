"""
$(TYPEDSIGNATURES)

Return the cross-section of the one region that stands for a conductor. The conductor has
the outline `O`, and `inner` is the outline of the part inside it in its stack. Its metal has
the cross-sections `S`.

Inside a circle, the cross-section is a tube about the center of `O`. It reaches the
farthest extent of the metal from that center. It begins where the part inside it ends, and
it is a circle when there is none. Any other outline is the cross-section itself.
"""
function flatten(O::Circle, ::EmptyBoundary, S::AbstractVector{<:AbstractShape})
    return Circle(maximum(shape -> support(resolve(inv(O.at), shape)), S), O.at)
end
function flatten(O::Circle, inner::Circle, S::AbstractVector{<:AbstractShape})
    return Annulus(inner.r, maximum(shape -> support(resolve(inv(O.at), shape)), S), O.at)
end
flatten(O::AbstractPrimitive, ::EmptyBoundary, ::AbstractVector{<:AbstractShape}) = O

"""
$(TYPEDSIGNATURES)

Reduce a terminal to one region with the material that stands for its metal, under the name
of the terminal. `inner` is the outline of the part inside the terminal in its stack, and an
`EmptyBoundary` when there is none.

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
region is what [`flatten`](@ref) gives the outline of the terminal outside `inner`. It is a
circle or an annulus inside a round outline, and any other outline itself. The region is
the [`Region`](@ref) of the cross-section ``S`` with the resistance ``R`` and the geometric
mean radius ``g``. The method sums the differences
``\\log \\operatorname{gmr}(S, 0) - \\log g_i`` and
``\\log \\operatorname{gmr}(S, 0) - \\log \\operatorname{gmd}(S_i,S_j)``, because ``g`` lies
close to ``\\operatorname{gmr}(S, 0)`` for a thin tube.

The current of the terminal makes ``N = \\sum_i w_i N_i`` turns per unit length, where
``N_i`` sums the inverse pitches of the helices of region ``i``, signed by their direction.
A terminal with turns has a [`Lay`](@ref) over its region. The helix has the lay ratio
``1 / (2 r_e \\lvert N \\rvert)`` and the direction of ``N``, where ``r_e`` is the outer
radius of the region. Its pitch at that radius is ``1 / \\lvert N \\rvert``.

# Returns

- A [`Terminal`](@ref) of the same name over one [`Region`](@ref), or over a [`Lay`](@ref)
  of that region when its current has turns.

# Errors

- Throws `ArgumentError` when a part of the stack of the terminal lists only dielectric
  regions.
- Throws `ArgumentError` when the conductive regions of the terminal have different
  reference temperatures.
"""
function flatten(inner::Union{EmptyBoundary, AbstractShape}, terminal::Terminal)
    strands = validate(
        resolve(NamedTuple[], Placement(0, 0), 1.0, 0.0, terminal), flatten, terminal)
    S = [strand.S for strand in strands]
    G = [area(strand.S) / (strand.material.rho * strand.ℓ) for strand in strands]
    w = G ./ sum(G)
    tube = flatten(boundary(terminal.part), inner, S)
    b = log(gmr(tube, 0))
    gap = sum(w[i]^2 * (b - log(strands[i].g)) for i in eachindex(strands)) +
          2sum(
        (w[i] * w[j] * (b - log(gmd(S[i], S[j])))
        for i in eachindex(S) for j in (i + 1):length(S));
        init = zero(eltype(w)))
    region = Region(terminal.name, tube, 1 / sum(G), exp(b - gap),
        first(strands).material.T0,
        sum(w .* [strand.material.alpha for strand in strands]))
    N = sum(w .* [strand.N for strand in strands])
    iszero(nominal(N)) && return Terminal(terminal.name, region)
    return Terminal(terminal.name,
        Lay(Helix(LayRatio(1 / (2 * r_ex(tube) * abs(N))); dir = nominal(N) > 0 ? 1 : -1),
            region))
end

# The conductive regions of a terminal share one reference temperature.
function validate(strands::Vector, ::typeof(flatten), terminal::Terminal)
    all(
        strand -> nominal(strand.material.T0) == nominal(first(strands).material.T0),
        strands) || throw(ArgumentError(
        "the materials of the terminal :$(terminal.name) must share one reference temperature"))
    return strands
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
conductor is a tube. The fill of its group beyond its metal is a layer of its own. The
layer has no width where the metal ends at the container of the group. A group around a
terminal is the conductor of that terminal.

Any other group is one layer of the composite material of its regions, by their fractions
of the area. A dielectric circle or annulus is its own layer, and any other region is its
declaration. An enclosure and an assembly keep their form around their flattened members. A
flattened stack is its own flattened form.

# Returns

- A [`Stack`](@ref) of regions, of terminals over one region each or over its lay, and of
  the enclosures and assemblies of the stack.

# Errors

- Throws `ArgumentError` when a part of the stack is neither a terminal, a group, an
  enclosure, an assembly nor a region.
- Throws `ArgumentError` when two parts of the stack declare one terminal.
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
    flat = flatten(isempty(items) ? EmptyBoundary() : boundary(last(items)), terminal)
    push!(items, flat)
    return flatten(items, flat, terminal.part)
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
