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
circle or an annulus inside a round outline, and any other outline itself. The region takes
the resistivity ``R A(S)`` and the relative permeability that gives ``S`` the radius ``g``:

```math
\\mu = \\frac{\\log \\operatorname{gmr}(S, 0) - \\log g}
    {\\log \\operatorname{gmr}(S, 0) - \\log \\operatorname{gmr}(S, 1)}.
```

Its turns per unit length are ``N = \\sum_i w_i N_i``, where ``N_i`` sums the inverse
pitches of the helices of region ``i``, signed by their direction.

# Returns

- The pair of a [`Terminal`](@ref) of the same name over one [`Region`](@ref) and the turns
  ``N`` \\[1/m\\].

# Errors

- Throws `ArgumentError` when a part of the stack of the terminal lists only dielectric
  regions.
"""
function flatten(inner::Union{EmptyBoundary, AbstractShape}, terminal::Terminal)
    strands = resolve(NamedTuple[], Placement(0, 0), 1.0, 0.0, terminal)
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
    material = Material(Val(:conductor), area(tube) / sum(G), 1,
        gap / (b - log(gmr(tube, 1))), first(strands).material.T0,
        sum(w .* [strand.material.alpha for strand in strands]), 0, 90, 0, 0)
    return Terminal(terminal.name,
        Region(terminal.name, tube, material)) => sum(w .* [strand.N for strand in strands])
end

"""
$(TYPEDSIGNATURES)

Return a terminal over one region as its own flattened form, with no turns. A circle or an
annulus is the region of its resolved cross-section. Any other region is its declaration.

# Returns

- The pair of the [`Terminal`](@ref) over one [`Region`](@ref) and zero turns \\[1/m\\].
"""
function flatten(terminal::Terminal{<:CablePart{<:Region, <:Union{Circle, Annulus}}})
    placed = terminal.part
    return Terminal(terminal.name,
        Region(placed.source.tag, placed.shape, placed.source.material)) => zero(eltype(placed.shape))
end
function flatten(terminal::Terminal{<:CablePart})
    return Terminal(terminal.name, terminal.part.source) => zero(eltype(terminal.part.shape))
end
flatten(terminal::Terminal{<:Region}) = terminal => zero(eltype(terminal.part.material))

"""
$(TYPEDSIGNATURES)

Reduce a stack to one region per terminal and per dielectric layer, from the inside out.
Each terminal becomes one region, as [`flatten`](@ref) of a terminal gives it. A round
conductor is a tube. The fill of its group beyond its metal is a layer of its own, where it
has a width. A group around a terminal is the conductor of that terminal.

Any other group is one layer of the composite material of its regions, by their fractions
of the area. A dielectric circle or annulus is its own layer, and any other region is its
declaration. An enclosure and an assembly keep their form around their flattened members. A
flattened stack is its own flattened form.

The interval of a round terminal begins at its tube, of outer radius ``r_e``. It ends at
the next terminal or the end of the stack, at the radius ``r_d``. The current of the terminal has
``N`` turns per unit length. Each layer of the interval takes its relative permeability
times

```math
1 + \\frac{2 \\pi^2 N^2 (r_d^2 - r_e^2)}{\\log(r_d / r_e)},
```

evaluated as ``1 + 4 \\pi^2 N^2 r_d r_e \\sinh(u) / u`` with ``u = \\log(r_d / r_e)``.

# Returns

- A [`Stack`](@ref) of regions, of terminals over one region each, and of the enclosures
  and assemblies of the stack.

# Errors

- Throws `ArgumentError` when a part of the stack is neither a terminal, a group, an
  enclosure, an assembly nor a region.
- Throws `ArgumentError` when two parts of the stack declare one terminal.
"""
function flatten(stack::Stack)
    items = AbstractCablePart[]
    tubes = Tuple{Int, Real}[]
    foreach(item -> flatten(items, tubes, item), stack.items)
    for (k, (index, N)) in enumerate(tubes)
        stop = k < length(tubes) ? first(tubes[k + 1]) - 1 : length(items)
        stop > index && !iszero(nominal(N)) || continue
        r_e, r_d = boundary(items[index]).r, boundary(items[stop]).r
        u = log(r_d / r_e)
        for i in (index + 1):stop
            medium = items[i].material
            items[i] = Region(items[i].tag, items[i].primitive,
                Material(medium, medium, medium.rho, medium.eps_r,
                    medium.mu_r * (1 + 4π^2 * N^2 * r_d * r_e * sinh(u) / u), medium.T0,
                    medium.alpha, medium.rho_thermal, medium.theta_max, medium.tan_delta,
                    medium.sigma_solar))
        end
    end
    return Stack(items)
end

# A terminal adds its one region to the flattened parts `items`, outside the last of them,
# then the fill beyond its metal. `tubes` lists the place of each such region in `items`
# with its turns. A terminal over one region does not read the outline inside it.
function flatten(items::Vector{AbstractCablePart}, tubes::Vector, terminal::Terminal)
    validate(terminal, flatten, items, tubes)
    flat = flatten(isempty(items) ? EmptyBoundary() : boundary(last(items)), terminal)
    push!(items, flat.first)
    push!(tubes, (length(items), flat.second))
    return flatten(items, flat.first, terminal.part)
end
function flatten(items::Vector{AbstractCablePart}, tubes::Vector,
        terminal::Terminal{<:Union{Region, CablePart}})
    flat = flatten(validate(terminal, flatten, items, tubes))
    push!(items, flat.first)
    push!(tubes, (length(items), flat.second))
    return items
end

# A stack declares the conductor of a terminal in one part.
function validate(terminal::Terminal, ::typeof(flatten), items::Vector{AbstractCablePart},
        tubes::Vector)
    any(tube -> items[first(tube)].name == terminal.name, tubes) && throw(ArgumentError(
        "the terminal :$(terminal.name) is declared in two parts of its stack; " *
        "flatten reads the conductor of a terminal from one Terminal"))
    return terminal
end

# A stack inside a stack adds its parts in their order.
function flatten(items::Vector{AbstractCablePart}, tubes::Vector, stack::Stack)
    foreach(item -> flatten(items, tubes, item), stack.items)
    return items
end

# A group around a terminal is the conductor of that terminal, which begins where the
# container of the group does. Any other group is a dielectric layer.
function flatten(items::Vector{AbstractCablePart}, tubes::Vector,
        group::Group{<:CablePart, <:Terminal})
    return flatten(items, tubes, Terminal(group.contents.name, group))
end
function flatten(items::Vector{AbstractCablePart}, tubes::Vector,
        group::Group{<:CablePart, <:Stack})
    k = findfirst(item -> item isa Terminal, group.contents.items)
    k === nothing && return flatten(items, group.fill.shape.outer, group)
    return flatten(items, tubes, Terminal(group.contents.items[k].name, group))
end
function flatten(items::Vector{AbstractCablePart}, ::Vector, group::Group{<:CablePart})
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
function flatten(items::Vector{AbstractCablePart}, ::Vector,
        placed::CablePart{<:Region, <:Union{Circle, Annulus}})
    push!(items, Region(placed.source.tag, placed.shape, placed.source.material))
    return items
end
function flatten(items::Vector{AbstractCablePart}, ::Vector, placed::CablePart)
    push!(items, placed.source)
    return items
end
function flatten(items::Vector{AbstractCablePart}, ::Vector, region::Region)
    push!(items, region)
    return items
end

# An enclosure and an assembly keep their form around their flattened members.
function flatten(items::Vector{AbstractCablePart}, ::Vector,
        part::Union{Enclosure, Assembly})
    push!(items, flatten(part))
    return items
end

function flatten(::Vector{AbstractCablePart}, ::Vector, part::AbstractCablePart)
    throw(ArgumentError(
        "flatten forms a layer of a stack from a group or a region; " *
        "a $(nameof(typeof(part))) of the stack is neither"))
end

# The fill of a round group reaches from the tube `core` of its terminal to its container.
# It is a layer where it has a width. A stack ends at its last part, and any other part
# ends at its metal.
function flatten(items::Vector{AbstractCablePart}, core::Terminal,
        group::Group{<:CablePart{
            <:Region, <:PerforatedShape{<:Real, <:Union{Circle, Annulus}}}})
    container = group.fill.shape.outer
    r_e, r_d = r_ex(core.part.primitive), r_ex(container)
    nominal(r_d - r_e) > tolerance(r_d) && push!(items,
        Region(group.fill.source.tag, Annulus(r_e, r_d, container.at),
            group.fill.source.material))
    return items
end
function flatten(items::Vector{AbstractCablePart}, core::Terminal, stack::Stack)
    return flatten(items, core, last(stack.items))
end
flatten(items::Vector{AbstractCablePart}, ::Terminal, ::AbstractCablePart) = items

"""
$(TYPEDSIGNATURES)

Flatten the contents of an enclosure. The enclosure keeps its form around them, and the
region of its fill takes the container as its primitive.

# Returns

- An [`Enclosure`](@ref) whose fill is a [`Region`](@ref) of the container.
"""
function flatten(enclosure::Enclosure{<:CablePart})
    fill = enclosure.fill
    return Enclosure(Region(fill.source.tag, fill.shape.outer, fill.source.material),
        flatten(enclosure.contents))
end
function flatten(enclosure::Enclosure{<:Region})
    return Enclosure(enclosure.fill, flatten(enclosure.contents))
end

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
