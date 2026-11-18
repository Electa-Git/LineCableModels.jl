"""
$(TYPEDSIGNATURES)

Return the cross-section of the equivalent tube of a conductor, a circle or an annulus about
the center. The conductor has the outline `O` and begins at the radius `r_i` \\[m\\]. Its
metal has the cross-sections `S`.

Inside a circle, the tube reaches from `r_i` to the farthest extent of the metal. Any other
outline, of perimeter ``l`` \\[m\\], takes Ametani's cylinder of the area of the metal:

```math
r_e = \\frac{l}{2\\pi}, \\qquad r_i = \\sqrt{r_e^2 - \\frac{1}{\\pi} \\sum_k A(S_k)}.
```

An inner radius within `tolerance(r_e)` of zero, decided on `nominal`, is zero. A tube
without an inner radius is a circle.
"""
function flatten(O::Circle, r_i::Real, S::AbstractVector{<:AbstractShape})
    r_e = maximum(support, S)
    return iszero(nominal(r_i)) ? Circle(r_e) : Annulus(r_i, r_e)
end
function flatten(O::AbstractShape, ::Real, S::AbstractVector{<:AbstractShape})
    r_e = perimeter(O) / 2π
    return nominal(r_e^2 - sum(area, S) / π) <= tolerance(r_e)^2 ? Circle(r_e) :
           Annulus(sqrt(r_e^2 - sum(area, S) / π), r_e)
end

"""
$(TYPEDSIGNATURES)

Reduce a resolved terminal to one region: the equivalent tube of its metal, under the name
of the terminal. A terminal over one circle or annulus is its own equivalent, and keeps its
region and its material.

The conductive regions ``i`` of any other terminal carry direct current in parallel. A region of
resistivity ``\\rho_i``, cross-section ``S_i`` and overlength ``\\ell_i`` along its helices
has the conductance ``G_i = A(S_i) / (\\rho_i \\ell_i)`` per unit length. Its share of the
current is ``w_i = G_i / \\sum_j G_j``. The terminal has the resistance
``R = 1 / \\sum_i G_i`` and the temperature coefficient ``\\alpha = \\sum_i w_i \\alpha_i``.
Its geometric mean radius ``g`` follows from

```math
\\log g = \\sum_i w_i^2 \\log g_i + 2 \\sum_{i < j} w_i w_j \\log \\operatorname{gmd}(S_i,S_j),
```

where ``g_i`` is the [`gmr`](@ref) of region ``i``. The tube is the equivalent of the outline
of the terminal, of inner radius ``r_i`` and outer radius ``r_e``. It takes the resistivity
``R \\pi (r_e^2 - r_i^2)`` and the relative permeability that gives it the radius ``g``:

```math
\\mu = \\frac{\\log r_e - \\log g}
    {\\log r_e - \\log \\operatorname{gmd}(S_\\mathrm{tube},S_\\mathrm{tube})}.
```

Its turns per unit length are ``N = \\sum_i w_i N_i``, where ``N_i`` sums the inverse
pitches of the helices of region ``i``, signed by their direction.

# Returns

- The pair of a [`Terminal`](@ref) of the same name over one [`Region`](@ref), a
  [`Circle`](@ref) or an [`Annulus`](@ref), and the turns ``N`` \\[1/m\\].

# Errors

- Throws `ArgumentError` when a part of the stack of the terminal lists only dielectric
  regions.
"""
function flatten(terminal::Terminal{<:CablePart{<:Region, <:Union{Circle, Annulus}}})
    placed = terminal.part
    return Terminal(terminal.name,
        Region(placed.source.tag, placed.shape, placed.source.material)) => zero(eltype(placed.shape))
end
function flatten(terminal::Terminal)
    strands = resolve(NamedTuple[], Placement(0, 0), 1.0, 0.0, terminal)
    S = [strand.S for strand in strands]
    G = [area(strand.S) / (strand.material.rho * strand.ℓ) for strand in strands]
    w = G ./ sum(G)
    tube = flatten(boundary(terminal.part), distance(terminal.part), S)
    r_e = r_ex(tube)
    gap = sum(w[i]^2 * (log(r_e) - log(strands[i].g)) for i in eachindex(strands)) +
          2sum(
        (w[i] * w[j] * (log(r_e) - log(gmd(S[i], S[j])))
        for i in eachindex(S) for j in (i + 1):length(S));
        init = zero(eltype(w)))
    material = Material(Val(:conductor), area(tube) / sum(G), 1,
        gap / (log(r_e) - log(gmd(tube, tube))), first(strands).material.T0,
        sum(w .* [strand.material.alpha for strand in strands]), 0, 90, 0, 0)
    return Terminal(terminal.name,
        Region(terminal.name, tube, material)) => sum(w .* [strand.N for strand in strands])
end

"""
$(TYPEDSIGNATURES)

Reduce a resolved stack to circles and annuli, from the inside out. Each terminal becomes
its equivalent tube, as [`flatten`](@ref) of a terminal gives it. The fill of a round group
beyond its metal is a layer of its own, where it has a width. A group around a terminal is
the conductor of that terminal.

Any other group is one layer of the composite material of its regions, by their fractions
of the area. A dielectric circle or annulus is its own layer. An offset around any other
outline is the annulus of its declared thickness outside the part inside it.

The interval of a terminal reaches from its tube, of outer radius ``r_e``, to the next
terminal or the end of the stack, at the radius ``r_d``. The current of the terminal has
``N`` turns per unit length. Each layer of the interval takes its relative permeability
times

```math
1 + \\frac{2 \\pi^2 N^2 (r_d^2 - r_e^2)}{\\log(r_d / r_e)},
```

evaluated as ``1 + 4 \\pi^2 N^2 r_d r_e \\sinh(u) / u`` with ``u = \\log(r_d / r_e)``.

# Returns

- A [`Stack`](@ref) of regions and of terminals over one region each.

# Errors

- Throws `ArgumentError` when a part of the stack is neither a terminal, nor a group, nor a
  dielectric circle, annulus or offset.
- Throws `ArgumentError` when two parts of the stack declare one terminal.
"""
function flatten(stack::Stack)
    items = AbstractCablePart[]
    tubes = Tuple{Int, Real}[]
    foreach(item -> flatten(items, tubes, item), stack.items)
    for (k, (index, N)) in enumerate(tubes)
        stop = k < length(tubes) ? first(tubes[k + 1]) - 1 : length(items)
        stop > index && !iszero(nominal(N)) || continue
        r_e, r_d = r_ex(boundary(items[index])), r_ex(boundary(items[stop]))
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

# A terminal adds its tube to the flattened parts `items`, then the fill beyond its metal.
# `tubes` lists the place of each tube in `items` with its turns.
function flatten(items::Vector{AbstractCablePart}, tubes::Vector, terminal::Terminal)
    flat = flatten(validate(terminal, flatten, items, tubes))
    push!(items, flat.first)
    push!(tubes, (length(items), flat.second))
    return flatten(items, r_ex(boundary(flat.first)), terminal.part)
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
        Region(group.fill.source.tag, Annulus(r_in(container), r_ex(container)),
            Material([placed.source.material for placed in pieces], A ./ sum(A))))
    return items
end

# A dielectric circle or annulus is its own layer.
function flatten(items::Vector{AbstractCablePart}, ::Vector,
        placed::CablePart{<:Region, <:Union{Circle, Annulus}})
    push!(items, Region(placed.source.tag, placed.shape, placed.source.material))
    return items
end

# An offset around any other outline is the annulus of its declared thickness outside the
# part inside it.
function flatten(items::Vector{AbstractCablePart}, ::Vector,
        placed::CablePart{<:Region{<:Offset}, <:PerforatedShape})
    r = r_ex(boundary(last(items)))
    push!(items,
        Region(placed.source.tag, Annulus(r, r + placed.source.primitive.t),
            placed.source.material))
    return items
end

function flatten(::Vector{AbstractCablePart}, ::Vector, part::AbstractCablePart)
    throw(ArgumentError(
        "flatten forms a layer of a stack from a group or a dielectric circle, annulus or offset; " *
        "a $(nameof(typeof(part))) of the stack is none of them"))
end

# The fill of a round group reaches from the metal of its terminal, at the radius `r_e`, to
# its container. It is a layer where it has a width. A stack ends at its last part, and any
# other part ends at its metal.
function flatten(items::Vector{AbstractCablePart}, r_e::Real,
        group::Group{<:CablePart{
            <:Region, <:PerforatedShape{<:Real, <:Union{Circle, Annulus}}}})
    r_d = r_ex(group.fill.shape.outer)
    nominal(r_d - r_e) > tolerance(r_d) && push!(items,
        Region(group.fill.source.tag, Annulus(r_e, r_d), group.fill.source.material))
    return items
end
function flatten(items::Vector{AbstractCablePart}, r_e::Real, stack::Stack)
    return flatten(items, r_e, last(stack.items))
end
flatten(items::Vector{AbstractCablePart}, ::Real, ::AbstractCablePart) = items
