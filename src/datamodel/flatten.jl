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

# Returns

- A [`Terminal`](@ref) of the same name over one [`Region`](@ref), a [`Circle`](@ref) or an
  [`Annulus`](@ref).
"""
function flatten(terminal::Terminal{<:CablePart{<:Region, <:Union{Circle, Annulus}}})
    placed = terminal.part
    return Terminal(terminal.name,
        Region(placed.source.tag, placed.shape, placed.source.material))
end
function flatten(terminal::Terminal)
    strands = resolve(NamedTuple[], Placement(0, 0), 1.0, 0.0, terminal.part)
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
    return Terminal(terminal.name, Region(terminal.name, tube, material))
end
