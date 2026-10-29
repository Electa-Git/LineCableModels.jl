"""
$(TYPEDEF)

Represent one terminal of a cable as the equivalent tube that coaxial formulations read. The
material of the tube reproduces the direct-current resistance and the geometric mean radius
of the conductive regions of the terminal.

$(TYPEDFIELDS)
"""
struct BlueprintConductor{T <: Real}
    "Row of the terminal in the terminals of the design."
    terminal::Int
    "Row of the concentric assembly that holds the conductor."
    assembly::Int
    "Inner radius of the tube \\[m\\]."
    r_in::T
    "Outer radius of the tube \\[m\\]."
    r_ex::T
    "Material of the tube."
    material::Material{T, :conductor}
    function BlueprintConductor{T}(terminal::Int, assembly::Int, r_in::T, r_ex::T,
            material::Material{T, :conductor}) where {T <: Real}
        return new{T}(terminal, assembly, r_in, r_ex, material)
    end
end

"""
$(TYPEDEF)

Represent one dielectric layer of the interval of a conductor, an annulus about the center of
its assembly.

$(TYPEDFIELDS)
"""
struct BlueprintDielectric{T <: Real, M <: Material}
    "Row of the conductor whose dielectric interval holds the layer."
    conductor::Int
    "Inner radius of the layer \\[m\\]."
    r_in::T
    "Outer radius of the layer \\[m\\]."
    r_ex::T
    "Material of the layer, as the design declares it."
    material::M
    function BlueprintDielectric(conductor::Int, r_in::T, r_ex::T,
            material::M) where {T <: Real, M <: Material}
        return new{T, M}(conductor, r_in, r_ex, material)
    end
end

"""
$(TYPEDEF)

Describe a cable for coaxial formulations. Each terminal is an equivalent tube with the
dielectric layers of its interval, and the tubes lie in concentric assemblies.
[`flatten`](@ref) builds it from a [`CableDesign`](@ref).

$(TYPEDFIELDS)
"""
struct CableBlueprint{T <: Real, V <: AbstractVector{<:BlueprintDielectric}}
    "Name of the cable."
    id::String
    "Equivalent tube of each terminal, from the inside out."
    conductors::Vector{BlueprintConductor{T}}
    "Dielectric layers, from the inside out."
    dielectrics::V
    "Range of the dielectric layers in the interval of each conductor."
    dielectric_ranges::Vector{UnitRange{Int}}
    """
    Relative permeability of the interval of each conductor, with the correction for the
    helical path of its current \\[dimensionless\\].
    """
    mu_r_ins::Vector{T}
    "Range of the conductors of each concentric assembly."
    assembly_ranges::Vector{UnitRange{Int}}
    "Center of each concentric assembly in the frame of the cable, in \\[m\\]."
    positions::Vector{Tuple{T, T}}
    function CableBlueprint{T}(id::String, conductors::Vector{BlueprintConductor{T}},
            dielectrics::V, dielectric_ranges::Vector{UnitRange{Int}}, mu_r_ins::Vector{T},
            assembly_ranges::Vector{UnitRange{Int}}, positions::Vector{Tuple{T, T}}
    ) where {T <: Real, V <: AbstractVector{<:BlueprintDielectric}}
        return new{T, V}(id, conductors, dielectrics, dielectric_ranges, mu_r_ins,
            assembly_ranges, positions)
    end
end

"""
$(TYPEDSIGNATURES)

Reduce a cable design to its [`CableBlueprint`](@ref). Each terminal becomes a tube about the
center of the cable. The tube spans the conductive regions of the terminal and the fills
between them, from their nearest to their farthest extent. The regions outside the tube, up to
the next terminal, become the layers of its dielectric interval. Regions inside the first
terminal lie in no interval.

The conductive regions ``i`` of a terminal carry direct current in parallel. A region of
resistivity ``\\rho_i``, cross-section ``A_i`` and overlength ``\\ell_i`` along its helices has
the conductance ``G_i = A_i / (\\rho_i \\ell_i)`` per unit length. Its share of the current is
``w_i = G_i / \\sum_j G_j``. The terminal has the resistance ``R = 1 / \\sum_i G_i`` and the
temperature coefficient ``\\alpha = \\sum_i w_i \\alpha_i``. Its turns per unit length are
``N = \\sum_i w_i N_i``, where ``N_i`` sums the inverse pitches of the helices of region ``i``,
signed by their direction. Its geometric mean radius ``g`` follows from

```math
\\log g = \\sum_i w_i^2 \\log g_i + 2 \\sum_{i < j} w_i w_j \\log d_{ij}.
```

An annulus about the center of inner radius ``a`` and outer radius ``b`` has the self term

```math
\\log g_i = \\log b - \\mu_i \\left[\\frac{a^4}{(b^2 - a^2)^2} \\log\\frac{b}{a}
    - \\frac{3a^2 - b^2}{4 (b^2 - a^2)}\\right],
```

which is ``\\log b - \\mu_i / 4`` for a circle and ``\\log b`` for a thin shell,
``a = b``. Any other region has the self term
``\\log \\sqrt{A_i / \\pi} - \\mu_i / 4`` of the circle of its area. The mutual distance
``d_{ij}`` is the distance between the centroids. For a region inside the hole of an annulus, it
is the mean distance of Rosa (1908), Eq. (57):

```math
\\log d_{ij} = \\frac{b^2 \\log b - a^2 \\log a}{b^2 - a^2} - \\frac{1}{2}.
```

The tube of inner radius ``r_i`` and outer radius ``r_e`` takes the resistivity
``R \\pi (r_e^2 - r_i^2)``. Its relative permeability gives the self term of an annulus the
value ``\\log g``:

```math
\\mu_\\mathrm{tube} = \\frac{\\log r_e - \\log g}
    {\\frac{r_i^4}{(r_e^2 - r_i^2)^2} \\log\\frac{r_e}{r_i}
        - \\frac{3 r_i^2 - r_e^2}{4 (r_e^2 - r_i^2)}},
```

with the first term of the denominator zero for ``r_i = 0``. The interval of the tube
reaches from ``r_e`` to ``r_d``, and its layer ``k`` of relative permeability ``\\mu_k``
reaches from ``a_k`` to ``b_k``. The interval has

```math
\\mu_{ins} = \\frac{\\sum_k \\mu_k \\log(b_k / a_k)}{\\log(r_d / r_e)}
    \\left[1 + \\frac{2 \\pi^2 N^2 (r_d^2 - r_e^2)}{\\log(r_d / r_e)}\\right],
```

and ``\\mu_{ins} = 1`` without layers.

# Returns

- A [`CableBlueprint`](@ref) with one concentric assembly, the cable.

# Errors

- Throws `ArgumentError` when a terminal lies in a copy of a repeat, or when a fill
  surrounds the regions of several terminals or of none.
- Throws `ArgumentError` when the regions of a terminal are not contiguous from the inside out
  or do not tile an annulus about the center, or when a dielectric layer is not such an annulus.
- Throws `ArgumentError` when the conductive regions of one terminal have different reference
  temperatures.
"""
function flatten(design::CableDesign)
    T = eltype(design)
    regions = design.regions
    overlengths = ones(T, length(regions))
    turns = zeros(T, length(regions))
    for row in design.lays,
        (member, indices) in enumerate(design.repeats[row.repeat].copies),
        index in indices
        overlengths[index] *= row.overlengths[member]
        turns[index] += inv(row.pitches[member])
    end
    assigned = [index for row in design.terminals for index in row.regions]
    blocks = [sort!([row.regions;
                     [fill.region
                      for fill in design.fills
                      if !isempty(intersect(fill.contents, assigned)) &&
                         issubset(intersect(fill.contents, assigned), row.regions)]])
              for row in design.terminals]
    validate(design, flatten, blocks)
    stops = [[first(block) - 1 for block in blocks[2:end]]; length(regions)]
    reduced = map(eachindex(design.terminals)) do k
        members = design.terminals[k].regions
        shapes = [regions[index].shape for index in members]
        materials = [regions[index].source.material for index in members]
        conductances = [area(shapes[i]) / (materials[i].rho * overlengths[members[i]])
                        for i in eachindex(members)]
        shares = conductances ./ sum(conductances)
        rings = [isapprox(area(shape), π * (support(shape)^2 - distance(shape)^2))
                 for shape in shapes]
        own = map(eachindex(shapes)) do i
            a, b, μ = distance(shapes[i]), support(shapes[i]), materials[i].mu_r
            rings[i] || return log(area(shapes[i]) / π) / 2 - μ / 4
            iszero(a) && return log(b) - μ / 4
            isapprox(a, b) && return log(b)
            return log(b) -
                   μ * (a^4 / (b^2 - a^2)^2 * log(b / a) - (3a^2 - b^2) / (4(b^2 - a^2)))
        end
        cross = sum(
            (begin
                 a, b = distance(shapes[i]), support(shapes[i])
                 c, d = distance(shapes[j]), support(shapes[j])
                 mean = if rings[i] && a > zero(a) &&
                           support(shapes[j]) <= a + sqrt(eps(a)) * max(one(a), a)
                     (b^2 * log(b) - a^2 * log(a)) / (b^2 - a^2) - one(a) / 2
                 elseif rings[j] && c > zero(c) &&
                        support(shapes[i]) <= c + sqrt(eps(c)) * max(one(c), c)
                     (d^2 * log(d) - c^2 * log(c)) / (d^2 - c^2) - one(c) / 2
                 else
                     log(hypot((centroid(shapes[i]) .- centroid(shapes[j]))...))
                 end
                 shares[i] * shares[j] * mean
             end
            for i in eachindex(shapes) for j in (i + 1):length(shapes));
            init = zero(eltype(shares)))
        gmr = sum(shares .^ 2 .* own) + 2cross
        contents = [index for fill in design.fills if fill.region in blocks[k]
                    for index in fill.contents]
        tube = [regions[index].shape for index in blocks[k] if index ∉ contents]
        inner = convert(T, minimum(distance, tube))
        outer = convert(T, maximum(support, tube))
        term = iszero(inner) ? zero(outer) :
               inner^4 / (outer^2 - inner^2)^2 * log(outer / inner)
        mu = (log(outer) - gmr) / (term - (3inner^2 - outer^2) / (4(outer^2 - inner^2)))
        material = Material(Val(:conductor),
            convert(T, π * (outer^2 - inner^2) / sum(conductances)), one(T), convert(T, mu),
            convert(T, first(materials).T0),
            convert(T, sum(shares .* [item.alpha for item in materials])), zero(T),
            convert(T, 90), zero(T), zero(T))
        interval = [BlueprintDielectric(k, convert(T, distance(regions[index].shape)),
                        convert(T, support(regions[index].shape)),
                        regions[index].source.material)
                    for index in (last(blocks[k]) + 1):stops[k]]
        N = sum(shares .* [turns[index] for index in members])
        relative = if isempty(interval)
            one(T)
        else
            width = log(last(interval).r_ex / outer)
            sum(layer -> layer.material.mu_r * log(layer.r_ex / layer.r_in), interval) /
            width * (1 + 2π^2 * N^2 * (last(interval).r_ex^2 - outer^2) / width)
        end
        return (conductor = BlueprintConductor{T}(k, 1, inner, outer, material),
            interval = interval, relative = convert(T, relative))
    end
    layers = [layer for item in reduced for layer in item.interval]
    counts = cumsum([0; [length(item.interval) for item in reduced]])
    at = first(design.assemblies).at
    return CableBlueprint{T}(
        design.id, BlueprintConductor{T}[item.conductor for item in reduced],
        collect(Union{map(typeof, layers)...}, layers),
        [(counts[k] + 1):counts[k + 1] for k in eachindex(reduced)],
        T[item.relative for item in reduced], [1:length(reduced)], [(at.x, at.y)])
end

# The design reduces to coaxial tubes. Each terminal lies in the scope of the cable. Its
# regions and the fills between them are contiguous in stack order and tile an annulus about
# the center. Its materials share one reference temperature. Past the first tube, every
# region outside the tubes is an annulus about the center.
function validate(design::CableDesign, ::typeof(flatten),
        blocks::AbstractVector{<:AbstractVector{Int}})
    for row in design.terminals
        isempty(row.scope) || throw(ArgumentError(
            "flatten reduces the terminals about the center of a cable; the terminal " *
            ":$(row.name) of the cable $(design.id) lies in a copy of a repeat"))
    end
    tubes = reduce(vcat, blocks; init = Int[])
    for row in design.fills
        row.region in tubes || throw(ArgumentError(
            "flatten reduces the fills between the regions of one terminal; the fill " *
            ":$(design.regions[row.region].source.tag) of the cable $(design.id) surrounds " *
            "the regions of several terminals or of none"))
    end
    for (row, block) in zip(design.terminals, blocks)
        block == first(block):last(block) || throw(ArgumentError(
            "flatten reduces a terminal whose regions are contiguous from the inside out; " *
            "other regions lie between those of :$(row.name) in the cable $(design.id)"))
        shapes = [design.regions[index].shape for index in block]
        inner, outer = minimum(distance, shapes), maximum(support, shapes)
        isapprox(sum(area, shapes), π * (outer^2 - inner^2)) || throw(ArgumentError(
            "flatten reduces a terminal whose regions tile an annulus about the center; " *
            ":$(row.name) of the cable $(design.id) does not"))
        temperature = design.regions[first(row.regions)].source.material.T0
        all(index -> isapprox(design.regions[index].source.material.T0, temperature),
            row.regions) || throw(ArgumentError(
            "the materials of the terminal :$(row.name) of the cable $(design.id) must " *
            "share one reference temperature"))
    end
    for index in eachindex(design.regions)
        (index in tubes || all(block -> first(block) > index, blocks)) && continue
        shape = design.regions[index].shape
        isapprox(area(shape), π * (support(shape)^2 - distance(shape)^2)) ||
            throw(ArgumentError(
                "flatten reduces dielectric layers that are annuli about the center; " *
                ":$(design.regions[index].source.tag) of the cable $(design.id) is not"))
    end
    return design
end
