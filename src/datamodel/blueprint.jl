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

Represent one dielectric layer of a cable, about the axis of its assembly. A layer of the
interval of a conductor is an annulus. The first layer of the interior of an assembly, its
filler, surrounds its members, and the others are annuli outside it.

$(TYPEDFIELDS)
"""
struct BlueprintDielectric{T <: Real, M <: Material}
    "Inner radius of the layer \\[m\\]."
    r_in::T
    "Outer radius of the layer \\[m\\]."
    r_ex::T
    "Material of the layer, as the design declares it."
    material::M
    function BlueprintDielectric(r_in::T, r_ex::T,
            material::M) where {T <: Real, M <: Material}
        return new{T, M}(r_in, r_ex, material)
    end
end

"""
    BlueprintDielectric(r_in, r_ex, media, bands)

Build the layer from `r_in` to `r_ex` of the dielectric regions of materials `media` and
cross-sections `bands` as one composite material. Each region weighs by its area fraction
``f_k``:

```math
\\frac{1}{\\rho} = \\sum_k \\frac{f_k}{\\rho_k}, \\qquad
\\varepsilon_r = \\sum_k f_k \\varepsilon_k, \\qquad
\\varepsilon_r \\tan\\delta = \\sum_k f_k \\varepsilon_k \\tan\\delta_k.
```

The other properties are those of the last region. The layer is an insulator when every
region is one, and a semiconductor otherwise. One region keeps its material.
"""
function BlueprintDielectric(r_in::T, r_ex::T, media::AbstractVector{<:Material},
        bands::AbstractVector{<:AbstractShape}) where {T <: Real}
    fractions = area.(bands) ./ sum(area, bands)
    permittivity = sum(fractions .* [item.eps_r for item in media])
    base = last(media)
    values = (
        inv(sum(fractions ./ [item.rho for item in media])), permittivity,
        base.mu_r, base.T0, base.alpha, base.rho_thermal, base.theta_max,
        sum(fractions .* [item.eps_r * item.tan_delta for item in media]) /
        permittivity, base.sigma_solar)
    return BlueprintDielectric(r_in, r_ex,
        foldl((merged, medium) -> Material(merged, medium, values...), media))
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
    """
    Range of the dielectric layers in the interval of each conductor. A conductor without
    one has the empty range `(n + 1):n`, where `n` counts the layers before its place.
    """
    dielectric_ranges::Vector{UnitRange{Int}}
    """
    Relative permeability of the interval of each conductor, with the correction for the
    helical path of its current \\[dimensionless\\].
    """
    mu_r_ins::Vector{T}
    """
    Range of the conductors of each concentric assembly. An assembly without one has the
    empty range `(n + 1):n`, where `n` counts the conductors before its place.
    """
    assembly_ranges::Vector{UnitRange{Int}}
    """
    Axis that [`flatten`](@ref) gives each concentric assembly in the frame of the cable, in
    \\[m\\].
    """
    positions::Vector{Tuple{T, T}}
    """
    Range of the interior layers of each concentric assembly: its filler and the annuli
    outside it. An assembly without members has the empty range `(n + 1):n`, where `n`
    counts the layers before its place.
    """
    interior_ranges::Vector{UnitRange{Int}}
    """
    Rows of the members of each assembly. A sub-cable without a terminal gives its members
    to its container.
    """
    members::Vector{Vector{Int}}
    function CableBlueprint{T}(id::String, conductors::Vector{BlueprintConductor{T}},
            dielectrics::V, dielectric_ranges::Vector{UnitRange{Int}}, mu_r_ins::Vector{T},
            assembly_ranges::Vector{UnitRange{Int}}, positions::Vector{Tuple{T, T}},
            interior_ranges::Vector{UnitRange{Int}}, members::Vector{Vector{Int}}
    ) where {T <: Real, V <: AbstractVector{<:BlueprintDielectric}}
        return new{T, V}(id, conductors, dielectrics, dielectric_ranges, mu_r_ins,
            assembly_ranges, positions, interior_ranges, members)
    end
end

"""
$(TYPEDSIGNATURES)

Return whether `shape`, or the set `shapes`, fills an annulus about the center. Its area
is then that of the annulus between its least `distance` and its greatest `support`. The
comparison takes `nominal` values and the relative tolerance of `isapprox`. A centered
circle or annulus fills one, and a course of wires does not.
"""
function isannular(shapes::AbstractVector{<:AbstractShape})
    return isapprox(nominal(sum(area, shapes)),
        nominal(π * (maximum(support, shapes)^2 - minimum(distance, shapes)^2)))
end
isannular(shape::AbstractShape) = isannular([shape])

"""
$(TYPEDSIGNATURES)

Reduce a cable design to its [`CableBlueprint`](@ref). Each row of the assemblies of the
design is a concentric assembly about its axis. Each terminal becomes a tube about the axis
of its assembly. The cable lies about its origin. A member lies about the centroid of the
regions of its first terminal, or about its placement when it has no terminal. The tube
spans the conductive regions of the terminal and the containers of the fills between them,
from their nearest to their farthest extent.

The regions of the same row outside the tube, up to its next terminal or its last region,
become the layers of its dielectric interval. The regions of a row around its members form
the interior of the assembly. The annuli about the axis outside its other regions are its
outer layers, from the inside out. The other regions form its filler, one layer from their
nearest to their farthest extent. A sub-cable without a terminal gives its members to its
container and its regions to the filler of its container. Regions inside the first terminal
of a row without members lie in no interval.

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

A block that tiles an annulus takes that tube. Any other block takes the tube
of its metal area ``S`` and of the perimeter ``l`` of its outline, Ametani's equivalent:

```math
r_e = \\frac{l}{2\\pi}, \\qquad r_i = \\sqrt{r_e^2 - \\frac{S}{\\pi}}.
```

Each layer of its interval then lies around the previous radius, with the distance across
the metal of its band as its thickness.

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

Several dielectric regions that tile an annulus, such as an open tape and its fill, form one
layer. Under a radial field they act in parallel. The layer takes their composite material
by their area fractions, the rule of [`BlueprintDielectric`](@ref). The filler of an
interior homogenizes its regions by the same area fractions.

# Returns

- A [`CableBlueprint`](@ref) with one concentric assembly per row of the assemblies of the
  design, the root first.

# Errors

- Throws `ArgumentError` when a fill surrounds the regions of several terminals of its own
  row, or of none.
- Throws `ArgumentError` when the regions of a terminal are not contiguous from the inside
  out.
- Throws `ArgumentError` when a layer of an interval is no annulus about the axis of its
  assembly and tiles none.
- Throws `ArgumentError` when the conductive regions of one terminal have different reference
  temperatures.
"""
function flatten(design::CableDesign)
    T = eltype(design)
    regions = design.regions
    overlengths = ones(T, length(regions))
    turns = zeros(T, length(regions))
    for row in design.lays,
        (member, indices) in enumerate(row.regions),
        index in indices
        overlengths[index] *= row.overlengths[member]
        turns[index] += inv(row.pitches[member])
    end
    assigned = [index for row in design.terminals for index in row.regions]
    owners = [isempty(row.scope) ? 1 : last(row.scope) for row in design.terminals]
    blocks = [sort!([row.regions;
                     [fill.region
                      for fill in design.fills
                      if fill.region in design.assemblies[owners[k]].regions &&
                         !isempty(intersect(fill.contents, assigned)) &&
                         issubset(intersect(fill.contents, assigned), row.regions)]])
              for (k, row) in enumerate(design.terminals)]
    positions = map(enumerate(design.assemblies)) do (k, row)
        (k == 1 || k ∉ owners) && return (row.at.x, row.at.y)
        terminal = design.terminals[findfirst(==(k), owners)]
        return centroid([regions[index].shape for index in terminal.regions])
    end
    centered = [[k == 1 ? placed.shape :
                 resolve(inv(Placement(positions[k]...)), placed.shape)
                 for placed in regions] for k in eachindex(design.assemblies)]
    validate(design, flatten, blocks, centered)
    stops = [maximum([index
                      for index in design.assemblies[owners[k]].regions
                      if all(
                 j -> owners[j] != owners[k] || j <= k || index < first(blocks[j]),
                 eachindex(blocks))]) for k in eachindex(blocks)]
    reduced = map(eachindex(design.terminals)) do k
        members = design.terminals[k].regions
        shapes = [centered[owners[k]][index] for index in members]
        materials = [regions[index].source.material for index in members]
        conductances = [area(shapes[i]) / (materials[i].rho * overlengths[members[i]])
                        for i in eachindex(members)]
        shares = conductances ./ sum(conductances)
        rings = [isannular(shape) for shape in shapes]
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
                           support(shapes[j]) <= a + tolerance(a)
                     (b^2 * log(b) - a^2 * log(a)) / (b^2 - a^2) - one(a) / 2
                 elseif rings[j] && c > zero(c) &&
                        support(shapes[i]) <= c + tolerance(c)
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
        containers = [fill.region for fill in design.fills if fill.region in blocks[k]]
        tube = [index in containers ? centered[owners[k]][index].outer :
                centered[owners[k]][index]
                for index in blocks[k] if index ∉ contents]
        annular = isannular([centered[owners[k]][index] for index in blocks[k]])
        outer = convert(T,
            annular ? maximum(support, tube) :
            perimeter(boundary(centered[owners[k]][last(blocks[k])])) / 2π)
        inner = convert(T, annular ? minimum(distance, tube) :
                           sqrt(outer^2 - sum(area, shapes) / π))
        term = iszero(inner) ? zero(outer) :
               inner^4 / (outer^2 - inner^2)^2 * log(outer / inner)
        mu = (log(outer) - gmr) / (term - (3inner^2 - outer^2) / (4(outer^2 - inner^2)))
        material = Material(Val(:conductor),
            convert(T, π * (outer^2 - inner^2) / sum(conductances)), one(T), convert(T, mu),
            convert(T, first(materials).T0),
            convert(T, sum(shares .* [item.alpha for item in materials])), zero(T),
            convert(T, 90), zero(T), zero(T))
        indices = [index
                   for index in design.assemblies[owners[k]].regions
                   if last(blocks[k]) < index <= stops[k]]
        grouped = [index for fill in design.fills if fill.region in indices
                   for index in fill.contents]
        interval = BlueprintDielectric[]
        for index in indices
            index in grouped && continue
            tiled = [fill for fill in design.fills if fill.region == index]
            pieces = isempty(tiled) ? [index] : [only(tiled).contents; index]
            bands = [centered[owners[k]][i] for i in pieces]
            media = [regions[i].source.material for i in pieces]
            lower = convert(T,
                annular ? minimum(distance, bands) :
                (isempty(interval) ? outer : last(interval).r_ex))
            upper = convert(
                T, annular ? maximum(support, bands) :
                   lower + dimension(centered[owners[k]][index]))
            push!(interval, BlueprintDielectric(lower, upper, media, bands))
        end
        N = sum(shares .* [turns[index] for index in members])
        relative = if isempty(interval)
            one(T)
        else
            width = log(last(interval).r_ex / outer)
            sum(layer -> layer.material.mu_r * log(layer.r_ex / layer.r_in), interval) /
            width * (1 + 2π^2 * N^2 * (last(interval).r_ex^2 - outer^2) / width)
        end
        return (conductor = BlueprintConductor{T}(k, owners[k], inner, outer, material),
            interval = interval, relative = convert(T, relative))
    end
    hosts = collect(eachindex(design.assemblies))
    contained = [Int[] for _ in design.assemblies]
    for (k, row) in enumerate(design.assemblies)
        for member in row.members
            if member in owners
                push!(contained[hosts[k]], member)
            else
                hosts[member] = hosts[k]
            end
        end
    end
    heads = [minimum([first(blocks[j]) for j in eachindex(blocks) if owners[j] == k];
                 init = length(regions) + 1) for k in eachindex(design.assemblies)]
    interiors = map(eachindex(design.assemblies)) do k
        shapes = centered[k]
        indices = isempty(contained[k]) ? Int[] :
                  sort!([index
                         for m in eachindex(design.assemblies) if hosts[m] == k
                         for index in design.assemblies[m].regions if index < heads[m]])
        tilings = [fill
                   for fill in design.fills
                   if fill.region in indices && all(in(indices), fill.contents)]
        grouped = [index for fill in tilings for index in fill.contents]
        bands = map([head for head in indices if head ∉ grouped]) do head
            tiled = [fill for fill in tilings if fill.region == head]
            return isempty(tiled) ? [head] : [only(tiled).contents; head]
        end
        annuli = Vector{Int}[]
        while length(bands) > 1
            outside = shapes[last(bands)]
            reach = nominal(maximum(support, shapes[reduce(vcat, bands[1:(end - 1)])]))
            nearest = nominal(minimum(distance, outside))
            isannular(outside) && (nearest >= reach || isapprox(nearest, reach)) || break
            pushfirst!(annuli, pop!(bands))
        end
        isempty(bands) && return BlueprintDielectric[]
        return BlueprintDielectric[BlueprintDielectric(
                                       convert(T, minimum(distance, shapes[band])),
                                       convert(T, maximum(support, shapes[band])),
                                       [placed.source.material for placed in regions[band]],
                                       shapes[band])
                                   for band in [[sort!(reduce(vcat, bands))]; annuli]]
    end
    ends = [maximum(row.regions; init = 0) for row in design.assemblies]
    for k in reverse(eachindex(design.assemblies))
        ends[k] = maximum([ends[k]; ends[design.assemblies[k].members]])
    end
    order = sortperm(eachindex(design.assemblies); by = k -> (ends[k], -k))
    layers = BlueprintDielectric[]
    dielectric_ranges = Vector{UnitRange{Int}}(undef, length(reduced))
    interior_ranges = Vector{UnitRange{Int}}(undef, length(design.assemblies))
    assembly_ranges = Vector{UnitRange{Int}}(undef, length(design.assemblies))
    count = 0
    for k in order
        start = length(layers) + 1
        append!(layers, interiors[k])
        interior_ranges[k] = start:length(layers)
        terminals = findall(==(k), owners)
        assembly_ranges[k] = (count + 1):(count + length(terminals))
        count += length(terminals)
        for j in terminals
            offset = length(layers)
            append!(layers, reduced[j].interval)
            dielectric_ranges[j] = (offset + 1):length(layers)
        end
    end
    return CableBlueprint{T}(
        design.id, BlueprintConductor{T}[item.conductor for item in reduced],
        collect(Union{map(typeof, layers)...}, layers), dielectric_ranges,
        T[item.relative for item in reduced], assembly_ranges, positions, interior_ranges,
        contained)
end

function validate(design::CableDesign, ::typeof(flatten),
        blocks::AbstractVector{<:AbstractVector{Int}},
        centered::AbstractVector{<:AbstractVector{<:AbstractShape}})
    owners = [isempty(row.scope) ? 1 : last(row.scope) for row in design.terminals]
    tubes = reduce(vcat, blocks; init = Int[])
    for row in design.fills
        own = only([k
                    for (k, assembly) in enumerate(design.assemblies)
                    if row.region in assembly.regions])
        row.region in tubes ||
            (!isempty(design.assemblies[own].members) &&
             all(index -> index ∉ design.assemblies[own].regions, row.contents)) ||
            (all(
                index -> index ∉ [i for terminal in design.terminals
                          for i in terminal.regions],
                row.contents) &&
             any(k -> owners[k] == own && last(blocks[k]) < row.region, eachindex(blocks))) ||
            throw(ArgumentError(
                "flatten reduces the fills between the regions of one terminal; the fill " *
                ":$(design.regions[row.region].source.tag) of the cable $(design.id) surrounds " *
                "the regions of several terminals or of none"))
    end
    for (k, (row, block)) in enumerate(zip(design.terminals, blocks))
        block == [index
         for index in design.assemblies[owners[k]].regions
         if first(block) <= index <= last(block)] || throw(ArgumentError(
            "flatten reduces a terminal whose regions are contiguous from the inside out; " *
            "other regions lie between those of :$(row.name) in the cable $(design.id)"))
        temperature = design.regions[first(row.regions)].source.material.T0
        all(index -> isapprox(design.regions[index].source.material.T0, temperature),
            row.regions) || throw(ArgumentError(
            "the materials of the terminal :$(row.name) of the cable $(design.id) must " *
            "share one reference temperature"))
    end
    for (k, assembly) in enumerate(design.assemblies)
        shapes = centered[k]
        for index in assembly.regions
            (index in tubes || all(j -> owners[j] != k || first(blocks[j]) > index,
                eachindex(blocks))) && continue
            previous = last([j
                             for j in eachindex(blocks)
                             if owners[j] == k && last(blocks[j]) < index])
            metal = [shapes[i] for i in blocks[previous]]
            isannular(metal) || continue
            grouped = [fill
                       for fill in design.fills
                       if (fill.region == index || index in fill.contents) &&
                          all(i -> i ∉ tubes, fill.contents)]
            bands = isempty(grouped) ? [shapes[index]] :
                    [shapes[i] for i in [only(grouped).contents; only(grouped).region]]
            isannular(bands) ||
                throw(ArgumentError(
                    "flatten reduces dielectric layers that are annuli about the center; " *
                    ":$(design.regions[index].source.tag) of the cable $(design.id) is not"))
        end
    end
    return design
end
