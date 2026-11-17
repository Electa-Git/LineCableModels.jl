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
    "Material of one region, or the composite material of several regions."
    material::M
    function BlueprintDielectric(r_in::T, r_ex::T,
            material::M) where {T <: Real, M <: Material}
        return new{T, M}(r_in, r_ex, material)
    end
end

"""
$(TYPEDSIGNATURES)

Build the layer from `r_in` to `r_ex`, both in \\[m\\], of the dielectric regions of
materials `media` and cross-sections `bands` as one composite material. Each region weighs
by its area fraction ``f_k``:

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
    The relative permeability of the interval of each conductor, with the correction for the
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
    Rows of the members of each assembly, in the order of their conductors. A sub-cable with
    no conductor of its own gives its members to its container.
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
regions of its first terminal, or about its placement when it has no terminal.

For an annular block, the tube starts at the nearest extent of the terminal's conductive
regions and the containers of the fills between them. It ends at the farthest extent of
its conductive regions. The fill beyond the metal forms the interval's first layer, with
the fill's material. This layer has zero width where the metal meets the fill's outer boundary.
A block that is not annular takes Ametani's equivalent below, with the fill's outline as
its outline. No fill layer forms.

The regions of the same row outside the tube, up to its next terminal or its last region,
become the layers of its dielectric interval. The regions of a row around its members form
the interior of the assembly. The annuli about the axis outside its other regions are its
outer layers, from the inside out. The other regions form its filler, one layer from their
nearest to their farthest extent. Regions inside the first terminal of a row without members
lie in no interval.

A sub-cable with no conductor of its own gives its members to its container, and its regions
join the container's interior.

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

A disk or annulus ``S_i`` about the center, of outer radius ``b``, has the own term

```math
\\log g_i = \\log b - \\mu_i (\\log b - \\log \\operatorname{gmd}(S_i,S_i)).
```

The difference in parentheses is the tube's internal term, equal to ``1/4``
for a disk. A thin shell at equal inner and outer radii has the own term ``\\log b``.
Any other region ``S_i`` has the own term
``\\log g_i = \\log \\operatorname{gmd}(S_i,S_i) - (\\mu_i-1)/4``.
Its geometric mean distance comes from its outline through [`gmd`](@ref), with the
circle's convention for the permeability correction. Every pair of regions contributes the cross term
``\\log d_{ij} = \\log \\operatorname{gmd}(S_i,S_j)``.

The calculation evaluates the own term once for each position within a repeated part.
That position in each copy takes the value from the first copy. The copies are rigid images under the repeat's placements.
A region outside every repeat takes its own term.

For a block that is not annular, Ametani's equivalent gives the tube's radii from its metal
area ``S`` and outline perimeter ``l``:

```math
r_e = \\frac{l}{2\\pi}, \\qquad r_i = \\sqrt{r_e^2 - \\frac{S}{\\pi}}.
```

An inner radius within `tolerance(r_e)` of zero, decided on `nominal`, is zero.
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
    \\left[1 + \\frac{2 \\pi^2 N^2 (r_d^2 - r_e^2)}{\\log(r_d / r_e)}\\right].
```

The calculation uses ``(r_d^2 - r_e^2) / \\log(r_d / r_e) =
2 r_d r_e \\sinh(\\log(r_d / r_e)) / \\log(r_d / r_e)``.
The last factor takes its exact limit of 1 at ``r_d = r_e``.
For the first factor of ``\\mu_{ins}``, the calculation uses
``\\mu_1 + \\sum_{k \\ge 2} (\\mu_k - \\mu_1) \\log(b_k / a_k) / \\log(r_d / r_e)``,
equal to the quotient because layers tile ``r_e`` to ``r_d``, and ``\\mu_1`` at ``r_d = r_e``.
Without layers, ``\\mu_{ins} = 1``.

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
- Throws `ArgumentError` when a layer of an interval would leave part of its annulus about
  its assembly's axis uncovered and take that part as its material.
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
    bands = [first([[[fill.contents; index]
                     for fill in design.fills
                     if fill.region == index && isdisjoint(fill.contents, assigned)];
                    [[index]]])
             for index in eachindex(regions)
             if !any(fill -> index in fill.contents && isdisjoint(fill.contents, assigned),
        design.fills)]
    validate(design, flatten, owners, blocks, centered, bands)
    stops = [maximum([index
                      for index in design.assemblies[owners[k]].regions
                      if all(
                 j -> owners[j] != owners[k] || j <= k || index < first(blocks[j]),
                 eachindex(blocks))]) for k in eachindex(blocks)]
    representatives = collect(eachindex(regions))
    for row in design.repeats, copy in row.copies

        for (index, firstindex) in zip(copy, first(row.copies))
            representatives[index] = min(representatives[index], firstindex)
        end
    end
    for index in eachindex(representatives)
        representatives[index] = representatives[representatives[index]]
    end
    terms = zeros(T, length(regions))
    for index in unique(representatives[assigned])
        terminal = findfirst(row -> index in row.regions, design.terminals)
        shape = centered[owners[terminal]][index]
        μ = regions[index].source.material.mu_r
        terms[index] = if isannular(shape)
            a, b = distance(shape), support(shape)
            if isapprox(a, b)
                log(b)
            else
                log(b) - μ * (log(b) - log(gmd(shape, shape)))
            end
        else
            log(gmd(shape, shape)) - (μ - 1) / 4
        end
    end
    reduced = map(eachindex(design.terminals)) do k
        members = design.terminals[k].regions
        shapes = [centered[owners[k]][index] for index in members]
        materials = [regions[index].source.material for index in members]
        conductances = [area(shapes[i]) / (materials[i].rho * overlengths[members[i]])
                        for i in eachindex(members)]
        shares = conductances ./ sum(conductances)
        own = [terms[representatives[index]] for index in members]
        cross = sum(
            (shares[i] * shares[j] * log(gmd(shapes[i], shapes[j]))
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
            annular ? maximum(support, shapes) :
            perimeter(boundary(centered[owners[k]][last(blocks[k])])) / 2π)
        inner = convert(T,
            annular ? minimum(distance, tube) :
            (nominal(outer^2 - sum(area, shapes) / π) <= tolerance(outer)^2 ? zero(outer) :
             sqrt(outer^2 - sum(area, shapes) / π)))
        term = iszero(nominal(inner)) ? zero(outer) :
               inner^4 / (outer^2 - inner^2)^2 * log(outer / inner)
        mu = (log(outer) - gmr) / (term - (3inner^2 - outer^2) / (4(outer^2 - inner^2)))
        material = Material(Val(:conductor),
            convert(T, π * (outer^2 - inner^2) / sum(conductances)), one(T), convert(T, mu),
            convert(T, first(materials).T0),
            convert(T, sum(shares .* [item.alpha for item in materials])), zero(T),
            convert(T, 90), zero(T), zero(T))
        indices = [index
                   for index in design.assemblies[owners[k]].regions
                   if (annular && index in containers &&
                       all(i -> i <= index || i in contents, blocks[k])) ||
                      last(blocks[k]) < index <= stops[k]]
        interval = BlueprintDielectric[]
        for band in bands
            last(band) in indices || continue
            pieces = [centered[owners[k]][i] for i in band]
            media = [regions[i].source.material for i in band]
            lower = convert(T,
                annular && last(band) ∉ containers ? minimum(distance, pieces) :
                (isempty(interval) ? outer : last(interval).r_ex))
            upper = convert(
                T, annular ? maximum(support, pieces) :
                   lower + dimension(centered[owners[k]][last(band)]))
            push!(interval, BlueprintDielectric(lower, upper, media, pieces))
        end
        N = sum(shares .* [turns[index] for index in members])
        relative = if isempty(interval)
            one(T)
        else
            width = log(last(interval).r_ex / outer)
            (first(interval).material.mu_r +
             sum(
                layer -> (layer.material.mu_r - first(interval).material.mu_r) *
                         log(layer.r_ex / layer.r_in) / width,
                interval[2:end]; init = zero(T))) *
            (1 +
             4π^2 * N^2 * last(interval).r_ex * outer *
             (iszero(nominal(width)) ? one(width) : sinh(width) / width))
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
        filler = [band for band in bands if last(band) in indices]
        annuli = Vector{Int}[]
        while length(filler) > 1
            outside = shapes[last(filler)]
            reach = nominal(maximum(support, shapes[reduce(vcat, filler[1:(end - 1)])]))
            nearest = nominal(minimum(distance, outside))
            isannular(outside) && (nearest >= reach || isapprox(nearest, reach)) || break
            pushfirst!(annuli, pop!(filler))
        end
        isempty(filler) && return BlueprintDielectric[]
        return BlueprintDielectric[BlueprintDielectric(
                                       convert(T, minimum(distance, shapes[band])),
                                       convert(T, maximum(support, shapes[band])),
                                       [placed.source.material for placed in regions[band]],
                                       shapes[band])
                                   for band in [[sort!(reduce(vcat, filler))]; annuli]]
    end
    ends = [maximum(row.regions; init = 0) for row in design.assemblies]
    for k in reverse(eachindex(design.assemblies))
        ends[k] = maximum([ends[k]; ends[design.assemblies[k].members]])
    end
    order = sortperm(eachindex(design.assemblies); by = k -> (ends[k], -k))
    # Each row lists its members in the order of their conductors.
    foreach(rows -> sort!(rows; by = m -> (ends[m], -m)), contained)
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

function validate(design::CableDesign, ::typeof(flatten), owners::AbstractVector{Int},
        blocks::AbstractVector{<:AbstractVector{Int}},
        centered::AbstractVector{<:AbstractVector{<:AbstractShape}},
        bands::AbstractVector{<:AbstractVector{Int}})
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
            isannular([shapes[i] for i in only([band for band in bands if index in band])]) ||
                throw(ArgumentError(
                    "flatten assigns a dielectric layer an annulus about its assembly's axis, " *
                    "between the layer's nearest and farthest radii. " *
                    ":$(design.regions[index].source.tag) of the cable $(design.id) leaves " *
                    "part of that annulus uncovered. " *
                    "Its layer would take that part as its material."))
        end
    end
    return design
end
