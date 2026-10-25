"""
$(TYPEDSIGNATURES)

Apply the multiplicative overlength of every longitudinal path to conductor
resistance \\[Ω/m\\]. An empty path tuple leaves the resistance unchanged.
"""
function path_corrected_resistance(resistance::T, ::Tuple{}) where {T}
    resistance
end

function path_corrected_resistance(resistance::T, paths::Tuple) where {T}
    inner = path_corrected_resistance(resistance, Base.front(paths))
    entry = last(paths)
    return inner * convert(T, overlength(entry.path, entry.radius))
end

"""
$(TYPEDSIGNATURES)

Add the turns per unit length contributed by a tuple of longitudinal paths
\\[1/m\\].
"""
function turns_per_length(paths::Tuple, ::Type{T}) where {T <: Real}
    return sum(
        entry -> inv(pitch(entry.path, entry.radius)),
        paths;
        init = zero(T)
    )
end

"""
$(TYPEDSIGNATURES)

Calculate the equivalent radial permeability of dielectric layers and apply
the helical-solenoid correction associated with the enclosed conductor.

# Arguments

- `layers`: Ordered dielectric layers with radii \\[m\\] and material
  permeability \\[dimensionless\\].
- `turns`: Equivalent conductor turns per unit length \\[1/m\\].
- `conductor_radius`: Outer conductor radius \\[m\\].
- `dielectric_radius`: Outer dielectric radius \\[m\\].

# Returns

- Corrected relative permeability \\[dimensionless\\].
"""
function equivalent_dielectric_permeability(
        layers,
        turns::Real,
        conductor_radius::Real,
        dielectric_radius::Real
)
    isempty(layers) && throw(ArgumentError(
        "equivalent dielectric permeability requires at least one layer"
    ))
    weighted = sum(layers) do layer
        layer.material.mu_r * log(layer.r_ex / layer.r_in)
    end
    radial = weighted / log(dielectric_radius / conductor_radius)
    return radial * solenoid_factor(
        turns,
        conductor_radius,
        dielectric_radius
    )
end

"""
$(TYPEDSIGNATURES)

Reduce one homogeneous dielectric region to concentric radii and material
properties. Primitive-specific methods retain circular and equivalent-area
cable-sector assumptions.
"""
function dielectric_layer(
        primitive::Annulus,
        source::PlacedRegion,
        ::Type{T}
) where {T <: Real}
    isempty(source.paths) || throw(ArgumentError(
        "flatten does not support helical dielectric layers"
    ))
    centre = (convert(T, primitive.at.x), convert(T, primitive.at.y))
    material = source.source.material
    material.kind === :conductor && throw(ArgumentError(
        "a dielectric interval cannot contain a conductor material"
    ))
    return (
        r_in = convert(T, r_in(source.primitive)),
        r_ex = convert(T, r_ex(source.primitive)),
        position = centre,
        material
    )
end

function dielectric_layer(
        primitive::ShellShape{
            <:Any,
            <:SectorShape,
            <:SectorShape
        },
        source::PlacedRegion,
        ::Type{T}
) where {T <: Real}
    isempty(source.paths) || throw(ArgumentError(
        "flatten does not support a helical sector dielectric"
    ))
    material = source.source.material
    material.kind === :conductor && throw(ArgumentError(
        "a dielectric interval cannot contain a conductor material"
    ))
    inner_area = convert(T, area(primitive.inner))
    outer_area = convert(T, area(primitive.outer))
    inner_radius = sqrt(inner_area / (one(T) * pi))
    outer_radius = sqrt(outer_area / (one(T) * pi))
    centre = convert.(T, centroid(primitive.inner))
    return (
        r_in = inner_radius,
        r_ex = outer_radius,
        position = centre,
        material
    )
end

function dielectric_layer(
        primitive::AbstractShape,
        source::PlacedRegion,
        ::Type
)
    throw(ArgumentError(
        "flatten does not support dielectric primitive " *
        "$(nameof(typeof(primitive)))"
    ))
end

function radial_extent(shape::Disk, centre)
    return hypot(shape.at.x - centre[1], shape.at.y - centre[2]) + shape.r
end

function radial_extent(shape::Annulus, centre)
    return hypot(shape.at.x - centre[1], shape.at.y - centre[2]) + shape.ro
end

function radial_extent(shape::Polygon, centre)
    cosine = cos(shape.at.φ)
    sine = sin(shape.at.φ)
    return maximum(shape.points) do point
        x = shape.at.x + cosine * point[1] - sine * point[2]
        y = shape.at.y + sine * point[1] + cosine * point[2]
        hypot(x - centre[1], y - centre[2])
    end
end

function radial_extent(shape::BentStrip, centre)
    dx = shape.at.x - centre[1]
    dy = shape.at.y - centre[2]
    offset = hypot(dx, dy)
    iszero(offset) && return shape.ro
    direction = atan(dy, dx)
    period = oftype(direction, 2pi)
    relative = mod(direction - shape.at.φ + pi, period) - pi
    nearest = clamp(relative, -shape.span / 2, shape.span / 2)
    cosine = cos(relative - nearest)
    return sqrt(max(
        offset^2 + shape.ri^2 + 2offset * shape.ri * cosine,
        offset^2 + shape.ro^2 + 2offset * shape.ro * cosine
    ))
end

function radial_extent(shape::Rectangle, centre)
    cosine = cos(shape.at.φ)
    sine = sin(shape.at.φ)
    half_width = shape.w / 2
    half_height = shape.h / 2
    return maximum(
        ((-half_width, -half_height), (half_width, -half_height),
            (half_width, half_height), (-half_width, half_height))
    ) do point
        x = shape.at.x + cosine * point[1] - sine * point[2]
        y = shape.at.y + sine * point[1] + cosine * point[2]
        hypot(x - centre[1], y - centre[2])
    end
end

# A sector conductor becomes the circle of its area about its centroid, as its
# sector dielectric layers do.
radial_extent(shape::SectorShape, centre) = sqrt(area(shape) / pi)

"""
$(TYPEDSIGNATURES)

Calculate the mutual GMD of two resolved strand sections \\[m\\].
General shapes retain the centroid-distance approximation. An annulus
surrounding the other section uses the uniform-current area average from
Rosa [rosa1908](@cite), Eq. (57):

```math
\\log D = \\frac{r_o^2\\log r_o-r_i^2\\log r_i}{r_o^2-r_i^2}-\\frac12.
```

# Arguments

- `left`, `right`: Resolved strand cross-sections, with coordinates in \\[m\\].

# Returns

- Mutual geometric mean distance \\[m\\].

# Notes

The annular expression is independent of the enclosed section's shape or
eccentricity while that section remains entirely inside the annular hole.
It does not replace the separate annular self-GMR or imply an exact integral
for arbitrary deformed strands.

# Errors

- `ArgumentError`: Coincident centroids without an enclosing annular section.
"""
function geometric_mean_distance(left::AbstractShape, right::AbstractShape)
    left_centre, right_centre = centroid(left), centroid(right)
    distance = hypot(left_centre[1] - right_centre[1], left_centre[2] - right_centre[2])
    distance > zero(distance) || throw(ArgumentError(
        "nested conductor members must have distinct centroids unless one is an enclosing annulus"
    ))
    return distance
end

function geometric_mean_distance(left::Annulus, right::AbstractShape)
    tolerance = 64eps(float(nominal(left.ro)))
    if left.ri > zero(left.ri) && radial_extent(right, centroid(left)) <= left.ri + tolerance
        # Dimensionless radius ratio avoids subtraction of nearly equal
        # squared radii in Rosa's area average for a thin annulus.
        ratio = (left.ro - left.ri) / left.ri
        return left.ri * exp(log1p(ratio) * (1 + ratio)^2 / (ratio * (2 + ratio)) - one(ratio) / 2)
    end
    return invoke(geometric_mean_distance, Tuple{AbstractShape, AbstractShape}, left, right)
end

geometric_mean_distance(left::AbstractShape, right::Annulus) =
    geometric_mean_distance(right, left)

function geometric_mean_distance(left::Annulus, right::Annulus)
    outer, inner = left.ro >= right.ro ? (left, right) : (right, left)
    return invoke(geometric_mean_distance, Tuple{Annulus, AbstractShape}, outer, inner)
end

"""
$(TYPEDSIGNATURES)

Reduce the conductive regions of one terminal to parallel resistance,
current-share-weighted GMR, equivalent temperature coefficient, turns, and
radial extent. Each region has the resistance of its resistivity over its
area, times the overlength of its paths.
"""
function nested_conductor_zone(
        sources,
        ::Type{T},
        expected_inner::T
) where {T <: Real}
    isempty(sources) && throw(ArgumentError(
        "a nested conductor requires at least one resolved primitive"
    ))
    all(sources) do source
        source.primitive isa Union{Disk, Rectangle, Polygon, BentStrip, Annulus, SectorShape}
    end || throw(ArgumentError(
        "flatten requires conductor regions resolved as Disk, Rectangle, Polygon, " *
        "BentStrip, Annulus, or SectorShape geometry"
    ))

    materials = Material{T}[convert(Material{T}, source.source.material)
                            for source in sources]
    reference = first(materials).T0
    all(material -> isapprox(material.T0, reference), materials) || throw(
        ArgumentError("all nested conductor materials must share one reference temperature")
    )

    # Compaction keeps the area a strand declares; a shell or a sector has the
    # area of its resolved shape.
    areas = T[convert(T, area(source.source.primitive isa Union{Disk, Rectangle} ?
                             source.source.primitive : source.primitive))
              for source in sources]
    all(eachindex(sources)) do index
        isapprox(
            convert(T, area(sources[index].primitive)),
            areas[index];
            rtol = 5.0e-6,
            atol = zero(T)
        )
    end || throw(ArgumentError(
        "resolved compacted strands must preserve every declared source area"
    ))
    equivalent_radii = sqrt.(areas ./ (one(T) * π))
    total_area = sum(areas)
    centre = convert.(T, conductor_zone_position(sources))
    resistances = T[path_corrected_resistance(
                        materials[index].rho / areas[index],
                        sources[index].paths
                    ) for index in eachindex(sources)]

    resistance = first(resistances)
    alpha = first(materials).alpha
    for index in Iterators.drop(eachindex(sources), 1)
        alpha = equivalent_alpha(
            alpha,
            resistance,
            materials[index].alpha,
            resistances[index]
        )
        resistance = parallel(resistance, resistances[index])
    end

    conductances = inv.(resistances)
    weights = conductances ./ sum(conductances)
    log_gmr = zero(T)
    for left in eachindex(sources)
        primitive = sources[left].primitive
        self_gmr = primitive isa Annulus ?
                   tubular_gmr(primitive.ro, primitive.ri, materials[left].mu_r) :
                   equivalent_radii[left] * exp(-materials[left].mu_r / 4)
        log_gmr += weights[left]^2 * log(self_gmr)
        for right in (left + 1):length(sources)
            distance = geometric_mean_distance(sources[left].primitive, sources[right].primitive)
            log_gmr += 2weights[left] * weights[right] * log(distance)
        end
    end

    turn_values = T[turns_per_length(source.paths, T) for source in sources]
    turns = sum(weights .* turn_values)
    bounded = map(sources) do source
        index = findfirst(
            entry -> entry.pattern isa BoundedPlacement,
            source.placement.patterns
        )
        index === nothing && return nothing
        any(entry -> entry.owner === Group, source.placement.patterns[(index + 1):end]) &&
            return nothing
        return source.placement.patterns[index].pattern.boundary
    end
    outer = if all(!isnothing, bounded)
        envelope = first(bounded)
        all(==(envelope), bounded) || throw(ArgumentError(
            "one nested conductor cannot combine different bounded formations"
        ))
        iszero(expected_inner) || throw(ArgumentError(
            "a bounded stranded formation must be the innermost conductor"
        ))
        sqrt(convert(T, area(envelope)) / (one(T) * π))
    else
        maximum(
            index -> radial_extent(sources[index].primitive, centre),
            eachindex(sources)
        )
    end
    outer > expected_inner || throw(ArgumentError(
        "nested conductor envelope must exceed its preceding radial boundary"
    ))
    return (
        r_in = expected_inner,
        r_ex = outer,
        cross_section = total_area,
        num_wires = length(sources),
        num_turns = turns,
        resistance,
        alpha,
        gmr = exp(log_gmr),
        reference_temperature = reference,
        position = centre
    )
end

radial_position(shape::Union{Disk, Annulus}) = (shape.at.x, shape.at.y)
radial_position(shape::SectorShape) = centroid(shape)
function radial_position(
        shape::ShellShape{
        <:Any,
        <:SectorShape,
        <:SectorShape
}
)
    centroid(shape.inner)
end
radial_position(region::PlacedRegion) = radial_position(region.primitive)

function conductor_zone_position(sources)
    isempty(sources) && throw(ArgumentError(
        "a conductor zone requires at least one resolved region"
    ))
    bounded = map(sources) do source
        index = findfirst(
            entry -> entry.pattern isa BoundedPlacement,
            source.placement.patterns
        )
        index === nothing && return nothing
        shape = source.placement.patterns[index].pattern.boundary
        return shape isa Disk ? (shape.at.x, shape.at.y) : centroid(shape)
    end
    if all(!isnothing, bounded)
        reference = first(bounded)
        all(position -> same_radial_position(position, reference), bounded) &&
            return reference
    end
    patterned = map(sources) do source
        entries = filter(entry -> entry.owner === Group, source.placement.patterns)
        length(entries) == 1 || return nothing
        primitive = source.primitive
        hasproperty(primitive, :at) || return nothing
        intrinsic = source.source.primitive
        hasproperty(intrinsic, :at) || return nothing
        local_pose = only(entries).pose * intrinsic.at
        absolute = primitive.at
        parent_angle = absolute.φ - local_pose.φ
        return (
            absolute.x - cos(parent_angle) * local_pose.x +
            sin(parent_angle) * local_pose.y,
            absolute.y - sin(parent_angle) * local_pose.x -
            cos(parent_angle) * local_pose.y
        )
    end
    if all(!isnothing, patterned)
        reference = first(patterned)
        all(position -> same_radial_position(position, reference), patterned) &&
            return reference
    end
    length(sources) == 1 && return centroid(only(sources).primitive)
    areas = map(source -> area(source.primitive), sources)
    centres = map(source -> centroid(source.primitive), sources)
    total = sum(areas)
    return (
        sum(index -> areas[index] * centres[index][1], eachindex(sources)) / total,
        sum(index -> areas[index] * centres[index][2], eachindex(sources)) / total
    )
end

function same_radial_position(left, right)
    coordinates = promote(
        float(left[1]),
        float(left[2]),
        float(right[1]),
        float(right[2])
    )
    T = typeof(first(coordinates))
    scale = max(one(T), maximum(abs, coordinates))
    tolerance = sqrt(eps(T)) * scale
    return hypot(
        coordinates[1] - coordinates[3],
        coordinates[2] - coordinates[4]
    ) <= tolerance
end

function same_radial_chain(left::PlacedRegion, right_position)
    return same_radial_position(radial_position(left), right_position)
end

"""
$(TYPEDSIGNATURES)

Convert a completed scalar conductor reduction into the homogeneous material
that reproduces its DC resistance and GMR.
"""
function equivalent_conductor_material(conductor)
    return Material(
        :conductor,
        equivalent_rho(
            conductor.resistance,
            conductor.r_ex,
            conductor.r_in
        ),
        zero(conductor.r_in),
        equivalent_mu(
            conductor.gmr,
            conductor.r_ex,
            conductor.r_in
        ),
        conductor.reference_temperature,
        conductor.alpha
    )
end

"""
$(TYPEDSIGNATURES)

Calculate the shunt capacitance \\[F/m\\] and conductance \\[S/m\\] of one radial
dielectric layer at `angular_frequency` \\[rad/s\\]. Polarization loss tangent
excludes the conduction already represented by material resistivity.
"""
function dielectric_circuit(layer, angular_frequency::Real)
    capacitance = shunt_capacitance(layer.r_in, layer.r_ex, layer.material.eps_r)
    return (
        capacitance = capacitance,
        conductance = shunt_conductance(
            layer.r_in,
            layer.r_ex,
            layer.material.rho
        ) + angular_frequency * capacitance * layer.material.tan_delta
    )
end

function dielectric_circuit(layers, angular_frequency::Real, ::Type{T}) where {T <: Real}
    isempty(layers) && return (
        capacitance = zero(T),
        conductance = zero(T)
    )
    combined = dielectric_circuit(first(layers), angular_frequency)
    @inbounds for layer in Iterators.drop(layers, 1)
        circuit = dielectric_circuit(layer, angular_frequency)
        combined = series_shunt_admittance(
            combined.conductance,
            combined.capacitance,
            circuit.conductance,
            circuit.capacitance,
            angular_frequency
        )
    end
    return combined
end

"""
$(TYPEDSIGNATURES)

Initialize the frequency-independent radial dielectric description from its
innermost layer.
"""
function initialize_dielectric(layer, conductor)
    same_radial_position(layer.position, conductor.position) || throw(ArgumentError(
        "dielectric layers must share the conductor radial centre"
    ))
    isapprox(layer.r_in, conductor.r_ex) || throw(ArgumentError(
        "the first dielectric layer must begin at the conductor boundary"
    ))
    isapprox(layer.material.T0, conductor.reference_temperature) || throw(
        ArgumentError("all cable materials must share one reference temperature")
    )
    return (
        r_in = layer.r_in,
        r_ex = layer.r_ex,
        cross_section = (one(layer.r_in) * pi) *
                        (layer.r_ex^2 - layer.r_in^2),
        position = layer.position,
        reference_temperature = layer.material.T0,
        layers = dielectric_layer(layer.material, layer.r_in, layer.r_ex)
    )
end

"""
$(TYPEDSIGNATURES)

Add one radial dielectric layer to a frequency-independent radial description.
"""
function add_dielectric_layer(dielectric, layer)
    same_radial_position(layer.position, dielectric.position) || throw(ArgumentError(
        "dielectric layers must share the conductor radial centre"
    ))
    isapprox(layer.r_in, dielectric.r_ex) || throw(ArgumentError(
        "dielectric layers must be radially contiguous"
    ))
    isapprox(layer.material.T0, dielectric.reference_temperature) || throw(
        ArgumentError("all cable materials must share one reference temperature")
    )
    layers = copy(dielectric.layers)
    append!(layers, dielectric_layer(layer.material, layer.r_in, layer.r_ex))
    return (
        r_in = dielectric.r_in,
        r_ex = layer.r_ex,
        cross_section = dielectric.cross_section +
                        (one(layer.r_in) * pi) * (layer.r_ex^2 - layer.r_in^2),
        position = dielectric.position,
        reference_temperature = dielectric.reference_temperature,
        layers
    )
end

function dielectric_layer(material::Material, inner::T, outer::T) where {T <: Real}
    return [(r_in = inner, r_ex = outer, material = convert(Material{T}, material))]
end

function dielectric_layer(material::RadialDielectric, inner::T, outer::T) where {T <: Real}
    # Reconstruct virtual radial intervals only for computation. The homogeneous
    # design still contains one physical shell, not these constituent surfaces.
    total = convert(T, sum(material.weights))
    width = log(outer / inner)
    radial_mu = sum(w * m.mu_r for (w, m) in zip(material.weights, material.materials)) / total
    mu_scale = material.mu_r / radial_mu
    accumulated = zero(T)
    previous = inner
    layers = NamedTuple{(:r_in, :r_ex, :material), Tuple{T, T, Material{T}}}[]
    for index in eachindex(material.materials)
        accumulated += convert(T, material.weights[index])
        next = index == lastindex(material.materials) ? outer :
               inner * exp(width * accumulated / total)
        source = material.materials[index]
        physical = Material(source.kind, source.rho, source.eps_r, source.mu_r * mu_scale,
            source.T0, source.alpha; rho_thermal=source.rho_thermal,
            theta_max=source.theta_max, tan_delta=source.tan_delta, sigma_solar=source.sigma_solar)
        push!(layers, (r_in = previous, r_ex = next,
            material = convert(Material{T}, physical)))
        previous = next
    end
    return layers
end

function empty_dielectric(conductor, ::Type{T}) where {T <: Real}
    Layer = NamedTuple{
        (:r_in, :r_ex, :material),
        Tuple{T, T, Material{T}}
    }
    return (
        r_in = conductor.r_ex,
        r_ex = conductor.r_ex,
        cross_section = zero(T),
        position = conductor.position,
        reference_temperature = conductor.reference_temperature,
        layers = Layer[]
    )
end

"""
$(TYPEDSIGNATURES)

Convert a completed scalar dielectric reduction into a homogeneous material
that matches its capacitance and conductance at `angular_frequency` \\[rad/s\\].
Resistivity retains the series DC conduction. The remaining conductance is
represented by polarization loss tangent, including interfacial polarization
from heterogeneous layers. The equivalent is a reference-frequency export,
not a broadband replacement for the original physical layers.

Relative permeability retains the radial material average and the conductor's
helical-solenoid correction.
"""
function equivalent_dielectric_material(dielectric, conductor, terminal::Symbol,
        angular_frequency::Real)
    isempty(dielectric.layers) && return Material(
        :insulator,
        oftype(dielectric.r_in, Inf),
        zero(dielectric.r_in),
        one(dielectric.r_in),
        conductor.reference_temperature,
        zero(dielectric.r_in)
    )
    relative_mu = equivalent_dielectric_permeability(
        dielectric.layers,
        conductor.num_turns,
        conductor.r_ex,
        dielectric.r_ex
    )
    isfinite(relative_mu) || throw(ArgumentError(
        "dielectric reduction produced a non-finite permeability " *
        "for terminal :$terminal (turns=$(conductor.num_turns), " *
        "r_con=$(conductor.r_ex), r_ins=$(dielectric.r_ex))"
    ))
    dc_conductance = inv(sum(dielectric.layers) do layer
        inv(shunt_conductance(layer.r_in, layer.r_ex, layer.material.rho))
    end)
    polarization_conductance = dielectric.shunt_conductance - dc_conductance
    tolerance = 100 * eps(typeof(float(nominal(dc_conductance)))) *
                max(abs(dielectric.shunt_conductance), abs(dc_conductance))
    polarization_conductance >= -tolerance || throw(DomainError(
        polarization_conductance,
        "equivalent dielectric conductance must not be smaller than its DC limit"))
    polarization_conductance = max(zero(polarization_conductance), polarization_conductance)
    loss_tangent = iszero(polarization_conductance) ? zero(polarization_conductance) :
                   polarization_conductance / (angular_frequency * dielectric.shunt_capacitance)
    return Material(
        :insulator,
        inv(equivalent_conductivity(
            dc_conductance,
            dielectric.r_in,
            dielectric.r_ex
        )),
        equivalent_eps(
            dielectric.shunt_capacitance,
            dielectric.r_ex,
            dielectric.r_in
        ),
        relative_mu,
        conductor.reference_temperature,
        zero(dielectric.r_in);
        tan_delta = loss_tangent
    )
end

function _radial_regions(regions)
    radial = PlacedRegion[]
    sizehint!(radial, length(regions))
    for source in regions
        source.primitive isa DifferenceShape || begin
            push!(radial, source)
            continue
        end
        material = source.source.material
        source.terminal === nothing || throw(ArgumentError(
            "a non-radial enclosure fill cannot own a terminal"
        ))
        material.kind === :insulator || throw(ArgumentError(
            "a non-radial enclosure fill must use an insulator material"
        ))
        isapprox(material.mu_r, one(material.mu_r)) || throw(ArgumentError(
            "the coaxial reduction requires a nonmagnetic non-radial enclosure fill"
        ))
    end
    return radial
end

"""
$(TYPEDSIGNATURES)

Reduce each radial terminal section of a completed cable design to equivalent
conductor geometry and its ordered physical dielectric layers.

The reduction combines conductor resistance and GMR but does not evaluate a
dielectric circuit. It performs no constitutive, frequency, line-parameter,
mutual-coupling, earth-return, or matrix calculation.

# Arguments

- `design`: Completed physical cable design.
# Returns

An ordered vector of named tuples containing equivalent conductor fields and
the uncombined physical dielectric layers for each retained terminal.
"""
function radial_components(design::CableDesign, ::Type{T}) where {T <: Real}
    regions = _radial_regions(design.geometry.regions)
    terminals = design.terminal_order
    isempty(terminals) && throw(ArgumentError(
        "flatten requires at least one retained terminal"
    ))
    terminal_starts = Int[]
    for terminal in terminals
        index = findfirst(region -> region.terminal === terminal, regions)
        index === nothing && throw(ArgumentError(
            "terminal :$terminal has no resolved conductive region"
        ))
        push!(terminal_starts, index)
    end
    issorted(terminal_starts) && allunique(terminal_starts) || throw(ArgumentError(
        "flatten requires nonreappearing radial terminal blocks"
    ))

    Layer = NamedTuple{
        (:r_in, :r_ex, :material),
        Tuple{T, T, Material{T}}
    }
    ConductorInput = NamedTuple{
        (:r_in, :r_ex, :cross_section, :num_wires, :num_turns,
            :resistance, :alpha, :gmr, :position, :material),
        Tuple{T, T, T, Int, T, T, T, T, Tuple{T, T}, Material{T}}
    }
    DielectricInput = NamedTuple{
        (:r_in, :r_ex, :cross_section, :layers),
        Tuple{T, T, T, Vector{Layer}}
    }
    ComponentInput = NamedTuple{
        (:name, :conductor, :dielectric),
        Tuple{Symbol, ConductorInput, DielectricInput}
    }
    effective = ComponentInput[]

    for (terminal_index, terminal) in enumerate(terminals)
        first_index = terminal_starts[terminal_index]
        last_index = terminal_index == length(terminals) ? length(regions) :
                     terminal_starts[terminal_index + 1] - 1
        block = @view regions[first_index:last_index]

        cursor = something(
            findfirst(region -> region.terminal !== terminal, block),
            lastindex(block) + 1
        )
        any(index -> block[index].terminal !== nothing, cursor:lastindex(block)) &&
            throw(ArgumentError(
                "terminal :$terminal reappears after its dielectric interval"
            ))

        conductor_sources = @view block[firstindex(block):(cursor - 1)]
        zone_position = conductor_zone_position(conductor_sources)
        expected_inner = if first_index > firstindex(regions) &&
                            same_radial_chain(regions[first_index - 1], zone_position)
            convert(T, r_ex(regions[first_index - 1].primitive))
        elseif first(conductor_sources).primitive isa Annulus
            # Nothing fills the hole of an innermost tube.
            convert(T, first(conductor_sources).primitive.ri)
        else
            zero(T)
        end
        conductor = nested_conductor_zone(
            conductor_sources,
            T,
            expected_inner
        )

        dielectric_sources = cursor > lastindex(block) ?
                             @view(block[1:0]) : @view(block[cursor:lastindex(block)])
        if first_index > firstindex(regions) &&
           same_radial_chain(regions[first_index - 1], conductor.position)
            preceding = regions[first_index - 1]
            preceding.terminal === nothing || throw(ArgumentError(
                "a radial terminal must be preceded by a dielectric region"
            ))
            preceding_layer = dielectric_layer(
                preceding.primitive,
                preceding,
                T
            )
            conductor = merge(conductor, (r_in = preceding_layer.r_ex,))
        end
        if !isempty(dielectric_sources)
            first_dielectric = first(dielectric_sources)
            first_layer = dielectric_layer(
                first_dielectric.primitive,
                first_dielectric,
                T
            )
            conductor = merge(conductor, (r_ex = first_layer.r_in,))
        end

        dielectric = empty_dielectric(conductor, T)
        for (layer_index, source) in enumerate(dielectric_sources)
            source.terminal === nothing || throw(ArgumentError(
                "a dielectric region cannot own a retained terminal"
            ))
            layer = dielectric_layer(
                source.primitive,
                source,
                T
            )
            dielectric = layer_index == 1 ?
                         initialize_dielectric(layer, conductor) :
                         add_dielectric_layer(dielectric, layer)
        end

        conductor_material = equivalent_conductor_material(conductor)
        push!(effective,
            (
                name = terminal,
                conductor = (
                    r_in = conductor.r_in,
                    r_ex = conductor.r_ex,
                    cross_section = conductor.cross_section,
                    num_wires = conductor.num_wires,
                    num_turns = conductor.num_turns,
                    resistance = conductor.resistance,
                    alpha = conductor.alpha,
                    gmr = conductor.gmr,
                    position = conductor.position,
                    material = conductor_material
                ),
                dielectric = (
                    r_in = dielectric.r_in,
                    r_ex = dielectric.r_ex,
                    cross_section = dielectric.cross_section,
                    layers = dielectric.layers
                )
            ))
    end
    return effective
end

function radial_components(design::CableDesign)
    T = promote_type(
        (eltype(region.primitive) for region in design.geometry.regions)...,
        (eltype(region.source.material) for region in design.geometry.regions)...
    )
    return radial_components(design, T)
end

"""
$(TYPEDSIGNATURES)

Reduce a completed cable design to equivalent radial conductor and dielectric
components at one dielectric reference frequency.

This reference-frequency operation serves the scalar ATP and TRALIN export
adapters. It includes the supplied physical material losses at their reference
temperature and combines radial layers in series. The resulting scalar
material reproduces capacitance and conductance at the requested frequency;
it is not a broadband replacement. It does not calculate mutual coupling or
earth return. Reusable [`homogenize`](@ref) designs instead retain their
constituents without calling this frequency-dependent operation.

# Arguments

- `design`: Completed physical cable design.
- `dielectric_frequency`: Reference frequency used to match the equivalent
  dielectric [Hz].

# Returns

- An ordered vector of equivalent radial component records.
"""
function flatten(
        design::CableDesign,
        dielectric_frequency::Real
)
    T = promote_type(
        typeof(float(dielectric_frequency)),
        (eltype(region.primitive) for region in design.geometry.regions)...,
        (eltype(region.source.material) for region in design.geometry.regions)...
    )
    return flatten(
        design,
        convert(T, float(dielectric_frequency)),
        T
    )
end

function flatten(
        design::CableDesign,
        frequency::T,
        ::Type{T}
) where {T <: Real}
    frequency > zero(frequency) || throw(DomainError(
        frequency,
        "dielectric reference frequency must be positive"
    ))
    components = radial_components(design, T)
    angular_frequency = 2 * (one(T) * pi) * frequency
    Layer = NamedTuple{
        (:r_in, :r_ex, :material),
        Tuple{T, T, Material{T}}
    }
    ConductorInput = typeof(first(components).conductor)
    DielectricInput = NamedTuple{
        (:r_in, :r_ex, :cross_section, :shunt_capacitance,
            :shunt_conductance, :frequency, :layers, :material),
        Tuple{T, T, T, T, T, T, Vector{Layer}, Material{T}}
    }
    ComponentInput = NamedTuple{
        (:name, :conductor, :dielectric),
        Tuple{Symbol, ConductorInput, DielectricInput}
    }
    effective = ComponentInput[]
    sizehint!(effective, length(components))
    for component in components
        circuit = dielectric_circuit(
            component.dielectric.layers,
            angular_frequency,
            T
        )
        dielectric = merge(component.dielectric,
            (
                shunt_capacitance = circuit.capacitance,
                shunt_conductance = circuit.conductance
            ))
        material = equivalent_dielectric_material(
            dielectric,
            merge(component.conductor, (
                reference_temperature = component.conductor.material.T0,
            )),
            component.name,
            angular_frequency
        )
        push!(effective,
            (
                name = component.name,
                conductor = component.conductor,
                dielectric = (
                    r_in = dielectric.r_in,
                    r_ex = dielectric.r_ex,
                    cross_section = dielectric.cross_section,
                    shunt_capacitance = circuit.capacitance,
                    shunt_conductance = circuit.conductance,
                    frequency,
                    layers = dielectric.layers,
                    material
                )
            ))
    end
    return effective
end

"""
$(TYPEDSIGNATURES)

Build a homogeneous cable design from locally reduced radial components.

Each retained terminal becomes one solid or annular conductor followed by one
homogeneous dielectric interval. The conductor material reproduces the
terminal's parallel resistance and geometric-mean radius. The dielectric
material retains the physical constituents and radial weights in a
[`RadialDielectric`](@ref). Constitutive laws are selected only during computation
or explicit reference-frequency export. The source design is not modified.

# Arguments

- `original`: Completed physical cable design.

# Keywords

- `new_id`: Identifier for the returned design. An empty value appends
  `"_equivalent"` to the source identifier.

# Returns

- A completed homogeneous [`CableDesign`](@ref).
"""
function homogenize(
        original::CableDesign;
        new_id::AbstractString = ""
)
    target_id = isempty(new_id) ? original.cable_id * "_equivalent" : String(new_id)
    components = radial_components(original)
    chains = Vector{typeof(components)}()
    for component in components
        if isempty(chains) ||
           !same_radial_position(
            first(last(chains)).conductor.position,
            component.conductor.position
        )
            any(
                chain -> same_radial_position(
                    first(chain).conductor.position,
                    component.conductor.position
                ),
                chains) && throw(ArgumentError(
                "a radial assembly member cannot reappear after another member"
            ))
            push!(chains, eltype(components)[])
        end
        push!(last(chains), component)
    end

    flattened_members = AbstractCablePart[]
    member_positions = Tuple[]
    for chain in chains
        parts = AbstractCablePart[]
        for component in chain
            terminal = component.name
            conductor = component.conductor
            conductor_primitive = iszero(conductor.r_in) ?
                                  Disk(conductor.r_ex) :
                                  Annulus(conductor.r_in, conductor.r_ex)
            conductor_region = Region(
                Symbol(terminal, :_equivalent_conductor),
                conductor_primitive,
                conductor.material
            )
            push!(parts,
                Group(
                    terminal,
                    Pose2(0, 0, 0),
                    conductor_region,
                    nothing,
                    nothing,
                    nothing
                ))

            dielectric = component.dielectric
            if dielectric.r_ex > dielectric.r_in
                push!(parts,
                    Region(
                        Symbol(terminal, :_equivalent_dielectric),
                        Annulus(dielectric.r_in, dielectric.r_ex),
                        RadialDielectric(
                            [layer.material for layer in dielectric.layers],
                            [log(layer.r_ex / layer.r_in) for layer in dielectric.layers];
                            mu_r = equivalent_dielectric_permeability(dielectric.layers,
                                conductor.num_turns, conductor.r_ex, dielectric.r_ex))
                    ))
            end
        end
        push!(flattened_members, Stack(parts))
        push!(member_positions, first(chain).conductor.position)
    end

    origin = (zero(first(member_positions)[1]), zero(first(member_positions)[2]))
    singleton = length(flattened_members) == 1 &&
                same_radial_position(only(member_positions), origin)
    root = if singleton
        only(flattened_members)
    else
        members = Tuple(
            AssemblyMember(
                member,
                Pose2(position[1], position[2], 0)
            )
        for (member, position) in zip(flattened_members, member_positions)
        )
        Assembly(
            Pose2(0, 0, 0),
            members,
            nothing,
            nothing,
            nothing,
            nothing
        )
    end
    return build(CableDesign, target_id, root; nominal_data = original.nominal_data)
end
