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

Calculate the geometric mean distance between two conductor-zone element sets
\\[m\\]. Coincident element centres use the greater element radius as their
finite self-distance. Element areas weight heterogeneous zones by their
cross-sectional participation.

# Arguments

- `left_coordinates`: Element centres of the accumulated zone \\[m\\].
- `left_radius`: Element radius of the accumulated zone \\[m\\].
- `right_coordinates`: Element centres of the added zone \\[m\\].
- `right_radius`: Element radius of the added zone \\[m\\].
- `left_element_area`: Area represented by each accumulated element \\[m²\\].
- `right_element_area`: Area represented by each added element \\[m²\\].

# Returns

- Area-weighted geometric mean distance \\[m\\].
"""
function geometric_mean_distance(
        left_coordinates,
        left_radius::Real,
        right_coordinates,
        right_radius::Real,
        left_element_area::Real = one(float(left_radius)),
        right_element_area::Real = one(float(right_radius))
)
    isempty(left_coordinates) && throw(ArgumentError(
        "the first conductor zone requires at least one element coordinate"
    ))
    isempty(right_coordinates) && throw(ArgumentError(
        "the second conductor zone requires at least one element coordinate"
    ))
    logarithmic_sum, weights = promote(
        zero(float(left_radius)),
        zero(float(right_radius))
    )
    weight = left_element_area * right_element_area
    weight > zero(weight) || throw(DomainError(
        weight, "conductor-zone element areas must be positive"
    ))
    for left in left_coordinates, right in right_coordinates

        distance = hypot(left[1] - right[1], left[2] - right[2])
        effective_distance = iszero(distance) ? max(left_radius, right_radius) : distance
        effective_distance > zero(effective_distance) || throw(DomainError(
            effective_distance, "conductor-zone distance must be positive"
        ))
        logarithmic_sum += weight * log(effective_distance)
        weights += weight
    end
    return exp(logarithmic_sum / weights)
end

function geometric_mean_distance(
        left_coordinates,
        left_radii::AbstractVector,
        right_coordinates,
        right_radii::AbstractVector,
        left_element_areas::AbstractVector,
        right_element_areas::AbstractVector
)
    isempty(left_coordinates) && throw(ArgumentError(
        "the first conductor zone requires at least one element coordinate"
    ))
    isempty(right_coordinates) && throw(ArgumentError(
        "the second conductor zone requires at least one element coordinate"
    ))
    length(left_coordinates) == length(left_radii) == length(left_element_areas) ||
        throw(DimensionMismatch(
            "first conductor coordinates, radii, and areas must have equal lengths"
        ))
    length(right_coordinates) == length(right_radii) == length(right_element_areas) ||
        throw(DimensionMismatch(
            "second conductor coordinates, radii, and areas must have equal lengths"
        ))
    logarithmic_sum, weights = promote(
        zero(float(first(left_radii))),
        zero(float(first(right_radii)))
    )
    for left in eachindex(left_coordinates), right in eachindex(right_coordinates)

        weight = left_element_areas[left] * right_element_areas[right]
        weight > zero(weight) || throw(DomainError(
            weight, "conductor-zone element areas must be positive"
        ))
        distance = hypot(
            left_coordinates[left][1] - right_coordinates[right][1],
            left_coordinates[left][2] - right_coordinates[right][2]
        )
        effective_distance = iszero(distance) ?
                             max(left_radii[left], right_radii[right]) : distance
        effective_distance > zero(effective_distance) || throw(DomainError(
            effective_distance, "conductor-zone distance must be positive"
        ))
        logarithmic_sum += weight * log(effective_distance)
        weights += weight
    end
    return exp(logarithmic_sum / weights)
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

Reduce one geometrically uniform conductor zone to resistance, GMR, path, and
cross-sectional values. Primitive-specific methods preserve the implemented
disk, annulus, tape-sector, and equivalent-area rounded-sector assumptions.
"""

function conductor_zone(
        primitive::Disk,
        zone,
        ::Type{T},
        expected_inner
) where {T <: Real}
    source = first(zone)
    material = convert(Material{T}, source.source.material)
    radius = convert(T, primitive.r)
    all(item -> isapprox(
            area(item.primitive),
            (one(T) * pi) * radius^2
        ), zone) || throw(ArgumentError(
        "resolved disk strands do not preserve their declared area"
    ))
    coordinates = Tuple{T, T}[convert.(T, centroid(item.primitive))
                              for item in zone]
    centre = conductor_zone_position(zone)
    centre_radius = hypot(
        first(coordinates)[1] - centre[1],
        first(coordinates)[2] - centre[2]
    )
    circular_locus = all(
        point -> isapprox(
            hypot(point[1] - centre[1], point[2] - centre[2]),
            centre_radius
        ), coordinates)
    element_area = (one(T) * pi) * radius^2
    zone_area = length(zone) * element_area
    zone_r_in, zone_r_ex = if length(zone) == 1 && iszero(centre_radius)
        zero(T), radius
    elseif !circular_locus
        outer_radius = maximum(coordinates) do point
            hypot(point[1] - centre[1], point[2] - centre[2]) + radius
        end
        expected_outer = expected_inner + 2radius
        isapprox(outer_radius, expected_outer) || throw(ArgumentError(
            "noncircular disk strand course must extend one strand diameter " *
            "beyond the preceding course"
        ))
        expected_inner, outer_radius
    else
        declared_inner = centre_radius - radius
        isapprox(declared_inner, expected_inner) || throw(ArgumentError(
            "disk strand layer does not begin at the preceding radial boundary"
        ))
        expected_inner, expected_inner + 2radius
    end
    patterned = !isempty(source.placement.patterns)
    uniform_paths = all(item -> item.paths == source.paths, zone)
    turns = uniform_paths ? turns_per_length(source.paths, T) :
            sum(item -> turns_per_length(item.paths, T), zone) / length(zone)
    base_resistance = tubular_resistance(
        zero(T),
        radius,
        material.rho
    )
    resistance = if uniform_paths
        path_corrected_resistance(base_resistance, source.paths) / length(zone)
    else
        resistances = map(zone) do item
            path_corrected_resistance(base_resistance, item.paths)
        end
        reduce(parallel, resistances)
    end
    gmr = patterned ?
          (circular_locus ?
           strand_gmr(
        centre_radius,
        length(zone),
        radius,
        material.mu_r
    ) :
           strand_gmr(coordinates, radius, material.mu_r)) :
          tubular_gmr(zone_r_ex, zone_r_in, material.mu_r)
    return (
        r_in = zone_r_in,
        r_ex = zone_r_ex,
        area = zone_area,
        wires = patterned ? length(zone) : 0,
        turns = convert(T, turns),
        resistance,
        gmr,
        element_radius = radius,
        element_area,
        coordinates,
        position = centre,
        material,
        pairwise_gmd = any(
            entry -> entry.pattern isa Hexa,
            source.placement.patterns
        )
    )
end

function conductor_zone(
        primitive::Annulus,
        zone,
        ::Type{T},
        expected_inner
) where {T <: Real}
    length(zone) == 1 || throw(ArgumentError(
        "flatten requires one region per annular conductor layer"
    ))
    source = only(zone)
    material = convert(Material{T}, source.source.material)
    centre = (convert(T, primitive.at.x), convert(T, primitive.at.y))
    zone_r_in = convert(T, primitive.ri)
    zone_r_ex = convert(T, primitive.ro)
    zone_area = (one(T) * pi) * (zone_r_ex^2 - zone_r_in^2)
    isapprox(area(source.primitive), zone_area) || throw(ArgumentError(
        "resolved annular conductor does not preserve its declared area"
    ))
    turns = turns_per_length(source.paths, T)
    resistance = path_corrected_resistance(
        material.rho / zone_area,
        source.paths
    )
    return (
        r_in = zone_r_in,
        r_ex = zone_r_ex,
        area = zone_area,
        wires = isempty(source.paths) ? 0 : 1,
        turns = convert(T, turns),
        resistance,
        gmr = tubular_gmr(zone_r_ex, zone_r_in, material.mu_r),
        element_radius = zone_r_ex,
        element_area = zone_area,
        coordinates = Tuple{T, T}[centre],
        position = centre,
        material,
        pairwise_gmd = false
    )
end

function conductor_zone(
        primitive::Sector,
        zone,
        ::Type{T},
        expected_inner
) where {T <: Real}
    source = first(zone)
    length(source.paths) == 1 || throw(ArgumentError(
        "flatten supports a sector conductor only as one helical tape"
    ))
    length(source.placement.patterns) <= 1 || throw(ArgumentError(
        "flatten does not support nested sector repetition"
    ))
    poses = map(item -> item.primitive.at, zone)
    centre = (convert(T, first(poses).x), convert(T, first(poses).y))
    all(pose -> isapprox(pose.x, centre[1]) && isapprox(pose.y, centre[2]), poses) ||
        throw(ArgumentError(
            "a sector conductor zone requires one common centre"
        ))
    material = convert(Material{T}, source.source.material)
    zone_r_in = convert(T, primitive.ri)
    zone_r_ex = convert(T, primitive.ro)
    tape_width = convert(T, primitive.span) * (zone_r_in + zone_r_ex) / 2
    tape_thickness = isempty(source.placement.patterns) ?
                     zone_r_ex - zone_r_in :
                     (one(T) * pi) * (zone_r_ex^2 - zone_r_in^2) /
                     (length(zone) * tape_width)
    zone_area = length(zone) * tape_thickness * tape_width
    turns = turns_per_length(source.paths, T)
    resistance = path_corrected_resistance(
        strip_resistance(
            tape_thickness,
            tape_width,
            material.rho
        ),
        source.paths) / length(zone)
    return (
        r_in = zone_r_in,
        r_ex = zone_r_ex,
        area = zone_area,
        wires = isempty(source.placement.patterns) ? 1 : 0,
        turns = convert(T, turns),
        resistance,
        gmr = tubular_gmr(zone_r_ex, zone_r_in, material.mu_r),
        element_radius = zone_r_ex,
        element_area = zone_area,
        coordinates = Tuple{T, T}[centre],
        position = centre,
        material,
        pairwise_gmd = false
    )
end

function conductor_zone(
        primitive::RoundedSectorShape,
        zone,
        ::Type{T},
        expected_inner
) where {T <: Real}
    length(zone) == 1 || throw(ArgumentError(
        "flatten requires one region per rounded-sector conductor"
    ))
    source = only(zone)
    isempty(source.paths) || throw(ArgumentError(
        "flatten does not support a helical rounded-sector conductor"
    ))
    material = convert(Material{T}, source.source.material)
    zone_area = convert(T, area(primitive))
    equivalent_radius = sqrt(zone_area / (one(T) * pi))
    centre = convert.(T, centroid(primitive))
    return (
        r_in = zero(T),
        r_ex = equivalent_radius,
        area = zone_area,
        wires = 0,
        turns = zero(T),
        resistance = material.rho / zone_area,
        gmr = equivalent_radius * exp(-material.mu_r / 4),
        element_radius = equivalent_radius,
        element_area = zone_area,
        coordinates = Tuple{T, T}[centre],
        position = centre,
        material,
        pairwise_gmd = false
    )
end

function conductor_zone(
        primitive::AbstractShape,
        zone,
        ::Type,
        expected_inner
)
    throw(ArgumentError(
        "flatten does not support conductor primitive " *
        "$(nameof(typeof(primitive)))"
    ))
end

"""
$(TYPEDSIGNATURES)

Reduce one homogeneous dielectric region to concentric radii and material
properties. Primitive-specific methods retain circular and equivalent-area
rounded-sector assumptions.
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
    material = convert(Material{T}, source.source.material)
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
            <:RoundedSectorShape,
            <:RoundedSectorShape
        },
        source::PlacedRegion,
        ::Type{T}
) where {T <: Real}
    isempty(source.paths) || throw(ArgumentError(
        "flatten does not support a helical rounded-sector dielectric"
    ))
    material = convert(Material{T}, source.source.material)
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

"""
$(TYPEDSIGNATURES)

Reduce a nested repeated-wire conductor to parallel resistance, area-weighted
GMR, equivalent temperature coefficient, turns, and radial extent.
"""
function nested_conductor_zone(
        sources,
        ::Type{T},
        expected_inner::T
) where {T <: Real}
    isempty(sources) && throw(ArgumentError(
        "a nested conductor requires at least one resolved primitive"
    ))
    all(source -> source.primitive isa Disk, sources) || throw(
        ArgumentError(
        "flatten supports nested conductor paths only for disk primitives"
    )
    )

    materials = Material{T}[convert(Material{T}, source.source.material)
                            for source in sources]
    reference = first(materials).T0
    all(material -> isapprox(material.T0, reference), materials) || throw(
        ArgumentError("all nested conductor materials must share one reference temperature")
    )

    areas = T[convert(T, area(source.primitive)) for source in sources]
    radii = T[convert(T, source.primitive.r) for source in sources]
    coordinates = Tuple{T, T}[convert.(T, centroid(source.primitive))
                              for source in sources]
    total_area = sum(areas)
    centre = (
        sum(index -> areas[index] * coordinates[index][1], eachindex(sources)) /
        total_area,
        sum(index -> areas[index] * coordinates[index][2], eachindex(sources)) /
        total_area
    )
    resistances = T[path_corrected_resistance(
                        tubular_resistance(
                            zero(T),
                            radii[index],
                            materials[index].rho
                        ),
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

    weights = areas ./ total_area
    log_gmr = zero(T)
    for left in eachindex(sources)
        self_gmr = tubular_gmr(
            radii[left],
            zero(T),
            materials[left].mu_r
        )
        log_gmr += weights[left]^2 * log(self_gmr)
        for right in (left + 1):length(sources)
            distance = hypot(
                coordinates[left][1] - coordinates[right][1],
                coordinates[left][2] - coordinates[right][2]
            )
            distance > zero(distance) || throw(ArgumentError(
                "nested conductor primitives must have distinct centres"
            ))
            log_gmr += 2weights[left] * weights[right] * log(distance)
        end
    end

    turn_values = T[turns_per_length(source.paths, T) for source in sources]
    turns = sum(turn_values) / length(turn_values)
    outer = maximum(eachindex(sources)) do index
        hypot(
            coordinates[index][1] - centre[1],
            coordinates[index][2] - centre[2]
        ) + radii[index]
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

radial_position(shape::Union{Disk, Annulus, Sector}) = (shape.at.x, shape.at.y)
radial_position(shape::RoundedSectorShape) = centroid(shape)
function radial_position(
        shape::ShellShape{
        <:Any,
        <:RoundedSectorShape,
        <:RoundedSectorShape
}
)
    centroid(shape.inner)
end
radial_position(region::PlacedRegion) = radial_position(region.primitive)

function conductor_zone_position(sources)
    isempty(sources) && throw(ArgumentError(
        "a conductor zone requires at least one resolved region"
    ))
    patterned = map(sources) do source
        entries = source.placement.patterns
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

Initialize the scalar conductor reduction from its innermost zone.
"""
function initialize_conductor(zone)
    return (
        r_in = zone.r_in,
        r_ex = zone.r_ex,
        cross_section = zone.area,
        num_wires = zone.wires,
        num_turns = zone.turns,
        resistance = zone.resistance,
        alpha = zone.material.alpha,
        gmr = zone.gmr,
        reference_temperature = zone.material.T0,
        position = zone.position,
        coordinates = copy(zone.coordinates),
        element_radii = fill(zone.element_radius, length(zone.coordinates)),
        element_areas = fill(zone.element_area, length(zone.coordinates)),
        previous_coordinates = zone.coordinates,
        previous_radius = zone.element_radius,
        previous_element_area = zone.element_area,
        pairwise_gmd = zone.pairwise_gmd
    )
end

"""
$(TYPEDSIGNATURES)

Combine the strand-count-weighted turns per unit length of two conductor
zones \\[1/m\\]. A zone without explicit wires leaves the accumulated value
unchanged.
"""
function equivalent_turns(
        num_wires::Integer,
        num_turns::Real,
        added_wires::Integer,
        added_turns::Real
)
    iszero(added_wires) && return num_turns
    total_wires = num_wires + added_wires
    return (num_wires * num_turns + added_wires * added_turns) / total_wires
end

"""
$(TYPEDSIGNATURES)

Add one concentric conductor zone to an accumulated scalar reduction. The
resistances combine in parallel, while GMR and temperature coefficient are
updated from the properties of the two zones.
"""
function add_conductor_zone(conductor, zone)
    same_radial_position(zone.position, conductor.position) || throw(ArgumentError(
        "conductor zones must share one radial centre"
    ))
    isapprox(zone.r_in, conductor.r_ex) || throw(ArgumentError(
        "conductor zones must be radially contiguous"
    ))
    isapprox(zone.material.T0, conductor.reference_temperature) || throw(
        ArgumentError("all cable materials must share one reference temperature")
    )
    pairwise_gmd = conductor.pairwise_gmd || zone.pairwise_gmd
    distance = if pairwise_gmd
        geometric_mean_distance(
            conductor.coordinates,
            conductor.element_radii,
            zone.coordinates,
            fill(zone.element_radius, length(zone.coordinates)),
            conductor.element_areas,
            fill(zone.element_area, length(zone.coordinates))
        )
    else
        geometric_mean_distance(
            conductor.previous_coordinates,
            conductor.previous_radius,
            zone.coordinates,
            zone.element_radius,
            conductor.previous_element_area,
            zone.element_area
        )
    end
    return (
        r_in = conductor.r_in,
        r_ex = zone.r_ex,
        cross_section = conductor.cross_section + zone.area,
        num_wires = conductor.num_wires + zone.wires,
        num_turns = equivalent_turns(
            conductor.num_wires,
            conductor.num_turns,
            zone.wires,
            zone.turns
        ),
        resistance = parallel(
            conductor.resistance,
            zone.resistance
        ),
        alpha = equivalent_alpha(
            conductor.alpha,
            conductor.resistance,
            zone.material.alpha,
            zone.resistance
        ),
        gmr = equivalent_gmr(
            conductor.gmr,
            conductor.cross_section,
            zone.gmr,
            zone.area,
            distance
        ),
        reference_temperature = conductor.reference_temperature,
        position = conductor.position,
        coordinates = append!(copy(conductor.coordinates), zone.coordinates),
        element_radii = append!(
            copy(conductor.element_radii),
            fill(zone.element_radius, length(zone.coordinates))
        ),
        element_areas = append!(
            copy(conductor.element_areas),
            fill(zone.element_area, length(zone.coordinates))
        ),
        previous_coordinates = zone.coordinates,
        previous_radius = zone.element_radius,
        previous_element_area = zone.element_area,
        pairwise_gmd
    )
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

Calculate the shunt capacitance and conductance of one radial dielectric
layer.
"""
function dielectric_circuit(layer)
    return (
        capacitance = shunt_capacitance(
            layer.r_in,
            layer.r_ex,
            layer.material.eps_r
        ),
        conductance = shunt_conductance(
            layer.r_in,
            layer.r_ex,
            layer.material.rho
        )
    )
end

function dielectric_circuit(layers, angular_frequency::Real, ::Type{T}) where {T <: Real}
    isempty(layers) && return (
        capacitance = zero(T),
        conductance = zero(T)
    )
    combined = dielectric_circuit(first(layers))
    @inbounds for layer in Iterators.drop(layers, 1)
        circuit = dielectric_circuit(layer)
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
        layers = [(r_in = layer.r_in, r_ex = layer.r_ex, material = layer.material)]
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
    push!(layers, (r_in = layer.r_in, r_ex = layer.r_ex, material = layer.material))
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

Convert a completed scalar dielectric reduction into the homogeneous material
that reproduces its shunt capacitance and conductance. Relative permeability
retains the radial material average and the conductor's helical-solenoid
correction.
"""
function equivalent_dielectric_material(dielectric, conductor, terminal::Symbol)
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
    return Material(
        :insulator,
        inv(equivalent_conductivity(
            dielectric.shunt_conductance,
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
        zero(dielectric.r_in)
    )
end

function _radial_regions(regions)
    radial = PlacedRegion[]
    sizehint!(radial, length(regions))
    for source in regions
        source.primitive isa _DifferencePrimitive || begin
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
        isinf(material.rho) || throw(ArgumentError(
            "the coaxial reduction requires a lossless non-radial enclosure fill"
        ))
        isapprox(material.mu_r, one(material.mu_r)) || throw(ArgumentError(
            "the coaxial reduction requires a nonmagnetic non-radial enclosure fill"
        ))
        isempty(source.paths) || throw(ArgumentError(
            "a non-radial enclosure fill cannot carry a longitudinal path"
        ))
    end
    return radial
end

function same_path_definitions(left, right)
    length(left) == length(right) || return false
    return all(eachindex(left)) do index
        left[index].path == right[index].path
    end
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

        zone_ranges = UnitRange{Int}[]
        cursor = firstindex(block)
        while cursor <= lastindex(block) && block[cursor].terminal === terminal
            start = cursor
            source = block[cursor]
            signature = map(entry -> entry.pattern, source.placement.patterns)
            while cursor < lastindex(block)
                next = block[cursor + 1]
                next_signature = map(entry -> entry.pattern, next.placement.patterns)
                same_zone = next.terminal === terminal &&
                            next.source == source.source &&
                            same_path_definitions(next.paths, source.paths) &&
                            next_signature == signature
                same_zone || break
                cursor += 1
            end
            push!(zone_ranges, start:cursor)
            cursor += 1
        end
        isempty(zone_ranges) && throw(ArgumentError(
            "terminal :$terminal has no conductive zone"
        ))
        any(index -> block[index].terminal !== nothing, cursor:lastindex(block)) &&
            throw(ArgumentError(
                "terminal :$terminal reappears after its dielectric interval"
            ))

        conductor_sources = @view block[firstindex(block):(cursor - 1)]
        nested = any(conductor_sources) do source
            length(source.placement.patterns) > 1 || length(source.paths) > 1
        end
        conductor = if nested
            zone_position = conductor_zone_position(conductor_sources)
            expected_inner = if first_index > firstindex(regions) &&
                                same_radial_chain(regions[first_index - 1], zone_position)
                convert(T, r_ex(regions[first_index - 1].primitive))
            else
                zero(T)
            end
            nested_conductor_zone(
                conductor_sources,
                T,
                expected_inner
            )
        else
            accumulated = nothing
            for (zone_index, indices) in enumerate(zone_ranges)
                zone = @view block[indices]
                source = first(zone)
                source.source.material.kind === :conductor || throw(ArgumentError(
                    "terminal :$terminal contains a nonconductor physical region"
                ))
                all(item -> item.source.material == source.source.material, zone) ||
                    throw(ArgumentError("one conductor zone must use one material"))
                all(item -> item.source.primitive == source.source.primitive, zone) ||
                    throw(ArgumentError("one conductor zone must use one primitive"))
                zone_position = conductor_zone_position(zone)
                expected_inner = if zone_index > 1
                    accumulated.r_ex
                elseif first_index > firstindex(regions) &&
                       same_radial_chain(regions[first_index - 1], zone_position)
                    convert(T, r_ex(regions[first_index - 1].primitive))
                else
                    zero(T)
                end
                values = conductor_zone(
                    source.primitive,
                    zone,
                    T,
                    expected_inner
                )
                accumulated = zone_index == 1 ?
                              initialize_conductor(values) :
                              add_conductor_zone(accumulated, values)
            end
            accumulated
        end

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

This DataModel operation supports [`homogenize`](@ref). It combines the
physical radial dielectric layers in series and produces an artificial
homogeneous material that reproduces their capacitance and conductance at the
requested frequency. It does not calculate mutual coupling or earth return.

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
        dielectric = merge(component.dielectric, (
            shunt_capacitance = circuit.capacitance,
            shunt_conductance = circuit.conductance
        ))
        material = equivalent_dielectric_material(
            dielectric,
            merge(component.conductor, (
                reference_temperature = component.conductor.material.T0,
            )),
            component.name
        )
        push!(effective, (
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
material reproduces the series combination of the physical dielectric layers
at `dielectric_frequency` \\[Hz\\]. The source design and its resolved geometry
are not modified.

# Arguments

- `original`: Completed physical cable design.

# Keywords

- `new_id`: Identifier for the returned design. An empty value appends
  `"_equivalent"` to the source identifier.
- `dielectric_frequency`: Frequency used to match a lossy homogeneous
  dielectric \\[Hz\\]. Default: `50`.

# Returns

- A completed homogeneous [`CableDesign`](@ref).
"""
function homogenize(
        original::CableDesign;
        new_id::AbstractString = "",
        dielectric_frequency::Real = 50
)
    target_id = isempty(new_id) ? original.cable_id * "_equivalent" : String(new_id)
    components = flatten(
        original,
        dielectric_frequency
    )
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
                        dielectric.material
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
            _AssemblyMember(
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
