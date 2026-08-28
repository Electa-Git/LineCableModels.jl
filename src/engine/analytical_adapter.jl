constitutive(
    ::AnalyticalFormulation,
    ::Val{:conductor},
    material::Material
) = material

constitutive(
    ::AnalyticalFormulation,
    ::Val{:insulator},
    material::Material
) = material

constitutive(
    ::AnalyticalFormulation,
    ::Val{:semicon},
    material::Material
) = material

function constitutive(
        ::AnalyticalFormulation,
        ::Val{Kind},
        ::Material
) where {Kind}
    throw(ArgumentError(
        "the analytical formulation does not support material class :$Kind"
    ))
end

function _input(
        formulation::AnalyticalFormulation,
        primitive::DataModel.Disk,
        zone,
        ::Type{T},
        expected_inner
) where {T <: Real}
    source = first(zone)
    material = convert(
        Material{T},
        constitutive(formulation, Val(source.source.material.kind), source.source.material)
    )
    radius = convert(T, primitive.r)
    all(item -> isapprox(
        DataModel.area(item.shape),
        (one(T) * pi) * radius^2
    ), zone) || throw(ArgumentError(
        "resolved disk strands do not preserve their declared area"
    ))
    coordinates = Tuple{T, T}[
        convert.(T, DataModel.centroid(item.shape)) for item in zone
    ]
    centre_radius = hypot(first(coordinates)...)
    all(point -> isapprox(hypot(point...), centre_radius), coordinates) ||
        throw(ArgumentError(
            "the analytical formulation requires one circular locus per disk layer"
        ))
    element_area = (one(T) * pi) * radius^2
    zone_area = length(zone) * element_area
    zone_r_in, zone_r_ex = if length(zone) == 1 && iszero(centre_radius)
        zero(T), radius
    else
        declared_inner = centre_radius - radius
        isapprox(declared_inner, expected_inner) || throw(ArgumentError(
            "disk strand layer does not begin at the preceding radial boundary"
        ))
        expected_inner, expected_inner + 2radius
    end
    path_factor = prod(
        entry -> DataModel.overlength(entry.path, entry.radius),
        source.paths;
        init = one(T)
    )
    turns = sum(
        entry -> inv(DataModel.pitch(entry.path, entry.radius)),
        source.paths;
        init = zero(T)
    )
    patterned = !isempty(source.placement.patterns)
    resistance = DataModel.tubular_resistance(
        zero(T),
        radius,
        material.rho,
        material.alpha,
        material.T0,
        material.T0
    ) * convert(T, path_factor) / length(zone)
    gmr = patterned ? DataModel.strand_gmr(
        centre_radius,
        length(zone),
        radius,
        material.mu_r
    ) : DataModel.tubular_gmr(zone_r_ex, zone_r_in, material.mu_r)
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
        material
    )
end

_input_pose(
    ::AnalyticalFormulation,
    source::DataModel.PlacedRegion{
        <:DataModel.Region,
        <:DataModel.PlacedShape
    }
) = source.shape.at

function _input_pose(
        ::AnalyticalFormulation,
        source::DataModel.PlacedRegion
)
    throw(ArgumentError(
        "sector conductor :$(source.source.tag) requires resolved placement"
    ))
end

function _input(
        formulation::AnalyticalFormulation,
        primitive::DataModel.Annulus,
        zone,
        ::Type{T},
        expected_inner
) where {T <: Real}
    length(zone) == 1 || throw(ArgumentError(
        "the analytical formulation requires one region per annular conductor layer"
    ))
    source = only(zone)
    material = convert(
        Material{T},
        constitutive(formulation, Val(source.source.material.kind), source.source.material)
    )
    iszero(hypot(DataModel.centroid(source.shape)...)) || throw(ArgumentError(
        "the analytical formulation requires concentric annular conductors"
    ))
    zone_r_in = convert(T, primitive.ri)
    zone_r_ex = convert(T, primitive.ro)
    zone_area = (one(T) * pi) * (zone_r_ex^2 - zone_r_in^2)
    isapprox(DataModel.area(source.shape), zone_area) || throw(ArgumentError(
        "resolved annular conductor does not preserve its declared area"
    ))
    path_factor = prod(
        entry -> DataModel.overlength(entry.path, entry.radius),
        source.paths;
        init = one(T)
    )
    turns = sum(
        entry -> inv(DataModel.pitch(entry.path, entry.radius)),
        source.paths;
        init = zero(T)
    )
    resistance = material.rho * convert(T, path_factor) / zone_area
    return (
        r_in = zone_r_in,
        r_ex = zone_r_ex,
        area = zone_area,
        wires = isempty(source.paths) ? 0 : 1,
        turns = convert(T, turns),
        resistance,
        gmr = DataModel.tubular_gmr(zone_r_ex, zone_r_in, material.mu_r),
        element_radius = zone_r_ex,
        element_area = zone_area,
        coordinates = Tuple{T, T}[(zero(T), zero(T))],
        material
    )
end

function _input(
        formulation::AnalyticalFormulation,
        primitive::DataModel.Sector,
        zone,
        ::Type{T},
        expected_inner
) where {T <: Real}
    source = first(zone)
    length(source.paths) == 1 || throw(ArgumentError(
        "the analytical formulation supports a sector conductor only as one helical tape"
    ))
    length(source.placement.patterns) <= 1 || throw(ArgumentError(
        "the analytical formulation does not support nested sector repetition"
    ))
    poses = map(item -> _input_pose(formulation, item), zone)
    all(pose -> iszero(hypot(pose.x, pose.y)), poses) ||
        throw(ArgumentError(
            "the analytical formulation requires concentric sector conductors"
        ))
    material = convert(
        Material{T},
        constitutive(formulation, Val(source.source.material.kind), source.source.material)
    )
    zone_r_in = convert(T, primitive.ri)
    zone_r_ex = convert(T, primitive.ro)
    tape_width = convert(T, primitive.span) * (zone_r_in + zone_r_ex) / 2
    tape_thickness = isempty(source.placement.patterns) ?
                     zone_r_ex - zone_r_in :
                     (one(T) * pi) * (zone_r_ex^2 - zone_r_in^2) /
                     (length(zone) * tape_width)
    zone_area = length(zone) * tape_thickness * tape_width
    path_factor = prod(
        entry -> DataModel.overlength(entry.path, entry.radius),
        source.paths;
        init = one(T)
    )
    turns = sum(
        entry -> inv(DataModel.pitch(entry.path, entry.radius)),
        source.paths;
        init = zero(T)
    )
    resistance = DataModel.strip_resistance(
        tape_thickness,
        tape_width,
        material.rho,
        material.alpha,
        material.T0,
        material.T0
    ) * convert(T, path_factor) / length(zone)
    return (
        r_in = zone_r_in,
        r_ex = zone_r_ex,
        area = zone_area,
        wires = isempty(source.placement.patterns) ? 1 : 0,
        turns = convert(T, turns),
        resistance,
        gmr = DataModel.tubular_gmr(zone_r_ex, zone_r_in, material.mu_r),
        element_radius = zone_r_ex,
        element_area = zone_area,
        coordinates = Tuple{T, T}[(zero(T), zero(T))],
        material
    )
end

function _input(
        ::AnalyticalFormulation,
        primitive::DataModel.AbstractPrimitive,
        zone,
        ::Type,
        expected_inner
)
    throw(ArgumentError(
        "the analytical formulation does not support conductor primitive " *
        "$(nameof(typeof(primitive)))"
    ))
end

function _input(
        formulation::AnalyticalFormulation,
        primitive::Union{DataModel.Shell, DataModel.Annulus},
        source::DataModel.PlacedRegion,
        ::Val{:dielectric},
        ::Type{T}
) where {T <: Real}
    isempty(source.placement.patterns) || throw(ArgumentError(
        "the analytical formulation does not support repeated dielectric layers"
    ))
    isempty(source.paths) || throw(ArgumentError(
        "the analytical formulation does not support helical dielectric layers"
    ))
    iszero(hypot(DataModel.centroid(source.shape)...)) || throw(ArgumentError(
        "the analytical formulation requires concentric dielectric layers"
    ))
    material = convert(
        Material{T},
        constitutive(formulation, Val(source.source.material.kind), source.source.material)
    )
    material.kind === :conductor && throw(ArgumentError(
        "a dielectric interval cannot contain a conductor material"
    ))
    return (
        r_in = convert(T, DataModel.r_in(source.shape)),
        r_ex = convert(T, DataModel.r_ex(source.shape)),
        material
    )
end

function _input(
        ::AnalyticalFormulation,
        primitive::DataModel.AbstractPrimitive,
        source::DataModel.PlacedRegion,
        ::Val{:dielectric},
        ::Type
)
    throw(ArgumentError(
        "the analytical formulation does not support dielectric primitive " *
        "$(nameof(typeof(primitive)))"
    ))
end

"""
$(TYPEDSIGNATURES)

Adapt a completed cable design to the concentric terminal records required by
the analytical formulation.

The returned named tuples are contingent formulation input. They are rebuilt
from the authoritative physical geometry and are never stored on
[`CableDesign`](@ref).

# Arguments

- `formulation`: Analytical formulation that owns support checks and material
  interpretation.
- `design`: Completed physical cable design.
- `dielectric_frequency`: Frequency in hertz used to combine dielectric layers.

# Returns

An ordered vector of named tuples containing conductor and dielectric fields
for each retained terminal.
"""
function analytical_components(
        formulation::AnalyticalFormulation,
        design::CableDesign,
        dielectric_frequency::Real
)
    regions = design.geometry.regions
    terminals = design.terminal_order
    isempty(terminals) && throw(ArgumentError(
        "the analytical formulation requires at least one retained terminal"
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
        "the analytical formulation requires nonreappearing radial terminal blocks"
    ))

    T = promote_type(
        typeof(float(dielectric_frequency)),
        (eltype(region.shape) for region in regions)...,
        (eltype(region.source.material) for region in regions)...
    )
    frequency = convert(T, float(dielectric_frequency))
    frequency > zero(frequency) || throw(DomainError(
        frequency,
        "analytical dielectric reference frequency must be positive"
    ))
    angular_frequency = 2 * (one(T) * pi) * frequency
    effective = Any[]

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
                            next.paths == source.paths &&
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

        conductor_r_in = zero(T)
        conductor_r_ex = zero(T)
        conductor_area = zero(T)
        conductor_wires = 0
        conductor_turns = zero(T)
        conductor_resistance = zero(T)
        conductor_alpha = zero(T)
        conductor_gmr = zero(T)
        conductor_reference = zero(T)
        previous_coordinates = Tuple{T, T}[]
        previous_radius = zero(T)
        previous_element_area = zero(T)

        for (zone_index, indices) in enumerate(zone_ranges)
            zone = @view block[indices]
            source = first(zone)
            source.source.material.kind === :conductor || throw(ArgumentError(
                "terminal :$terminal contains a nonconductor physical region"
            ))
            length(source.placement.patterns) <= 1 || throw(ArgumentError(
                "the analytical formulation does not support nested conductor patterns"
            ))
            length(source.paths) <= 1 || throw(ArgumentError(
                "the analytical formulation does not support nested conductor paths"
            ))
            all(item -> item.source.material == source.source.material, zone) ||
                throw(ArgumentError("one conductor zone must use one material"))
            all(item -> item.source.primitive == source.source.primitive, zone) ||
                throw(ArgumentError("one conductor zone must use one primitive"))
            expected_inner = if zone_index > 1
                conductor_r_ex
            elseif first_index > firstindex(regions)
                convert(T, DataModel.r_ex(regions[first_index - 1].shape))
            else
                zero(T)
            end
            values = _input(
                formulation,
                source.source.primitive,
                zone,
                T,
                expected_inner
            )

            if zone_index == 1
                conductor_r_in = values.r_in
                conductor_r_ex = values.r_ex
                conductor_area = values.area
                conductor_wires = values.wires
                conductor_turns = values.turns
                conductor_resistance = values.resistance
                conductor_alpha = values.material.alpha
                conductor_gmr = values.gmr
                conductor_reference = values.material.T0
            else
                isapprox(values.r_in, conductor_r_ex) || throw(ArgumentError(
                    "conductor zones must be radially contiguous"
                ))
                isapprox(values.material.T0, conductor_reference) || throw(ArgumentError(
                    "all cable materials must share one reference temperature"
                ))
                log_sum = zero(T)
                weights = zero(T)
                weight = previous_element_area * values.element_area
                for left in previous_coordinates, right in values.coordinates
                    distance = hypot(left[1] - right[1], left[2] - right[2])
                    effective_distance = iszero(distance) ?
                                         max(previous_radius, values.element_radius) :
                                         distance
                    log_sum += weight * log(effective_distance)
                    weights += weight
                end
                distance = exp(log_sum / weights)
                beta = conductor_area / (conductor_area + values.area)
                conductor_gmr = conductor_gmr^(beta^2) *
                                values.gmr^((one(beta) - beta)^2) *
                                distance^(2 * beta * (one(beta) - beta))
                conductor_alpha = DataModel.equivalent_alpha(
                    conductor_alpha,
                    conductor_resistance,
                    values.material.alpha,
                    values.resistance
                )
                conductor_resistance = DataModel.parallel(
                    conductor_resistance,
                    values.resistance
                )
                next_wires = conductor_wires + values.wires
                conductor_turns = iszero(values.wires) ? conductor_turns :
                                  (conductor_wires * conductor_turns +
                                   values.wires * values.turns) / next_wires
                conductor_wires = next_wires
                conductor_area += values.area
                conductor_r_ex = values.r_ex
            end
            previous_coordinates = values.coordinates
            previous_radius = values.element_radius
            previous_element_area = values.element_area
        end

        dielectric_sources = cursor > lastindex(block) ?
                             @view(block[1:0]) : @view(block[cursor:lastindex(block)])
        if first_index > firstindex(regions)
            preceding = regions[first_index - 1]
            preceding.terminal === nothing || throw(ArgumentError(
                "a radial terminal must be preceded by a dielectric region"
            ))
            preceding_layer = _input(
                formulation,
                preceding.source.primitive,
                preceding,
                Val(:dielectric),
                T
            )
            conductor_r_in = preceding_layer.r_ex
        end
        if !isempty(dielectric_sources)
            first_dielectric = first(dielectric_sources)
            first_layer = _input(
                formulation,
                first_dielectric.source.primitive,
                first_dielectric,
                Val(:dielectric),
                T
            )
            conductor_r_ex = first_layer.r_in
        end

        conductor_material = Material(
            :conductor,
            DataModel.equivalent_rho(
                conductor_resistance,
                conductor_r_ex,
                conductor_r_in
            ),
            zero(T),
            DataModel.equivalent_mu(
                conductor_gmr,
                conductor_r_ex,
                conductor_r_in
            ),
            conductor_reference,
            conductor_alpha
        )

        layers = Any[]
        dielectric_r_in = conductor_r_ex
        dielectric_r_ex = conductor_r_ex
        dielectric_area = zero(T)
        dielectric_capacitance = zero(T)
        dielectric_conductance = zero(T)
        weighted_mu = zero(T)

        for (layer_index, source) in enumerate(dielectric_sources)
            source.terminal === nothing || throw(ArgumentError(
                "a dielectric region cannot own a retained terminal"
            ))
            layer = _input(
                formulation,
                source.source.primitive,
                source,
                Val(:dielectric),
                T
            )
            expected_inner = layer_index == 1 ? conductor_r_ex : dielectric_r_ex
            isapprox(layer.r_in, expected_inner) || throw(ArgumentError(
                "dielectric layers must be radially contiguous"
            ))
            isapprox(layer.material.T0, conductor_reference) || throw(ArgumentError(
                "all cable materials must share one reference temperature"
            ))
            capacitance = DataModel.shunt_capacitance(
                layer.r_in,
                layer.r_ex,
                layer.material.eps_r
            )
            conductance = DataModel.shunt_conductance(
                layer.r_in,
                layer.r_ex,
                layer.material.rho
            )
            if layer_index == 1
                dielectric_capacitance = capacitance
                dielectric_conductance = conductance
            else
                accumulated = complex(
                    dielectric_conductance,
                    angular_frequency * dielectric_capacitance
                )
                layer_value = complex(conductance, angular_frequency * capacitance)
                combined = DataModel.parallel(accumulated, layer_value)
                dielectric_conductance = real(combined)
                dielectric_capacitance = imag(combined) / angular_frequency
            end
            dielectric_area += (one(T) * pi) * (layer.r_ex^2 - layer.r_in^2)
            weighted_mu += layer.material.mu_r * log(layer.r_ex / layer.r_in)
            dielectric_r_ex = layer.r_ex
            push!(layers, (r_in = layer.r_in, r_ex = layer.r_ex, material = layer.material))
        end

        dielectric_material = if isempty(layers)
            Material(:insulator, T(Inf), zero(T), one(T), conductor_reference, zero(T))
        else
            relative_mu = weighted_mu / log(dielectric_r_ex / dielectric_r_in)
            relative_mu *= DataModel.solenoid_factor(
                conductor_turns,
                conductor_r_ex,
                dielectric_r_ex
            )
            isfinite(relative_mu) || throw(ArgumentError(
                "analytical dielectric reduction produced a non-finite permeability " *
                "for terminal :$terminal (turns=$conductor_turns, " *
                "r_con=$conductor_r_ex, r_ins=$dielectric_r_ex, " *
                "weighted_mu=$weighted_mu)"
            ))
            Material(
                :insulator,
                inv(DataModel.equivalent_conductivity(
                    dielectric_conductance,
                    dielectric_r_in,
                    dielectric_r_ex
                )),
                DataModel.equivalent_eps(
                    dielectric_capacitance,
                    dielectric_r_ex,
                    dielectric_r_in
                ),
                relative_mu,
                conductor_reference,
                zero(T)
            )
        end
        typed_layers = isempty(layers) ?
                       NamedTuple{(:r_in, :r_ex, :material),
                                  Tuple{T, T, Material{T}}}[] :
                       typeof(first(layers))[layers...]
        push!(effective, (
            name = terminal,
            conductor = (
                r_in = conductor_r_in,
                r_ex = conductor_r_ex,
                cross_section = conductor_area,
                num_wires = conductor_wires,
                num_turns = conductor_turns,
                resistance = conductor_resistance,
                alpha = conductor_alpha,
                gmr = conductor_gmr,
                material = conductor_material
            ),
            dielectric = (
                r_in = dielectric_r_in,
                r_ex = dielectric_r_ex,
                cross_section = dielectric_area,
                shunt_capacitance = dielectric_capacitance,
                shunt_conductance = dielectric_conductance,
                frequency = frequency,
                layers = typed_layers,
                material = dielectric_material
            )
        ))
    end
    return typeof(first(effective))[effective...]
end
