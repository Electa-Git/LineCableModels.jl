constitutive(
    ::LineParametersFormulation,
    ::Val{:conductor},
    material::Material
) = material

constitutive(
    ::LineParametersFormulation,
    ::Val{:insulator},
    material::Material
) = material

constitutive(
    ::LineParametersFormulation,
    ::Val{:semicon},
    material::Material
) = material

function constitutive(
        ::LineParametersFormulation,
        ::Val{Kind},
        ::Material
) where {Kind}
    throw(ArgumentError(
        "the native formulation does not support material class :$Kind"
    ))
end

function _path_resistance(
        ::LineParametersFormulation, resistance::T, ::Tuple{}
) where {T}
    resistance
end

function _path_resistance(
        formulation::LineParametersFormulation,
        resistance::T,
        paths::Tuple
) where {T}
    inner = _path_resistance(formulation, resistance, Base.front(paths))
    entry = last(paths)
    return inner * convert(T, DataModel.overlength(entry.path, entry.radius))
end

function _input(
        formulation::LineParametersFormulation,
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
            DataModel.area(item.primitive),
            (one(T) * pi) * radius^2
        ), zone) || throw(ArgumentError(
        "resolved disk strands do not preserve their declared area"
    ))
    coordinates = Tuple{T, T}[convert.(T, DataModel.centroid(item.primitive))
                              for item in zone]
    centre_radius = hypot(first(coordinates)...)
    all(point -> isapprox(hypot(point...), centre_radius), coordinates) ||
        throw(ArgumentError(
            "the native formulation requires one circular locus per disk layer"
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
    turns = sum(
        entry -> inv(DataModel.pitch(entry.path, entry.radius)),
        source.paths;
        init = zero(T)
    )
    patterned = !isempty(source.placement.patterns)
    resistance = _path_resistance(formulation,
        DataModel.tubular_resistance(
            zero(T),
            radius,
            material.rho,
            material.alpha,
            material.T0,
            material.T0
        ),
        source.paths) / length(zone)
    gmr = patterned ?
          DataModel.strand_gmr(
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

function _input(
        formulation::LineParametersFormulation,
        primitive::DataModel.Annulus,
        zone,
        ::Type{T},
        expected_inner
) where {T <: Real}
    length(zone) == 1 || throw(ArgumentError(
        "the native formulation requires one region per annular conductor layer"
    ))
    source = only(zone)
    material = convert(
        Material{T},
        constitutive(formulation, Val(source.source.material.kind), source.source.material)
    )
    iszero(hypot(DataModel.centroid(source.primitive)...)) || throw(ArgumentError(
        "the native formulation requires concentric annular conductors"
    ))
    zone_r_in = convert(T, primitive.ri)
    zone_r_ex = convert(T, primitive.ro)
    zone_area = (one(T) * pi) * (zone_r_ex^2 - zone_r_in^2)
    isapprox(DataModel.area(source.primitive), zone_area) || throw(ArgumentError(
        "resolved annular conductor does not preserve its declared area"
    ))
    turns = sum(
        entry -> inv(DataModel.pitch(entry.path, entry.radius)),
        source.paths;
        init = zero(T)
    )
    resistance = _path_resistance(
        formulation,
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
        gmr = DataModel.tubular_gmr(zone_r_ex, zone_r_in, material.mu_r),
        element_radius = zone_r_ex,
        element_area = zone_area,
        coordinates = Tuple{T, T}[(zero(T), zero(T))],
        material
    )
end

function _input(
        formulation::LineParametersFormulation,
        primitive::DataModel.Sector,
        zone,
        ::Type{T},
        expected_inner
) where {T <: Real}
    source = first(zone)
    length(source.paths) == 1 || throw(ArgumentError(
        "the native formulation supports a sector conductor only as one helical tape"
    ))
    length(source.placement.patterns) <= 1 || throw(ArgumentError(
        "the native formulation does not support nested sector repetition"
    ))
    poses = map(item -> item.primitive.at, zone)
    all(pose -> iszero(hypot(pose.x, pose.y)), poses) ||
        throw(ArgumentError(
            "the native formulation requires concentric sector conductors"
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
    turns = sum(
        entry -> inv(DataModel.pitch(entry.path, entry.radius)),
        source.paths;
        init = zero(T)
    )
    resistance = _path_resistance(formulation,
        DataModel.strip_resistance(
            tape_thickness,
            tape_width,
            material.rho,
            material.alpha,
            material.T0,
            material.T0
        ),
        source.paths) / length(zone)
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
        formulation::LineParametersFormulation,
        primitive::DataModel.RoundedSectorShape,
        zone,
        ::Type{T},
        expected_inner
) where {T <: Real}
    length(zone) == 1 || throw(ArgumentError(
        "the native formulation requires one region per rounded-sector conductor"
    ))
    source = only(zone)
    isempty(source.paths) || throw(ArgumentError(
        "the native formulation does not support a helical rounded-sector conductor"
    ))
    material = convert(
        Material{T},
        constitutive(formulation, Val(source.source.material.kind), source.source.material)
    )
    zone_area = convert(T, DataModel.area(primitive))
    equivalent_radius = sqrt(zone_area / (one(T) * pi))
    centre = convert.(T, DataModel.centroid(primitive))
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
        material
    )
end

function _input(
        ::LineParametersFormulation,
        primitive::DataModel.AbstractPrimitive,
        zone,
        ::Type,
        expected_inner
)
    throw(ArgumentError(
        "the native formulation does not support conductor primitive " *
        "$(nameof(typeof(primitive)))"
    ))
end

function _input(
        formulation::LineParametersFormulation,
        primitive::DataModel.Annulus,
        source::DataModel.PlacedRegion,
        ::Val{:dielectric},
        ::Type{T}
) where {T <: Real}
    isempty(source.placement.patterns) || throw(ArgumentError(
        "the native formulation does not support repeated dielectric layers"
    ))
    isempty(source.paths) || throw(ArgumentError(
        "the native formulation does not support helical dielectric layers"
    ))
    iszero(hypot(DataModel.centroid(source.primitive)...)) || throw(ArgumentError(
        "the native formulation requires concentric dielectric layers"
    ))
    material = convert(
        Material{T},
        constitutive(formulation, Val(source.source.material.kind), source.source.material)
    )
    material.kind === :conductor && throw(ArgumentError(
        "a dielectric interval cannot contain a conductor material"
    ))
    return (
        r_in = convert(T, DataModel.r_in(source.primitive)),
        r_ex = convert(T, DataModel.r_ex(source.primitive)),
        material
    )
end

function _input(
        formulation::LineParametersFormulation,
        primitive::DataModel.ShellShape{
            <:Any,
            <:DataModel.RoundedSectorShape,
            <:DataModel.RoundedSectorShape
        },
        source::DataModel.PlacedRegion,
        ::Val{:dielectric},
        ::Type{T}
) where {T <: Real}
    isempty(source.paths) || throw(ArgumentError(
        "the native formulation does not support a helical rounded-sector dielectric"
    ))
    material = convert(
        Material{T},
        constitutive(formulation, Val(source.source.material.kind), source.source.material)
    )
    material.kind === :conductor && throw(ArgumentError(
        "a dielectric interval cannot contain a conductor material"
    ))
    inner_area = convert(T, DataModel.area(primitive.inner))
    outer_area = convert(T, DataModel.area(primitive.outer))
    inner_radius = sqrt(inner_area / (one(T) * pi))
    outer_radius = sqrt(outer_area / (one(T) * pi))
    return (
        r_in = inner_radius,
        r_ex = outer_radius,
        material
    )
end

function _input(
        ::LineParametersFormulation,
        primitive::DataModel.AbstractPrimitive,
        source::DataModel.PlacedRegion,
        ::Val{:dielectric},
        ::Type
)
    throw(ArgumentError(
        "the native formulation does not support dielectric primitive " *
        "$(nameof(typeof(primitive)))"
    ))
end

function _nested_conductor_input(
        formulation::LineParametersFormulation,
        sources,
        ::Type{T},
        expected_inner::T
) where {T <: Real}
    isempty(sources) && throw(ArgumentError(
        "a nested conductor requires at least one resolved primitive"
    ))
    all(source -> source.primitive isa DataModel.Disk, sources) || throw(
        ArgumentError(
        "the native formulation supports nested conductor paths only for disk primitives"
    )
    )

    materials = Material{T}[convert(
                                Material{T},
                                constitutive(
                                    formulation,
                                    Val(source.source.material.kind),
                                    source.source.material
                                )
                            ) for source in sources]
    reference = first(materials).T0
    all(material -> isapprox(material.T0, reference), materials) || throw(
        ArgumentError("all nested conductor materials must share one reference temperature")
    )

    areas = T[convert(T, DataModel.area(source.primitive)) for source in sources]
    radii = T[convert(T, source.primitive.r) for source in sources]
    coordinates = Tuple{T, T}[convert.(T, DataModel.centroid(source.primitive))
                              for source in sources]
    resistances = T[_path_resistance(
                        formulation,
                        DataModel.tubular_resistance(
                            zero(T),
                            radii[index],
                            materials[index].rho,
                            materials[index].alpha,
                            materials[index].T0,
                            materials[index].T0
                        ),
                        sources[index].paths
                    ) for index in eachindex(sources)]

    resistance = first(resistances)
    alpha = first(materials).alpha
    for index in Iterators.drop(eachindex(sources), 1)
        alpha = DataModel.equivalent_alpha(
            alpha,
            resistance,
            materials[index].alpha,
            resistances[index]
        )
        resistance = DataModel.parallel(resistance, resistances[index])
    end

    total_area = sum(areas)
    weights = areas ./ total_area
    log_gmr = zero(T)
    for left in eachindex(sources)
        self_gmr = DataModel.tubular_gmr(
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

    turn_values = T[sum(
                        entry -> inv(DataModel.pitch(entry.path, entry.radius)),
                        source.paths;
                        init = zero(T)
                    ) for source in sources]
    turns = sum(turn_values) / length(turn_values)
    outer = maximum(source -> convert(T, DataModel.support(source.primitive)), sources)
    outer > expected_inner || throw(ArgumentError(
        "nested conductor envelope must exceed its preceding radial boundary"
    ))
    return (
        r_in = expected_inner,
        r_ex = outer,
        area = total_area,
        wires = length(sources),
        turns,
        resistance,
        alpha,
        gmr = exp(log_gmr),
        reference,
        position = (zero(T), zero(T))
    )
end

_radial_chain_key(::DataModel.AbstractShape, placement) = nothing
_radial_chain_key(::DataModel.RoundedSectorShape, placement) = placement.patterns
_radial_chain_key(
    ::DataModel.ShellShape{
        <:Any,
        <:DataModel.RoundedSectorShape,
        <:DataModel.RoundedSectorShape
    },
    placement
) = placement.patterns

function _same_radial_chain(left::DataModel.PlacedRegion, right::DataModel.PlacedRegion)
    left_key = _radial_chain_key(left.primitive, left.placement)
    right_key = _radial_chain_key(right.primitive, right.placement)
    return left_key === nothing && right_key === nothing || left_key == right_key
end

"""
$(TYPEDSIGNATURES)

Adapt a completed cable design to the concentric terminal records required by
the native formulation.

The returned named tuples are contingent formulation input. They are rebuilt
from the authoritative physical geometry and are never stored on
[`CableDesign`](@ref).

# Arguments

- `formulation`: Native formulation that owns support checks and material
  interpretation.
- `design`: Completed physical cable design.
- `dielectric_frequency`: Frequency in hertz used to combine dielectric layers.

# Returns

An ordered vector of named tuples containing conductor and dielectric fields
for each retained terminal.
"""
function homogeneous_components(
        formulation::LineParametersFormulation,
        design::CableDesign,
        dielectric_frequency::Real
)
    T = promote_type(
        typeof(float(dielectric_frequency)),
        (eltype(region.primitive) for region in design.geometry.regions)...,
        (eltype(region.source.material) for region in design.geometry.regions)...
    )
    return homogeneous_components(
        formulation,
        design,
        convert(T, float(dielectric_frequency)),
        T
    )
end

function homogeneous_components(
        formulation::LineParametersFormulation,
        design::CableDesign,
        frequency::T,
        ::Type{T}
) where {T <: Real}
    regions = design.geometry.regions
    terminals = design.terminal_order
    isempty(terminals) && throw(ArgumentError(
        "the native formulation requires at least one retained terminal"
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
        "the native formulation requires nonreappearing radial terminal blocks"
    ))

    frequency > zero(frequency) || throw(DomainError(
        frequency,
        "native dielectric reference frequency must be positive"
    ))
    angular_frequency = 2 * (one(T) * pi) * frequency
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
        (:r_in, :r_ex, :cross_section, :shunt_capacitance,
            :shunt_conductance, :frequency, :layers, :material),
        Tuple{T, T, T, T, T, T, Vector{Layer}, Material{T}}
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
        conductor_position = (zero(T), zero(T))
        previous_coordinates = Tuple{T, T}[]
        previous_radius = zero(T)
        previous_element_area = zero(T)

        conductor_sources = @view block[firstindex(block):(cursor - 1)]
        nested = any(conductor_sources) do source
            length(source.placement.patterns) > 1 || length(source.paths) > 1
        end
        if nested
            expected_inner = if first_index > firstindex(regions) &&
                    _same_radial_chain(regions[first_index - 1], first(conductor_sources))
                convert(T, DataModel.r_ex(regions[first_index - 1].primitive))
            else
                zero(T)
            end
            values = _nested_conductor_input(
                formulation,
                conductor_sources,
                T,
                expected_inner
            )
            conductor_r_in = values.r_in
            conductor_r_ex = values.r_ex
            conductor_area = values.area
            conductor_wires = values.wires
            conductor_turns = values.turns
            conductor_resistance = values.resistance
            conductor_alpha = values.alpha
            conductor_gmr = values.gmr
            conductor_reference = values.reference
            conductor_position = values.position
        else
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
                expected_inner = if zone_index > 1
                    conductor_r_ex
                elseif first_index > firstindex(regions) &&
                        _same_radial_chain(regions[first_index - 1], source)
                    convert(T, DataModel.r_ex(regions[first_index - 1].primitive))
                else
                    zero(T)
                end
                values = _input(
                    formulation,
                    source.primitive,
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
                    conductor_position = hasproperty(values, :position) ?
                                         values.position : (zero(T), zero(T))
                else
                    isapprox(values.r_in, conductor_r_ex) || throw(ArgumentError(
                        "conductor zones must be radially contiguous"
                    ))
                    isapprox(values.material.T0, conductor_reference) ||
                        throw(ArgumentError(
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
        end

        dielectric_sources = cursor > lastindex(block) ?
                             @view(block[1:0]) : @view(block[cursor:lastindex(block)])
        if first_index > firstindex(regions) &&
                _same_radial_chain(regions[first_index - 1], first(block))
            preceding = regions[first_index - 1]
            preceding.terminal === nothing || throw(ArgumentError(
                "a radial terminal must be preceded by a dielectric region"
            ))
            preceding_layer = _input(
                formulation,
                preceding.primitive,
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
                first_dielectric.primitive,
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

        layers = Layer[]
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
                source.primitive,
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
                "native dielectric reduction produced a non-finite permeability " *
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
        push!(effective,
            (
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
                    position = conductor_position,
                    material = conductor_material
                ),
                dielectric = (
                    r_in = dielectric_r_in,
                    r_ex = dielectric_r_ex,
                    cross_section = dielectric_area,
                    shunt_capacitance = dielectric_capacitance,
                    shunt_conductance = dielectric_conductance,
                    frequency = frequency,
                    layers,
                    material = dielectric_material
                )
            ))
    end
    return effective
end

"""
$(TYPEDSIGNATURES)

Build a homogeneous cable design from the effective radial components admitted
by the native formulation.

Each retained terminal becomes one solid or annular conductor followed by one
homogeneous dielectric interval. The conductor material reproduces the
terminal's parallel resistance and geometric-mean radius. The dielectric
material reproduces the series combination of the physical dielectric layers
at `dielectric_frequency` \\[Hz\\]. The source design and its resolved geometry
are not modified.

# Arguments

- `formulation`: Native formulation that owns geometry and material support.
- `original`: Completed physical cable design.

# Keywords

- `new_id`: Identifier for the returned design. An empty value appends
  `"_equivalent"` to the source identifier.
- `dielectric_frequency`: Frequency used to match a lossy homogeneous
  dielectric \\[Hz\\]. Default: `50`.

# Returns

- A completed homogeneous [`CableDesign`](@ref).
"""
function flatten(
        formulation::LineParametersFormulation,
        original::CableDesign;
        new_id::AbstractString = "",
        dielectric_frequency::Real = 50
)
    target_id = isempty(new_id) ? original.cable_id * "_equivalent" : String(new_id)
    components = homogeneous_components(
        formulation,
        original,
        dielectric_frequency
    )
    parts = DataModel.AbstractCablePart[]
    for component in components
        terminal = component.name
        conductor = component.conductor
        conductor_primitive = iszero(conductor.r_in) ?
                              DataModel.Disk(conductor.r_ex) :
                              DataModel.Annulus(conductor.r_in, conductor.r_ex)
        conductor_region = DataModel.Region(
            Symbol(terminal, :_equivalent_conductor),
            conductor_primitive,
            conductor.material
        )
        push!(parts,
            DataModel.Group(
                terminal,
                DataModel.Pose2(0, 0, 0),
                conductor_region,
                nothing,
                nothing,
                nothing
            ))

        dielectric = component.dielectric
        if dielectric.r_ex > dielectric.r_in
            push!(parts,
                DataModel.Region(
                    Symbol(terminal, :_equivalent_dielectric),
                    DataModel.Annulus(dielectric.r_in, dielectric.r_ex),
                    dielectric.material
                ))
        end
    end
    return build(
        CableDesign,
        target_id,
        DataModel.Stack(parts)
    )
end

function flatten(
        original::CableDesign;
        new_id::AbstractString = "",
        dielectric_frequency::Real = 50
)
    return flatten(
        Formulation(),
        original;
        new_id,
        dielectric_frequency
    )
end
