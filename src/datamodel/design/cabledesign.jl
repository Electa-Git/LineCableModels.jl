"""
$(TYPEDEF)

Store one authoritative physical cable declaration and its eager resolved state.

The `root` cable part is authoritative. Geometry, terminal membership, terminal
order, and the current radial analytical equivalence are derived together by
the constructor and cannot be supplied independently.

$(TYPEDFIELDS)
"""
struct CableDesign{
        T <: Real,
        R <: AbstractCablePart,
        G <: ResolvedPart,
        N,
        E
}
    "Stable cable identifier."
    cable_id::String
    "Optional catalogue data."
    nominal_data::N
    "Reference frequency used for dielectric equivalencing [Hz]."
    reference_frequency::T
    "Authoritative physical declaration."
    root::R
    "Eager resolved geometry."
    geometry::G
    "Retained terminals in physical order."
    terminal_order::Vector{Symbol}
    "Terminal index for every resolved region; zero denotes no terminal."
    terminal_map::Vector{Int}
    "Current radial analytical equivalence, or `nothing` when unsupported."
    effective::E
end

Base.eltype(::CableDesign{T}) where {T} = T
Base.eltype(::Type{<:CableDesign{T}}) where {T} = T

function _current_radial_equivalence(geometry::ResolvedPart, reference_frequency::Real)
    regions = geometry.regions
    terminals = geometry.terminals
    isempty(terminals) && return nothing

    # The current analytical engine accepts terminal blocks in increasing
    # radial order. A block starts with one or more conductive zones and may
    # end with zero or more dielectric regions. This scan is intentionally
    # here, beside CableDesign materialisation: it derives transitional solver
    # data without becoming another physical representation.
    terminal_starts = Int[]
    for terminal in terminals
        index = findfirst(region -> region.terminal === terminal, regions)
        index === nothing && return nothing
        push!(terminal_starts, index)
    end
    issorted(terminal_starts) || return nothing
    allunique(terminal_starts) || return nothing

    T = promote_type(
        typeof(float(reference_frequency)),
        (eltype(region.shape) for region in regions)...,
        (eltype(region.region.material) for region in regions)...
    )
    frequency = convert(T, float(reference_frequency))
    frequency > zero(frequency) || throw(DomainError(
        frequency,
        "reference frequency must be positive"
    ))
    angular_frequency = 2 * (one(T) * pi) * frequency
    effective = Any[]
    seen = Set{Symbol}()

    for (terminal_index, terminal) in enumerate(terminals)
        terminal in seen && return nothing
        push!(seen, terminal)
        first_index = terminal_starts[terminal_index]
        last_index = terminal_index == length(terminals) ? length(regions) :
                     terminal_starts[terminal_index + 1] - 1
        block = @view regions[first_index:last_index]

        # Conductive zones are consecutive repetitions of one Region
        # declaration. Keeping those zones explicit reproduces the existing
        # wire-layer and helical-tape equivalence without reinstating an
        # aggregate conductor-layer container as physical storage.
        zone_ranges = UnitRange{Int}[]
        cursor = firstindex(block)
        while cursor <= lastindex(block) && block[cursor].terminal === terminal
            start = cursor
            source = block[cursor]
            while cursor < lastindex(block)
                next = block[cursor + 1]
                same_zone = next.terminal === terminal &&
                            next.region == source.region &&
                            next.overlength == source.overlength &&
                            next.turns == source.turns &&
                            next.pattern_depth == source.pattern_depth &&
                            next.path_depth == source.path_depth
                same_zone || break
                cursor += 1
            end
            push!(zone_ranges, start:cursor)
            cursor += 1
        end
        isempty(zone_ranges) && return nothing
        any(index -> block[index].terminal !== nothing, cursor:lastindex(block)) &&
            return nothing

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
            source.region.material.kind === :conductor || return nothing
            source.pattern_depth <= 1 || return nothing
            source.path_depth <= 1 || return nothing
            material = convert(Material{T}, source.region.material)
            all(item -> item.region.material == source.region.material, zone) ||
                return nothing
            all(item -> item.region.primitive == source.region.primitive, zone) ||
                return nothing
            all(item -> item.overlength == source.overlength, zone) || return nothing
            all(item -> item.turns == source.turns, zone) || return nothing

            primitive = source.region.primitive
            coordinates = Tuple{T, T}[
                convert.(T, centroid(item.shape)) for item in zone
            ]
            zone_r_in = zero(T)
            zone_r_ex = zero(T)
            zone_area = zero(T)
            zone_wires = 0
            zone_turns = convert(T, source.turns)
            zone_resistance = zero(T)
            zone_gmr = zero(T)
            element_radius = zero(T)
            element_area = zero(T)

            if primitive isa Disk
                radius = convert(T, primitive.r)
                all(item -> isapprox(area(item.shape), (one(T) * pi) * radius^2), zone) ||
                    return nothing
                centre_radius = hypot(first(coordinates)...)
                all(point -> isapprox(hypot(point...), centre_radius), coordinates) ||
                    return nothing
                element_area = (one(T) * pi) * radius^2
                element_radius = radius
                zone_area = length(zone) * (one(T) * pi) * radius^2
                if length(zone) == 1 && iszero(centre_radius)
                    zone_r_in = zero(T)
                    zone_r_ex = radius
                else
                    declared_inner = centre_radius - radius
                    expected_inner = if zone_index > 1
                        conductor_r_ex
                    elseif first_index > firstindex(regions)
                        convert(T, r_ex(regions[first_index - 1].shape))
                    else
                        zero(T)
                    end
                    isapprox(declared_inner, expected_inner) || return nothing
                    zone_r_in = expected_inner
                    zone_r_ex = expected_inner + 2radius
                end
                zone_wires = source.pattern_depth > 0 ? length(zone) : 0
                zone_resistance = tubular_resistance(
                    zero(T),
                    radius,
                    material.rho,
                    material.alpha,
                    material.T0,
                    material.T0
                ) * convert(T, source.overlength) / length(zone)
                zone_gmr = source.pattern_depth > 0 ?
                           strand_gmr(centre_radius, length(zone), radius, material.mu_r) :
                           tubular_gmr(zone_r_ex, zone_r_in, material.mu_r)
            elseif primitive isa Annulus
                length(zone) == 1 || return nothing
                iszero(hypot(first(coordinates)...)) || return nothing
                zone_r_in = convert(T, primitive.ri)
                zone_r_ex = convert(T, primitive.ro)
                zone_area = (one(T) * pi) * (zone_r_ex^2 - zone_r_in^2)
                isapprox(area(source.shape), zone_area) || return nothing
                element_area = zone_area
                element_radius = zone_r_ex
                zone_wires = source.path_depth > 0 ? 1 : 0
                zone_resistance = material.rho * convert(T, source.overlength) /
                                  zone_area
                zone_gmr = tubular_gmr(zone_r_ex, zone_r_in, material.mu_r)
            elseif primitive isa Sector
                source.path_depth == 1 || return nothing
                source.pattern_depth in (0, 1) || return nothing
                all(item -> item.shape isa PlacedShape, zone) || return nothing
                all(item -> iszero(hypot(item.shape.at.x, item.shape.at.y)), zone) ||
                    return nothing
                zone_r_in = convert(T, primitive.ri)
                zone_r_ex = convert(T, primitive.ro)
                tape_width = convert(T, primitive.span) *
                             (zone_r_in + zone_r_ex) / 2
                tape_thickness = if source.pattern_depth == 0
                    zone_r_ex - zone_r_in
                else
                    (one(T) * pi) * (zone_r_ex^2 - zone_r_in^2) /
                    (length(zone) * tape_width)
                end
                zone_area = length(zone) * tape_thickness * tape_width
                element_area = zone_area
                element_radius = zone_r_ex
                coordinates = Tuple{T, T}[(zero(T), zero(T))]
                zone_wires = source.pattern_depth == 0 ? 1 : 0
                zone_resistance = strip_resistance(
                    tape_thickness,
                    tape_width,
                    material.rho,
                    material.alpha,
                    material.T0,
                    material.T0
                ) * convert(T, source.overlength) / length(zone)
                zone_gmr = tubular_gmr(zone_r_ex, zone_r_in, material.mu_r)
            else
                return nothing
            end

            if zone_index == 1
                conductor_r_in = zone_r_in
                conductor_r_ex = zone_r_ex
                conductor_area = zone_area
                conductor_wires = zone_wires
                conductor_turns = zone_turns
                conductor_resistance = zone_resistance
                conductor_alpha = material.alpha
                conductor_gmr = zone_gmr
                conductor_reference = material.T0
            else
                isapprox(zone_r_in, conductor_r_ex) || return nothing
                isapprox(material.T0, conductor_reference) || throw(ArgumentError(
                    "all cable materials must share one reference temperature"
                ))
                log_sum = zero(T)
                weights = zero(T)
                weight = previous_element_area * element_area
                for left in previous_coordinates, right in coordinates
                    distance = hypot(left[1] - right[1], left[2] - right[2])
                    effective_distance = iszero(distance) ?
                                         max(previous_radius, element_radius) : distance
                    log_sum += weight * log(effective_distance)
                    weights += weight
                end
                distance = exp(log_sum / weights)
                beta = conductor_area / (conductor_area + zone_area)
                conductor_gmr = conductor_gmr^(beta^2) *
                                zone_gmr^((one(beta) - beta)^2) *
                                distance^(2 * beta * (one(beta) - beta))
                conductor_alpha = equivalent_alpha(
                    conductor_alpha,
                    conductor_resistance,
                    material.alpha,
                    zone_resistance
                )
                conductor_resistance = parallel(
                    conductor_resistance,
                    zone_resistance
                )
                next_wires = conductor_wires + zone_wires
                conductor_turns = iszero(zone_wires) ? conductor_turns :
                                  (conductor_wires * conductor_turns +
                                   zone_wires * zone_turns) / next_wires
                conductor_wires = next_wires
                conductor_area += zone_area
                conductor_r_ex = zone_r_ex
            end
            previous_coordinates = coordinates
            previous_radius = element_radius
            previous_element_area = element_area
        end

        dielectric_sources = cursor > lastindex(block) ?
                             @view(block[1:0]) : @view(block[cursor:lastindex(block)])
        if first_index > firstindex(regions)
            preceding = regions[first_index - 1]
            preceding.terminal === nothing || return nothing
            preceding.region.primitive isa Union{Shell, Annulus} || return nothing
            iszero(hypot(centroid(preceding.shape)...)) || return nothing
            conductor_r_in = convert(T, r_ex(preceding.shape))
        end
        if !isempty(dielectric_sources)
            first_dielectric = first(dielectric_sources)
            first_dielectric.region.primitive isa Union{Shell, Annulus} || return nothing
            iszero(hypot(centroid(first_dielectric.shape)...)) || return nothing
            conductor_r_ex = convert(T, r_in(first_dielectric.shape))
        end

        conductor_material = Material(
            :conductor,
            equivalent_rho(
                conductor_resistance,
                conductor_r_ex,
                conductor_r_in
            ),
            zero(T),
            equivalent_mu(conductor_gmr, conductor_r_ex, conductor_r_in),
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
            source.terminal === nothing || return nothing
            source.region.material.kind === :conductor && return nothing
            source.pattern_depth == 0 || return nothing
            source.path_depth == 0 || return nothing
            primitive = source.region.primitive
            primitive isa Union{Shell, Annulus} || return nothing
            iszero(hypot(centroid(source.shape)...)) || return nothing
            layer_r_in = convert(T, r_in(source.shape))
            layer_r_ex = convert(T, r_ex(source.shape))
            expected_inner = layer_index == 1 ? conductor_r_ex : dielectric_r_ex
            isapprox(layer_r_in, expected_inner) || return nothing
            material = convert(Material{T}, source.region.material)
            isapprox(material.T0, conductor_reference) || throw(ArgumentError(
                "all cable materials must share one reference temperature"
            ))
            capacitance = shunt_capacitance(
                layer_r_in,
                layer_r_ex,
                material.eps_r
            )
            conductance = shunt_conductance(
                layer_r_in,
                layer_r_ex,
                material.rho
            )
            if layer_index == 1
                dielectric_capacitance = capacitance
                dielectric_conductance = conductance
            else
                accumulated = complex(
                    dielectric_conductance,
                    angular_frequency * dielectric_capacitance
                )
                layer_value = complex(
                    conductance,
                    angular_frequency * capacitance
                )
                combined = parallel(accumulated, layer_value)
                dielectric_conductance = real(combined)
                dielectric_capacitance = imag(combined) / angular_frequency
            end
            dielectric_area += (one(T) * pi) *
                               (layer_r_ex^2 - layer_r_in^2)
            weighted_mu += material.mu_r * log(layer_r_ex / layer_r_in)
            dielectric_r_ex = layer_r_ex
            push!(layers, (
                r_in = layer_r_in,
                r_ex = layer_r_ex,
                material = material
            ))
        end

        dielectric_material = if isempty(layers)
            Material(:insulator, T(Inf), zero(T), one(T), conductor_reference, zero(T))
        else
            relative_mu = weighted_mu / log(dielectric_r_ex / dielectric_r_in)
            relative_mu *= solenoid_factor(
                conductor_turns,
                conductor_r_ex,
                dielectric_r_ex
            )
            Material(
                :insulator,
                inv(equivalent_conductivity(
                    dielectric_conductance,
                    dielectric_r_in,
                    dielectric_r_ex
                )),
                equivalent_eps(
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
                reference_frequency = frequency,
                layers = typed_layers,
                material = dielectric_material
            )
        ))
    end
    return typeof(first(effective))[effective...]
end

function CableDesign(
        root::AbstractCablePart;
        cable_id::AbstractString = "cable",
        nominal_data::Union{Nothing, NominalData} = nothing,
        reference_frequency::Real = 50
)
    # 1. Normalize the unambiguous single-conductor shorthand.
    normalized = if root isa Region && root.material.kind === :conductor
        Group(root.tag, Pose2(0, 0, 0), root, nothing, nothing, nothing)
    else
        root
    end

    # 2. Check formulation-independent physical invariants.
    identifier = String(cable_id)
    isempty(identifier) && throw(ArgumentError("cable_id cannot be empty"))
    isfinite(reference_frequency) && reference_frequency > zero(reference_frequency) ||
        throw(DomainError(reference_frequency, "reference frequency must be positive"))

    # 3–4. Resolve contextual geometry, placement, path, and compaction through
    # the one physical resolution action owned by the cable-part grammar.
    geometry = resolve(EmptyBoundary(), normalized)
    isempty(geometry.regions) && throw(ArgumentError("a cable design requires one region"))

    # 5. Resolve retained terminal membership and order.
    terminal_order = copy(geometry.terminals)
    isempty(terminal_order) && throw(ArgumentError(
        "a cable design requires at least one retained terminal"
    ))
    allunique(terminal_order) || throw(ArgumentError(
        "retained terminal names must be unique"
    ))
    terminal_map = Int[]
    for source in geometry.regions
        if source.region.material.kind === :conductor && source.terminal === nothing
            throw(ArgumentError(
                "conductive region :$(source.region.tag) is not owned by a Group terminal"
            ))
        end
        index = source.terminal === nothing ? 0 :
                something(findfirst(==(source.terminal), terminal_order), 0)
        source.terminal === nothing || index > 0 || throw(ArgumentError(
            "resolved region references an unknown terminal :$(source.terminal)"
        ))
        push!(terminal_map, index)
    end

    reference = first(geometry.regions).region.material.T0
    all(source -> isapprox(source.region.material.T0, reference), geometry.regions) ||
        throw(ArgumentError("all cable materials must share one reference temperature"))

    # 6. Derive the current radial engine equivalence when the resolved cable
    # satisfies that formulation's narrow compatibility profile.
    effective = _current_radial_equivalence(geometry, reference_frequency)

    # 7. Freeze authoritative and derived state together.
    T = promote_type(
        (eltype(source.shape) for source in geometry.regions)...,
        (eltype(source.region.material) for source in geometry.regions)...,
        typeof(float(reference_frequency))
    )
    nominal = nominal_data === nothing ? nothing : convert(NominalData{T}, nominal_data)
    return CableDesign{
        T,
        typeof(normalized),
        typeof(geometry),
        typeof(nominal),
        typeof(effective)
    }(
        identifier,
        nominal,
        convert(T, float(reference_frequency)),
        normalized,
        geometry,
        terminal_order,
        terminal_map,
        effective
    )
end

function Base.show(io::IO, ::MIME"text/plain", design::CableDesign)
    print(
        io,
        "CableDesign \"$(design.cable_id)\": [regions=$(length(design.geometry.regions)), " *
        "terminals=($(join(design.terminal_order, ", "))), " *
        "outer_radius=$(round(outer_radius(design), sigdigits=5))]"
    )
end
