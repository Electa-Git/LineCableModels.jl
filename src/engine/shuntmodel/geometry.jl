# Dielectric group records owned by the local shunt formulations.
const ShuntWire{T} = NamedTuple{(:x, :y, :r, :terminal), Tuple{T, T, T, Int}}
const ShuntTape{T} = NamedTuple{(:ri, :ro, :phi, :span, :terminal), Tuple{T, T, T, T, Int}}
const ShuntLayer{T} = NamedTuple{(:ri, :ro, :epsilon), Tuple{T, T, T}}
# A host that the formula cannot represent, with its terminals and the reason.
const ShuntRejection = NamedTuple{(:terminals, :reason), Tuple{Vector{Int}, BoundarySolveError}}

"""
$(TYPEDEF)

Hold the dielectric layers between two closed equipotentials, the inner anchor and the
shield or jacket reference, with the open conductors inside them, solved as one group. The
first terminal is the inner anchor. A shield occupies the last real terminal, while a jacket
reference has a radius and does not own a terminal. Coordinates and radii use \\[m\\]. Relative
permittivities are those of the formulation's lossless laws. The extraction builds each
group once, and the boundary shunt calculation consumes it.

$(TYPEDFIELDS)
"""
struct DielectricGroup{T <: Real}
    "Source design index."
    design::Int
    "Global conductor range of the containing concentric assembly."
    assembly::UnitRange{Int}
    "Real conductors of the group within its assembly."
    terminals::UnitRange{Int}
    "Closed metallic shield or external jacket equipotential."
    reference::Union{ShuntShield{T}, ShuntJacket{T}}
    "Host dielectric inner radius \\[m\\]."
    a::T
    "Host dielectric outer radius \\[m\\]."
    b::T
    "Host relative permittivity \\[dimensionless\\]."
    epsilon::T
    "Ordered layers between inner conductor and host."
    left::Vector{ShuntLayer{T}}
    "Ordered layers between host and reference shield."
    right::Vector{ShuntLayer{T}}
    "Exposed circular wires in the host frame."
    wires::Vector{ShuntWire{T}}
    "Exposed finite-thickness tapes in the host frame."
    tapes::Vector{ShuntTape{T}}
end

_shunt_scalar(x) = float(nominal(x))
# Faces touch when their gap is below this tolerance, which scales with the host's outer
# radius.
_shunt_tol(x) = 128eps(typeof(_shunt_scalar(x))) * abs(_shunt_scalar(x))
_shunt_same(a, b, scale) = tolerance_norm(a - b) <= _shunt_tol(scale)
function _shunt_concentric(shape, x, y, scale)
    _shunt_same(shape.at.x, x, scale) && _shunt_same(shape.at.y, y, scale)
end

# Relative permittivity of a material under the formulation's lossless laws.
function _shunt_permittivity(material::Material{T}, formulas) where {T}
    frequency = one(T)
    kappa = constitutive((formulas.insulation_admittance, formulas.semicon_admittance),
        material, frequency, material.T0; temperature_dependence = nothing)
    return imag(kappa)/(2pi*frequency*vacuum_permittivity(T))
end

# Preserve ordering and dependency information. Nominal values decide ordering.
# Coincident interfaces and centers must also have matching uncertain dependencies.
function _shunt_layers(regions, inner, outer, x, y, ::Type{T}, formulas) where {T}
    inner, outer = T(inner), T(outer)
    layers = ShuntLayer{T}[]
    _shunt_same(inner, outer, outer) && return layers
    for placed in regions
        shape = placed.primitive
        placed.terminal === nothing && shape isa DataModel.Annulus || continue
        _shunt_concentric(shape, x, y, outer) || continue
        _shunt_scalar(shape.ri - inner) >= -_shunt_tol(outer) &&
        _shunt_scalar(shape.ro - outer) <= _shunt_tol(outer) || continue
        push!(layers,
            (ri = T(shape.ri), ro = T(shape.ro), epsilon = _shunt_permittivity(
                convert(Material{T}, placed.source.material), formulas)))
    end
    sort!(layers; by = layer -> _shunt_scalar(layer.ri))
    cursor = inner
    for layer in layers
        _shunt_same(cursor, layer.ri, outer) || return nothing
        cursor = layer.ro
    end
    _shunt_same(cursor, outer, outer) || return nothing
    return layers
end

function _shunt_closed(design, index, x, y, scale)
    indices = findall(==(index), design.terminal_map)
    length(indices) == 1 || return nothing
    shape = design.geometry.regions[only(indices)].primitive
    shape isa DataModel.Annulus && _shunt_concentric(shape, x, y, scale) ||
        return nothing
    return shape
end

# Add the group of one annular host to `groups`, or the reason that the boundary shunt
# formula cannot represent it to `rejected`. A host without terminal members is not a group.
function _shunt_host_group!(groups, rejected, design, blueprint, host_region, design_index,
        ::Type{T}, formulas) where {T}
    host = host_region.primitive.outer
    x, y, phi = host.at.x, host.at.y, host.at.φ
    scale = T(host.ro)
    members = Int[]
    uncertain_contact = false
    for (i, placed) in pairs(design.geometry.regions)
        placed.terminal === nothing && continue
        shape = placed.primitive
        inner, outer = if shape isa DataModel.Disk
            radius = hypot(shape.at.x - x, shape.at.y - y)
            radius-shape.r, radius+shape.r
        elseif shape isa DataModel.BentStrip && _shunt_concentric(shape, x, y, scale)
            shape.ri, shape.ro
        else
            continue
        end
        _shunt_scalar(inner - host.ri) >= -_shunt_tol(scale) &&
            _shunt_scalar(outer - host.ro) <= _shunt_tol(scale) || continue
        # A nominal contact with an independent uncertain displacement
        # crosses a dielectric interface, outside a fixed-topology tangent.
        for (edge, boundary) in ((inner, host.ri), (outer, host.ro))
            uncertain_contact |= abs(_shunt_scalar(edge-boundary)) <= _shunt_tol(scale) &&
                                 !_shunt_same(edge, boundary, scale)
        end
        push!(members, i)
    end
    isempty(members) && return groups
    terminals = sort!(unique(design.terminal_map[members]))
    reject(reason) = push!(rejected, (; terminals, reason = BoundarySolveError(:unsupported,
        (; cable = design.cable_id, terminals), "boundary shunt cannot represent this host: " * reason)))
    uncertain_contact && return reject(
        "an uncertain displacement changes whether a wire or tape touches the host boundary")
    # Coatings and dielectric islands invalidate the homogeneous-host Green function.
    holes = host_region.primitive.holes
    length(holes) == length(members) && all(
        index -> any(hole -> isequal(design.geometry.regions[index].primitive, hole), holes), members) ||
        return reject("it holds coatings or dielectric islands besides its exposed wires and tapes")
    first(terminals) > 1 || return reject(
        "it holds pieces of the innermost conductor, which keeps the equivalent-core assumption")
    terminals == collect(first(terminals):last(terminals)) ||
        return reject("its terminals are not consecutive in radial order")
    # Never resolve only some pieces of a retained terminal.
    all(i -> design.terminal_map[i] ∉ terminals || i in members,
        eachindex(design.terminal_map)) || return reject("a terminal has pieces outside it")
    first_index = first(terminals) - 1
    assembly_index = blueprint.conductors[first_index].assembly
    assembly = blueprint.assembly_ranges[assembly_index]
    last(terminals) in assembly || return reject("its terminals span several assemblies")
    inner = blueprint.conductors[first_index]
    center = blueprint.positions[assembly_index]
    _shunt_same(center[1], x, scale) && _shunt_same(center[2], y, scale) ||
        return reject("it is not concentric with its assembly")
    # Inner cores keep the equivalent-core assumption. Other anchors must be closed.
    first_index == first(assembly) ||
        _shunt_closed(design, first_index, x, y, scale) !== nothing ||
        return reject("the conductor inside it is not a closed concentric tube")
    reference, last_index = if last(terminals) < last(assembly)
        index = last(terminals)+1
        shield = _shunt_closed(design, index, x, y, scale)
        shield === nothing &&
            return reject("the next conductor outward is not a closed concentric shield")
        ShuntShield{T}(shield.ri), index
    else
        layers = blueprint[last(assembly)].layers
        isempty(layers) && return reject("the outermost conductor has no jacket to reference")
        ShuntJacket{T}(last(layers).r_ex), last(terminals)
    end
    left = _shunt_layers(design.geometry.regions, inner.r_ex, host.ri, x, y, T, formulas)
    left === nothing && return reject(
        "the dielectrics between the inner anchor and the host are not continuous concentric layers")
    right = _shunt_layers(design.geometry.regions, host.ro, reference.radius, x, y, T,
        formulas)
    right === nothing && return reject(
        "the dielectrics between the host and its reference are not continuous concentric layers")
    wires, tapes = ShuntWire{T}[], ShuntTape{T}[]
    for index in members
        shape = design.geometry.regions[index].primitive
        terminal = design.terminal_map[index] - first_index + 1
        if shape isa DataModel.Disk
            dx, dy = shape.at.x - x, shape.at.y - y
            push!(wires,
                (x = T(cos(phi)*dx + sin(phi)*dy),
                    y = T(-sin(phi)*dx + cos(phi)*dy), r = T(shape.r), terminal))
        else
            zero(shape.span) < shape.span < 2pi ||
                return reject("a tape spans no angle or a full turn")
            push!(tapes,
                (ri = T(shape.ri), ro = T(shape.ro),
                    phi = T(shape.at.φ - phi), span = T(shape.span), terminal))
        end
    end
    # Touching wires of one terminal use exposed surface points in collocation.
    # Buried tape faces still require an exposed-union adapter.
    validate(wires, DielectricGroup, design, scale)
    for tape in tapes
        # A tape face that touches a metal surface belongs to that conductor.
        isempty(left) && abs(_shunt_scalar(tape.ri - host.ri)) <= _shunt_tol(scale) &&
            return reject("a tape touches the inner anchor")
        isempty(right) && abs(_shunt_scalar(tape.ro - host.ro)) <= _shunt_tol(scale) &&
            return reject("a tape touches the reference")
    end
    for tape in tapes, wire in wires
        gap = tape.ri - hypot(wire.x, wire.y) - wire.r
        _shunt_scalar(gap) >= -_shunt_tol(scale) || return reject("a wire overlaps a tape")
        abs(_shunt_scalar(gap)) <= _shunt_tol(scale) &&
            !_shunt_same(gap, zero(gap), scale) &&
            return reject("an uncertain displacement changes whether a wire touches a tape")
        if abs(_shunt_scalar(gap)) <= _shunt_tol(scale) && wire.terminal != tape.terminal
            theta = atan(sin(atan(wire.y, wire.x)-tape.phi),
                cos(atan(wire.y, wire.x)-tape.phi))
            abs(_shunt_scalar(theta)) <= _shunt_scalar(tape.span)/2 &&
                return reject("a wire touches a tape of another terminal")
        end
    end
    for i in eachindex(tapes), j in 1:(i - 1)
        l, r = tapes[i], tapes[j]
        separation = abs(atan(sin(l.phi-r.phi), cos(l.phi-r.phi)))
        radial = _shunt_scalar(max(l.ri, r.ri)-min(l.ro, r.ro))
        angular = _shunt_scalar(separation - (l.span+r.span)/2)
        radial > _shunt_tol(scale) || angular > _shunt_tol(scale)/_shunt_scalar(scale) ||
            return reject("two tapes overlap")
    end
    return push!(groups, DielectricGroup{T}(design_index, assembly,
        first_index:last_index, reference, T(host.ri), T(host.ro),
        _shunt_permittivity(convert(Material{T}, host_region.source.material), formulas),
        left, right, wires, tapes))
end

function validate(wires::Vector{ShuntWire{T}}, ::Type{DielectricGroup}, design, scale) where {T}
    for i in eachindex(wires), j in 1:(i - 1)
        wire, other = wires[i], wires[j]
        gap = hypot(wire.x-other.x, wire.y-other.y)-wire.r-other.r
        _shunt_scalar(gap) > _shunt_tol(scale) ||
            wire.terminal == other.terminal || throw(BoundarySolveError(:unsupported,
                (; cable = design.cable_id), "boundary shunt: touching wire terminals"))
        # A contact that changes with uncertain inputs has no fixed exposed topology.
        abs(_shunt_scalar(gap)) > _shunt_tol(scale) || _shunt_same(gap, zero(gap), scale) ||
            throw(BoundarySolveError(:unsupported, (; cable = design.cable_id),
                "boundary shunt cannot differentiate an uncertain wire contact that changes topology"))
    end
    return wires
end

# Extract the annular hosts of one design: each one a group or a rejection.
function _shunt_extract(design::CableDesign, blueprint::CableBlueprint, ::Type{T},
        design_index, formulas) where {T}
    groups, rejected = DielectricGroup{T}[], ShuntRejection[]
    for placed in design.geometry.regions
        shape = placed.primitive
        placed.terminal === nothing && shape isa DataModel.DifferenceShape &&
        shape.outer isa DataModel.Annulus || continue
        _shunt_host_group!(groups, rejected, design, blueprint, placed, design_index, T,
            formulas)
    end
    return (; groups, rejected)
end

"""
$(TYPEDSIGNATURES)

Return the physical descriptors of `group` in one flat vector: host radii and permittivity,
each layer's radii and permittivity, each wire's center and radius, and each tape's radii,
angle and span. It is the differentiation input: a dense factorization cannot hold uncertain
scalars, so Measurements differentiates through these values. `DielectricGroup(group,
values)` is its inverse.
"""
function Base.Vector(group::DielectricGroup{T}) where {T}
    values = T[group.a, group.b, group.epsilon]
    for layers in (group.left, group.right), layer in layers

        append!(values, (layer.ri, layer.ro, layer.epsilon))
    end
    for wire in group.wires
        append!(values, (wire.x, wire.y, wire.r))
    end
    for tape in group.tapes
        append!(values, (tape.ri, tape.ro, tape.phi, tape.span))
    end
    return values
end

"""
$(TYPEDSIGNATURES)

Rebuild `group` from its physical descriptors `values`, in the order of `Vector(group)`, with
their element type. The reference keeps its nominal radius: the calculation does not read it.
"""
function DielectricGroup(group::DielectricGroup, values::AbstractVector{S}) where {S}
    left, right = ShuntLayer{S}[], ShuntLayer{S}[]
    cursor = 4
    for (destination, layers) in ((left, group.left), (right, group.right)), _ in layers

        push!(destination,
            (ri = values[cursor], ro = values[cursor + 1], epsilon = values[cursor + 2]))
        cursor += 3
    end
    wires, tapes = ShuntWire{S}[], ShuntTape{S}[]
    for wire in group.wires
        push!(wires, (x = values[cursor], y = values[cursor + 1], r = values[cursor + 2],
            terminal = wire.terminal))
        cursor += 3
    end
    for tape in group.tapes
        push!(tapes, (ri = values[cursor], ro = values[cursor + 1], phi = values[cursor + 2],
            span = values[cursor + 3], terminal = tape.terminal))
        cursor += 4
    end
    reference = convert(Union{ShuntShield{S}, ShuntJacket{S}}, nominal(group.reference))
    return DielectricGroup{S}(group.design, group.assembly, group.terminals, reference,
        values[1], values[2], values[3], left, right, wires, tapes)
end

_shunt_lossless(::Any) = false
function _shunt_lossless(::Union{InsulationAdmittance.Formula{:lossless},
        SemiconAdmittance.Formula{:lossless}})
    true
end
function _shunt_lossless(formulas::NamedTuple)
    _shunt_lossless(formulas.insulation_admittance) &&
        _shunt_lossless(formulas.semicon_admittance)
end
