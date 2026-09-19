# Numerical boundary-domain records owned by the local shunt formulations.
const ShuntWire{T} = NamedTuple{(:x, :y, :r, :terminal), Tuple{T, T, T, Int}}
const ShuntTape{T} = NamedTuple{(:ri, :ro, :phi, :span, :terminal), Tuple{T, T, T, T, Int}}
const ShuntLayer{T} = NamedTuple{(:ri, :ro, :material), Tuple{T, T, Material{T}}}

"""
$(TYPEDEF)

Describe an open-conductor domain inside a circular layered shield. The first
terminal is the equivalent inner conductor; the last is the closed reference
shield. Intermediate terminals retain their exposed circular wires and tapes.
Coordinates and radii are in \\[m\\]. Materials are not homogenized.

$(TYPEDFIELDS)
"""
struct InternalShuntDomain{T <: Real}
    "Source design index."
    design::Int
    "Global conductor range of the containing concentric assembly."
    assembly::UnitRange{Int}
    "Global terminal range, including inner conductor and reference shield."
    terminals::UnitRange{Int}
    "Host dielectric inner radius \\[m\\]."
    a::T
    "Host dielectric outer radius \\[m\\]."
    b::T
    "Unmodified physical host dielectric."
    material::Material{T}
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
_shunt_tol(x) = 128eps(typeof(_shunt_scalar(x))) * abs(_shunt_scalar(x))
_shunt_same(a, b, scale) = numerical_magnitude(a - b) <= _shunt_tol(scale)
function _shunt_concentric(shape, x, y, scale)
    _shunt_same(shape.at.x, x, scale) && _shunt_same(shape.at.y, y, scale)
end

# Preserve ordering and dependency information. Nominal values decide ordering;
# coincident interfaces/centres must also have matching uncertain dependencies.
function _shunt_layers(regions, inner, outer, x, y, ::Type{T}) where {T}
    layers = ShuntLayer{T}[]
    _shunt_same(inner, outer, outer) && return layers
    for placed in regions
        shape = placed.primitive
        placed.terminal === nothing && shape isa DataModel.Annulus || continue
        _shunt_concentric(shape, x, y, outer) || continue
        _shunt_scalar(shape.ri - inner) >= -_shunt_tol(outer) &&
        _shunt_scalar(shape.ro - outer) <= _shunt_tol(outer) || continue
        push!(layers,
            (ri = T(shape.ri), ro = T(shape.ro),
                material = convert(Material{T}, placed.source.material)))
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

function _shunt_host_domain(
        design, blueprint, host_region, design_index, offset, ::Type{T}) where {T}
    host = host_region.primitive.outer
    x, y, phi = host.at.x, host.at.y, host.at.φ
    scale = host.ro
    members = Int[]
    for (i, placed) in pairs(design.geometry.regions)
        placed.terminal === nothing && continue
        shape = placed.primitive
        inside = if shape isa DataModel.Disk
            radius = hypot(shape.at.x - x, shape.at.y - y)
            _shunt_scalar(radius - shape.r - host.ri) >= -_shunt_tol(scale) &&
                _shunt_scalar(radius + shape.r - host.ro) <= _shunt_tol(scale)
        elseif shape isa DataModel.BentStrip
            _shunt_concentric(shape, x, y, scale) &&
                _shunt_scalar(shape.ri - host.ri) >= -_shunt_tol(scale) &&
                _shunt_scalar(shape.ro - host.ro) <= _shunt_tol(scale)
        else
            false
        end
        if inside
            inner, outer = if shape isa DataModel.Disk
                radius = hypot(shape.at.x-x, shape.at.y-y)
                radius-shape.r, radius+shape.r
            else
                shape.ri, shape.ro
            end
            # A nominal contact with an independent uncertain displacement
            # crosses a dielectric interface, outside a fixed-topology tangent.
            for (edge, boundary) in ((inner, host.ri), (outer, host.ro))
                abs(_shunt_scalar(edge-boundary)) <= _shunt_tol(scale) &&
                    !_shunt_same(edge, boundary, scale) && return nothing
            end
            push!(members, i)
        end
    end
    isempty(members) && return nothing
    # A coating, void, or another dielectric island inside the host invalidates
    # the homogeneous-host Green function even if the metal itself fits radially.
    holes = host_region.primitive.holes
    length(holes) == length(members) && all(
        index -> any(hole -> isequal(design.geometry.regions[index].primitive, hole), holes), members) ||
        return nothing
    terminals = sort!(unique(design.terminal_map[members]))
    first(terminals) > 1 && last(terminals) < length(blueprint.conductors) || return nothing
    terminals == collect(first(terminals):last(terminals)) || return nothing
    # Never resolve only some pieces of a retained terminal.
    all(i -> design.terminal_map[i] ∉ terminals || i in members,
        eachindex(design.terminal_map)) || return nothing
    first_index, reference = first(terminals) - 1, last(terminals) + 1
    assembly_index = blueprint.conductors[first_index].assembly
    assembly = blueprint.assembly_ranges[assembly_index]
    reference in assembly || return nothing
    inner = blueprint.conductors[first_index]
    _shunt_same(inner.position[1], x, scale) &&
    _shunt_same(inner.position[2], y, scale) || return nothing
    # Inner cores keep the existing equivalent-core assumption. Elsewhere an
    # actual closed conductor is needed to separate independently solved gaps.
    first_index == first(assembly) ||
        _shunt_closed(design, first_index, x, y, scale) !== nothing || return nothing
    shield = _shunt_closed(design, reference, x, y, scale)
    shield === nothing && return nothing
    left = _shunt_layers(design.geometry.regions, inner.r_ex, host.ri, x, y, T)
    right = _shunt_layers(design.geometry.regions, host.ro, shield.ri, x, y, T)
    (left === nothing || right === nothing) && return nothing
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
            zero(shape.span) < shape.span < 2pi || return nothing
            push!(tapes,
                (ri = T(shape.ri), ro = T(shape.ro),
                    phi = T(shape.at.φ - phi), span = T(shape.span), terminal))
        end
    end
    # The whole-face element admits exposed faces only. Overlapping/touching
    # faces require an exposed-union adapter; never impose buried boundaries.
    for i in eachindex(wires), j in 1:(i - 1)

        l, r = wires[i], wires[j]
        gap = hypot(l.x-r.x, l.y-r.y)-l.r-r.r
        if _shunt_scalar(gap) <= _shunt_tol(scale)
            l.terminal == r.terminal || throw(ArgumentError(
                "internal shunt: touching wire terminals in $(design.cable_id)"))
            return nothing
        end
    end
    for tape in tapes, wire in wires

        gap = tape.ri - hypot(wire.x, wire.y) - wire.r
        _shunt_scalar(gap) >= -_shunt_tol(scale) || return nothing
        abs(_shunt_scalar(gap)) <= _shunt_tol(scale) &&
            !_shunt_same(gap, zero(gap), scale) && return nothing
        if abs(_shunt_scalar(gap)) <= _shunt_tol(scale) && wire.terminal != tape.terminal
            theta = atan(sin(atan(wire.y, wire.x)-tape.phi),
                cos(atan(wire.y, wire.x)-tape.phi))
            abs(_shunt_scalar(theta)) <= _shunt_scalar(tape.span)/2 &&
                throw(ArgumentError("internal shunt: touching wire and tape terminals in $(design.cable_id)"))
        end
    end
    for i in eachindex(tapes), j in 1:(i - 1)

        l, r = tapes[i], tapes[j]
        separation = abs(atan(sin(l.phi-r.phi), cos(l.phi-r.phi)))
        radial = _shunt_scalar(max(l.ri, r.ri)-min(l.ro, r.ro))
        angular = _shunt_scalar(separation - (l.span+r.span)/2)
        radial > _shunt_tol(scale) || angular > _shunt_tol(scale)/_shunt_scalar(scale) ||
            return nothing
    end
    return InternalShuntDomain{T}(design_index,
        (first(assembly) + offset):(last(assembly) + offset),
        (first_index + offset):(reference + offset), T(host.ri), T(host.ro),
        convert(Material{T}, host_region.source.material), left, right, wires, tapes)
end

"""
$(TYPEDSIGNATURES)

Extract eligible local shunt domains from completed physical designs. Concentric
dielectric continuity, closed reference shields and exposed faces are checked
before any numerical allocation. Unsupported geometries keep their existing
equivalent-coaxial treatment. No shapes, materials or terminals are changed.

# Returns

- Concrete vector of internal domains, with global terminal indices.
"""
function internal_shunt_domains(designs, blueprints::AbstractVector{<:CableBlueprint{T}}) where {T}
    domains = InternalShuntDomain{T}[]
    offset = 0
    for (design_index, (design, blueprint)) in enumerate(zip(designs, blueprints))
        append!(domains, internal_shunt_domains(design, blueprint, T; design_index, offset))
        offset += length(blueprint)
    end
    return domains
end

function internal_shunt_domains(
        design::CableDesign, geometry, ::Type{T}; design_index = 1, offset = 0) where {T}
    domains = InternalShuntDomain{T}[]
    for placed in design.geometry.regions
        shape = placed.primitive
        placed.terminal === nothing && shape isa DataModel.DifferenceShape &&
        shape.outer isa DataModel.Annulus || continue
        domain = _shunt_host_domain(design, geometry, placed, design_index, offset, T)
        domain === nothing && continue
        any(
            old -> max(first(old.terminals), first(domain.terminals)) <
                   min(last(old.terminals), last(domain.terminals)),
            domains) && continue
        push!(domains, domain)
    end
    return domains
end

function _shunt_domain_equal(a, b)
    length(a.terminals) == length(b.terminals) || return false
    for name in (:a, :b)
        same_physical_state(getproperty(a, name), getproperty(b, name)) || return false
    end
    for name in (:left, :right, :wires, :tapes)
        l, r = getproperty(a, name), getproperty(b, name)
        length(l) == length(r) || return false
        if name in (:left, :right)
            all(same_physical_state((x.ri, x.ro, x.material.eps_r), (
                    y.ri, y.ro, y.material.eps_r))
            for (x, y) in zip(l, r)) || return false
        else
            all(same_physical_state(x, y) for (x, y) in zip(l, r)) || return false
        end
    end
    return same_physical_state(a.material.eps_r, b.material.eps_r)
end

_shunt_lossless(::Any) = false
function _shunt_lossless(::Union{InsulationAdmittance.Formula{:lossless},
        SemiconAdmittance.Formula{:lossless}})
    true
end
function _shunt_lossless(methods::NamedTuple)
    _shunt_lossless(methods.insulation_admittance) &&
        _shunt_lossless(methods.semicon_admittance)
end
