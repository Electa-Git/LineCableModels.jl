# Dielectric group records owned by the local shunt formulations. A wire or tape records the
# index of the dielectric layer it lies in. A wire records its metal radius `r`, the outer
# radius `coat` of its coat and the coat's relative permittivity `epsilon`. A bare wire's coat
# radius is its metal radius, and its coat permittivity is that of its layer.
const ShuntWire{T} = NamedTuple{(:x, :y, :r, :coat, :epsilon, :terminal, :layer),
    Tuple{T, T, T, T, T, Int, Int}}
const ShuntTape{T} = NamedTuple{(:ri, :ro, :phi, :span, :terminal, :layer),
    Tuple{T, T, T, T, Int, Int}}
const ShuntLayer{T} = NamedTuple{(:ri, :ro, :epsilon), Tuple{T, T, T}}
# Fills that the formula cannot represent, with their terminals and the reason.
const ShuntRejection = NamedTuple{(:terminals, :reason), Tuple{Vector{Int}, BoundarySolveError}}

"""
$(TYPEDEF)

Hold the dielectric layers between two closed equipotentials, the inner anchor and the
shield or jacket reference, with the open conductors that lie in them, solved as one group.
Each open conductor lies in one layer. The first terminal is the inner anchor. A shield
occupies the last real terminal, while a jacket reference has a radius and does not own a
terminal. Coordinates and radii use \\[m\\]. Relative permittivities are those of the
formulation's lossless laws. The extraction builds each group once, and the boundary shunt
calculation consumes it.

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
    "Concentric dielectric layers from the inner anchor to the reference, in radial order."
    layers::Vector{ShuntLayer{T}}
    "Exposed circular wires in the group frame, each with its coat and the index of its layer."
    wires::Vector{ShuntWire{T}}
    "Exposed finite-thickness tapes in the group frame, each with the index of its layer."
    tapes::Vector{ShuntTape{T}}
end

_shunt_scalar(x) = float(nominal(x))
# Faces touch when their gap is below this tolerance, which scales with the outer radius of
# their layer.
_shunt_tol(x) = 128eps(typeof(_shunt_scalar(x))) * abs(_shunt_scalar(x))
_shunt_same(a, b, scale) = tolerance_norm(a - b) <= _shunt_tol(scale)
function _shunt_concentric(shape, x, y, scale)
    _shunt_same(shape.at.x, x, scale) && _shunt_same(shape.at.y, y, scale)
end

# The relative permittivity of `source` under the formulation's lossless laws.
function _shunt_permittivity(source, ::Type{T}, formulas) where {T}
    material = convert(Material{T}, source)
    frequency = one(T)
    kappa = constitutive((formulas.insulation_admittance, formulas.semicon_admittance),
        material, frequency, material.T0; temperature_dependence = nothing)
    return imag(kappa)/(2pi*frequency*vacuum_permittivity(T))
end

# The concentric dielectric layers between `inner` and `outer`, the fills of open conductors
# `fills` among them, in radial order, each with the relative permittivity of its material
# under the formulation's lossless laws. Preserve ordering and dependency information. Nominal
# values decide ordering. Coincident interfaces and centers must also have matching uncertain
# dependencies.
function _shunt_layers(regions, fills, inner, outer, x, y, ::Type{T}, formulas) where {T}
    inner, outer = T(inner), T(outer)
    layers = ShuntLayer{T}[]
    _shunt_same(inner, outer, outer) && return layers
    annuli = Tuple{DataModel.Annulus, Material}[(placed.primitive, placed.source.material)
                                                for placed in regions
              if placed.terminal === nothing && placed.primitive isa DataModel.Annulus &&
              _shunt_concentric(placed.primitive, x, y, outer) &&
              _shunt_scalar(placed.primitive.ri - inner) >= -_shunt_tol(outer) &&
              _shunt_scalar(placed.primitive.ro - outer) <= _shunt_tol(outer)]
    append!(annuli, [(regions[fill].primitive.outer, regions[fill].source.material)
                     for fill in fills])
    for (shape, source) in annuli
        push!(layers, (ri = T(shape.ri), ro = T(shape.ro),
            epsilon = _shunt_permittivity(source, T, formulas)))
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

# Add the group of the fills `chain` to `groups`, or the reason that the boundary shunt formula
# cannot represent it to `rejected`. The fills lie between one inner anchor and one reference,
# with no closed conductor between them. Each record contains the region index of its fill,
# the region indices of the open conductors in it and whether an uncertain displacement
# changes a contact with its inner or outer surface.
function _shunt_group!(groups, rejected, design, blueprint, chain, design_index,
        ::Type{T}, formulas) where {T}
    regions = design.geometry.regions
    innermost = regions[first(chain).region].primitive.outer
    x, y, phi = innermost.at.x, innermost.at.y, innermost.at.φ
    scale = T(regions[last(chain).region].primitive.outer.ro)
    members = mapreduce(current -> current.members, vcat, chain)
    terminals = sort!(unique(design.terminal_map[members]))
    reject(reason) = push!(rejected, (; terminals, reason = BoundarySolveError(:unsupported,
        (; cable = design.cable_id, terminals),
        "boundary shunt cannot represent this dielectric group: " * reason)))
    any(current -> current.uncertain, chain) && return reject(
        "an uncertain displacement changes whether a wire or tape touches a layer boundary")
    for current in chain
        holes = regions[current.region].primitive.holes
        for hole in holes
            any(member -> isequal(DataModel.fill_hole(regions, member), hole),
                current.members) && continue
            # A conductor that leaves its fill crosses a dielectric interface.
            any(placed -> placed.terminal !== nothing && isequal(placed.primitive, hole),
                regions) && return reject("an open conductor crosses a dielectric interface")
            return reject(
                "it holds coatings or dielectric islands besides its exposed wires and tapes")
        end
        length(holes) == length(current.members) && all(member -> any(
            hole -> isequal(DataModel.fill_hole(regions, member), hole), holes),
            current.members) ||
            return reject(
                "it holds coatings or dielectric islands besides its exposed wires and tapes")
    end
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
        jacket = blueprint[last(assembly)].layers
        isempty(jacket) && return reject("the outermost conductor has no jacket to reference")
        ShuntJacket{T}(last(jacket).r_ex), last(terminals)
    end
    layers = _shunt_layers(regions, [current.region for current in chain], inner.r_ex,
        reference.radius, x, y, T, formulas)
    layers === nothing && return reject(
        "the dielectrics between the inner anchor and the reference are not continuous concentric layers")
    wires, tapes = ShuntWire{T}[], ShuntTape{T}[]
    for current in chain
        boundary = regions[current.region].primitive.outer
        # The fill is one of the layers: its open conductors lie in it.
        layer = findfirst(layer -> _shunt_same(layer.ri, boundary.ri, scale), layers)::Int
        for index in current.members
            shape = regions[index].primitive
            terminal = design.terminal_map[index] - first_index + 1
            if shape isa DataModel.Disk
                dx, dy = shape.at.x - x, shape.at.y - y
                # A wire's hole is its own disk, or the outer circle of the coat that follows it.
                hole = DataModel.fill_hole(regions, index)
                coat, epsilon = isequal(hole, shape) ? (T(shape.r), layers[layer].epsilon) :
                                (T(hole.r),
                    _shunt_permittivity(regions[index + 1].source.material, T, formulas))
                push!(wires,
                    (x = T(cos(phi)*dx + sin(phi)*dy), y = T(-sin(phi)*dx + cos(phi)*dy),
                        r = T(shape.r), coat, epsilon, terminal, layer))
            else
                zero(shape.span) < shape.span < 2pi ||
                    return reject("a tape spans no angle or a full turn")
                push!(tapes,
                    (ri = T(shape.ri), ro = T(shape.ro), phi = T(shape.at.φ - phi),
                        span = T(shape.span), terminal, layer))
            end
        end
    end
    # Touching wires of one terminal use exposed surface points in collocation.
    # Buried tape faces still require an exposed-union adapter.
    validate(wires, DielectricGroup, design, layers, scale)
    for tape in tapes
        # A tape face that touches a metal surface belongs to that conductor.
        tape.layer == 1 && abs(_shunt_scalar(tape.ri - first(layers).ri)) <= _shunt_tol(scale) &&
            return reject("a tape touches the inner anchor")
        tape.layer == length(layers) &&
            abs(_shunt_scalar(tape.ro - last(layers).ro)) <= _shunt_tol(scale) &&
            return reject("a tape touches the reference")
    end
    for tape in tapes, wire in wires
        gap = tape.ri - hypot(wire.x, wire.y) - wire.coat
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
    return push!(groups, DielectricGroup{T}(design_index, assembly, first_index:last_index,
        reference, layers, wires, tapes))
end

function validate(wires::Vector{ShuntWire{T}}, ::Type{DielectricGroup}, design, layers,
        scale) where {T}
    for wire in wires
        layer, pitch = layers[wire.layer], hypot(wire.x, wire.y)
        # A coat lies in the layer of its wire. It can touch an interface but not cross it.
        _shunt_scalar(pitch - wire.coat - layer.ri) >= -_shunt_tol(scale) &&
            _shunt_scalar(layer.ro - pitch - wire.coat) >= -_shunt_tol(scale) ||
            throw(BoundarySolveError(:unsupported, (; cable = design.cable_id),
                "boundary shunt: a coat crosses a dielectric interface"))
    end
    for i in eachindex(wires), j in 1:(i - 1)
        wire, other = wires[i], wires[j]
        distance = hypot(wire.x-other.x, wire.y-other.y)
        gap = distance-wire.r-other.r
        _shunt_scalar(gap) > _shunt_tol(scale) ||
            wire.terminal == other.terminal || throw(BoundarySolveError(:unsupported,
                (; cable = design.cable_id), "boundary shunt: touching wire terminals"))
        # Coats can touch each other but not overlap. Bare wires of one terminal can overlap.
        outline = distance-wire.coat-other.coat
        wire.coat == wire.r && other.coat == other.r ||
            _shunt_scalar(outline) >= -_shunt_tol(scale) ||
            throw(BoundarySolveError(:unsupported, (; cable = design.cable_id),
                "boundary shunt: overlapping coats"))
        # A contact that changes with uncertain inputs has no fixed exposed topology.
        abs(_shunt_scalar(outline)) > _shunt_tol(scale) ||
            _shunt_same(outline, zero(outline), scale) ||
            throw(BoundarySolveError(:unsupported, (; cable = design.cable_id),
                "boundary shunt cannot differentiate an uncertain wire contact that changes topology"))
    end
    return wires
end

# Extract the dielectric groups of one design: each one a group or a rejection. The fills of
# open conductors around one center with no closed conductor between them form one group.
function _shunt_extract(design::CableDesign, blueprint::CableBlueprint, ::Type{T},
        design_index, formulas) where {T}
    groups, rejected = DielectricGroup{T}[], ShuntRejection[]
    regions = design.geometry.regions
    fills = NamedTuple{(:region, :members, :uncertain), Tuple{Int, Vector{Int}, Bool}}[]
    for (region, placed) in pairs(regions)
        shape = placed.primitive
        placed.terminal === nothing && shape isa DataModel.DifferenceShape &&
        shape.outer isa DataModel.Annulus || continue
        annulus = shape.outer
        x, y = annulus.at.x, annulus.at.y
        scale = T(annulus.ro)
        members = Int[]
        uncertain = false
        for (i, conductor) in pairs(regions)
            conductor.terminal === nothing && continue
            piece = conductor.primitive
            inner, outer = if piece isa DataModel.Disk
                radius = hypot(piece.at.x - x, piece.at.y - y)
                radius-piece.r, radius+piece.r
            elseif piece isa DataModel.BentStrip && _shunt_concentric(piece, x, y, scale)
                piece.ri, piece.ro
            else
                continue
            end
            _shunt_scalar(inner - annulus.ri) >= -_shunt_tol(scale) &&
                _shunt_scalar(outer - annulus.ro) <= _shunt_tol(scale) || continue
            # A nominal contact with an independent uncertain displacement
            # crosses a dielectric interface, outside a fixed-topology tangent.
            for (edge, boundary) in ((inner, annulus.ri), (outer, annulus.ro))
                uncertain |= abs(_shunt_scalar(edge-boundary)) <= _shunt_tol(scale) &&
                             !_shunt_same(edge, boundary, scale)
            end
            push!(members, i)
        end
        # A fill without open conductors is only a dielectric layer of a group.
        isempty(members) || push!(fills, (; region, members, uncertain))
    end
    sort!(fills; by = current -> minimum(design.terminal_map[current.members]))
    chain = similar(fills, 0)
    for current in fills
        if !isempty(chain)
            outer = regions[last(chain).region].primitive.outer
            annulus = regions[current.region].primitive.outer
            # A closed conductor between two fills takes the terminal between theirs.
            continues = minimum(design.terminal_map[current.members]) ==
                        maximum(design.terminal_map[last(chain).members]) + 1 &&
                        _shunt_concentric(annulus, outer.at.x, outer.at.y, annulus.ro)
            if !continues
                _shunt_group!(groups, rejected, design, blueprint, chain, design_index, T,
                    formulas)
                chain = similar(fills, 0)
            end
        end
        push!(chain, current)
    end
    isempty(chain) ||
        _shunt_group!(groups, rejected, design, blueprint, chain, design_index, T, formulas)
    return (; groups, rejected)
end

"""
$(TYPEDSIGNATURES)

Return the physical descriptors of `group` in one flat vector. The layers come first, in
radial order, with their radii and permittivity. The wires follow with their centers, metal
radii, coat radii and coat permittivities, then the tapes with their radii, angles and spans. It is the differentiation input: a
dense factorization cannot hold uncertain scalars, so Measurements differentiates through
these values. The layer of a wire or tape is discrete and is not a descriptor.
`DielectricGroup(group, values)` is its inverse.
"""
function Base.Vector(group::DielectricGroup{T}) where {T}
    values = T[]
    for layer in group.layers
        append!(values, (layer.ri, layer.ro, layer.epsilon))
    end
    for wire in group.wires
        append!(values, (wire.x, wire.y, wire.r, wire.coat, wire.epsilon))
    end
    for tape in group.tapes
        append!(values, (tape.ri, tape.ro, tape.phi, tape.span))
    end
    return values
end

"""
$(TYPEDSIGNATURES)

Rebuild `group` from its physical descriptors `values`, in the order of `Vector(group)`, with
their element type. The terminals and layers of the wires and tapes and the nominal radius of
the reference are those of `group`: the calculation does not read the reference radius.
"""
function DielectricGroup(group::DielectricGroup, values::AbstractVector{S}) where {S}
    layers = ShuntLayer{S}[]
    cursor = 1
    for _ in group.layers
        push!(layers,
            (ri = values[cursor], ro = values[cursor + 1], epsilon = values[cursor + 2]))
        cursor += 3
    end
    wires, tapes = ShuntWire{S}[], ShuntTape{S}[]
    for wire in group.wires
        push!(wires, (x = values[cursor], y = values[cursor + 1], r = values[cursor + 2],
            coat = values[cursor + 3], epsilon = values[cursor + 4], wire.terminal, wire.layer))
        cursor += 5
    end
    for tape in group.tapes
        push!(tapes, (ri = values[cursor], ro = values[cursor + 1], phi = values[cursor + 2],
            span = values[cursor + 3], tape.terminal, tape.layer))
        cursor += 4
    end
    reference = convert(Union{ShuntShield{S}, ShuntJacket{S}}, nominal(group.reference))
    return DielectricGroup{S}(group.design, group.assembly, group.terminals, reference,
        layers, wires, tapes)
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
