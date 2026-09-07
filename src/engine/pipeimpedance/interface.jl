"""
$(TYPEDEF)

Store a pipe-type impedance selection until the backend checks the topology.
This value provides no substitute numerical formula for an unsupported pipe.

$(TYPEDFIELDS)
"""
struct Formula{ID, A <: NamedTuple} <: PipeImpedanceFormulation
    "Declared assumptions owned by the selected pipe formula."
    assumptions::A
end

formula_id(::Formula{ID}) where {ID} = ID
assumptions(selected::Formula) = selected.assumptions

"""
$(TYPEDSIGNATURES)

Construct a registered pipe-type selection. Unknown identifiers or assumptions
raise `ArgumentError`; backend applicability is checked against the design.
"""
Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)

function Formula(::Val{ID}; kwargs...) where {ID}
    ID in FORMULAS || throw(ArgumentError("unknown pipe-impedance formula :$ID"))
    defaults = assumptions(Val(ID))
    overrides = (; kwargs...)
    isempty(setdiff(keys(overrides), keys(defaults))) || throw(ArgumentError(
        "unknown assumptions for pipe-impedance formula :$ID"))
    values = merge(defaults, overrides)
    return Formula{ID, typeof(values)}(values)
end

function Formulation(backend, selected::Formula, design::CableDesign)
    # Compare conductor axes, not wire positions. A concentric sheath remains
    # ordinary coaxial geometry even when declared through pipe(...).
    topology = Val(:coaxial)
    for wall in design.geometry.regions
        wall.source.material.kind === :conductor || continue
        annular_wall = wall.primitive isa DataModel.Annulus && wall.primitive.ri > 0
        enclosed_wall = any(entry -> entry.pattern isa DataModel.EnclosureBoundary,
            wall.placement.patterns)
        annular_wall || enclosed_wall || continue
        if !(wall.primitive isa DataModel.Annulus)
            topology = Val(:pipe)
            break
        end
        axis = DataModel.radial_position(wall)
        for terminal in design.terminal_order
            terminal === wall.terminal && continue
            sources = filter(region -> region.terminal === terminal, design.geometry.regions)
            centre = DataModel.conductor_zone_position(sources)
            distance = hypot(centre[1] - axis[1], centre[2] - axis[2])
            if distance < wall.primitive.ri && !DataModel.same_radial_position(centre, axis)
                topology = Val(:pipe)
                break
            end
        end
        topology === Val(:pipe) && break
    end
    return Formulation(backend, Val(formula_id(selected)), selected, topology)
end

function Formulation(backend, ::Val{ID}, ::Formula{ID}, ::Val{Topology}) where {ID, Topology}
    throw(ArgumentError(
        "pipe-impedance :$ID is not yet implemented for $Topology topology on $(nameof(typeof(backend)))"))
end
