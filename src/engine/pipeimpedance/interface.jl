"""
$(TYPEDEF)

Store a pipe-type impedance selection until the backend checks the topology.
This value provides no substitute numerical formula for an unsupported pipe.

"""
struct Formula{ID} <: PipeImpedanceFormulation end

formula_id(::Formula{ID}) where {ID} = ID

"""
$(TYPEDSIGNATURES)

Construct a registered pipe-type selection. Unknown identifiers or customizations
raise `ArgumentError`; backend applicability is checked against the design.
"""
Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)

function Formula(::Val{ID}; parameters::NamedTuple = (;), hooks::NamedTuple = (;)) where {ID}
    ID in FORMULAS || throw(ArgumentError("unknown pipe-impedance formula :$ID"))
    isempty(hooks) || throw(ArgumentError("pipe impedance has no configurable hooks"))
    isempty(parameters) ||
        throw(ArgumentError("pipe impedance has no configurable parameters"))
    return Formula{ID}()
end

Formula(selected::Formula) = selected

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

function Formulation(backend, ::Val{ID}, ::Formula{ID}, ::Val{Topology}) where {
        ID, Topology}
    throw(ArgumentError(
        "pipe-impedance :$ID is not yet implemented for $Topology topology on $(nameof(typeof(backend)))"))
end

function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    isempty(selection.options) ||
        throw(ArgumentError("deferred pipe contribution has no numerical options"))
    selection.equivalent_earth === nothing ||
        throw(ArgumentError("pipe contribution cannot consume equivalent_earth"))
    return Formula(Val(ID); parameters = selection.parameters, hooks = selection.hooks)
end

"""Expose the selected pipe equation as a native record."""
Base.NamedTuple(value::Formula) = (identifier=formula_id(value),)
