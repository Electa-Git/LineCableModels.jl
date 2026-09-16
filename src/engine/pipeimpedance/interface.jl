"""
$(TYPEDEF)

Store a pipe-type impedance selection until the backend checks the topology.
This value provides no substitute numerical formula for an unsupported pipe.

"""
struct Formula{ID} <: PipeImpedanceFormulation end

formula_id(::Formula{ID}) where {ID} = ID

"""
$(TYPEDSIGNATURES)

Construct a registered pipe-type selection. Unknown identifiers or controls
raise `ArgumentError`; backend applicability is checked against the design.
"""
Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)

function Formula(::Val{:none}; parameters::NamedTuple = (;),
        options::Union{NamedTuple, FormulationOptions} = FormulationOptions())
    options = options isa NamedTuple ? FormulationOptions(options) : options
    isempty(parameters) && isempty(options.data) ||
        throw(ArgumentError("pipe-impedance :none accepts no model parameters or numerical controls"))
    return Formula{:none}()
end

Formula(::Val{ID}; kwargs...) where {ID} = throw(ArgumentError("unknown pipe-impedance formula :$ID"))

Formula(selected::PipeImpedanceFormulation) = selected

function Formulation(backend, selected::PipeImpedanceFormulation, design::CableDesign)
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
    return Formulation(backend, selected, topology)
end

function Formulation(backend, selected::PipeImpedanceFormulation, ::Val{Topology}) where {Topology}
    throw(ArgumentError(
        "pipe-impedance :$(formula_id(selected)) is not yet implemented for $Topology topology on $(nameof(typeof(backend)))"))
end

function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    isempty(selection.options.data) ||
        throw(ArgumentError("deferred pipe contribution has no numerical options"))
    selection.equivalent_earth === nothing ||
        throw(ArgumentError("pipe contribution cannot consume equivalent_earth"))
    return Formula(Val(ID); parameters = selection.parameters)
end

"""Expose the selected pipe equation as a native record."""
Base.NamedTuple(value::Formula) = (identifier=formula_id(value),parameters=(;),options=(;))

# Identity-only dispatch also describes retained selections without constructors.
import ...Grammar: formulation_options
description(value::Formula; compact::Bool=false) = description(typeof(value); compact)

"""Iterate the independently selectable child slots admitted by this formula family."""
Base.pairs(::Type{<:Formula}; quantity=nothing) = pairs((;))
formula_id(::Type{<:Formula{ID}}) where {ID} = ID
formulation_options(::Formula) = FormulationOptions()
