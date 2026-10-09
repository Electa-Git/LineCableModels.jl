"""
$(TYPEDEF)

Store a pipe-type impedance selection until the backend checks the topology.
Unsupported pipes require a separate numerical formula. `Formula{ID}()` constructs a
registered pipe-type selection. Unknown identifiers or controls raise `ArgumentError`.
Backend applicability is checked against the design.
"""
struct Formula{ID} <: PipeImpedanceFormulation
    function Formula{ID}(; parameters::NamedTuple = (;),
            options::Union{NamedTuple, FormulationOptions} = FormulationOptions()) where {ID}
        ID in formulas(Formula) || throw(ArgumentError("unknown pipe-impedance formula :$ID"))
        options = options isa NamedTuple ? FormulationOptions(options) : options
        isempty(parameters) && isempty(options.data) ||
            throw(ArgumentError("pipe-impedance :$ID accepts no model parameters or numerical controls"))
        return new{ID}()
    end
end

"""
$(TYPEDSIGNATURES)

Check that `backend` can compute `design` with the pipe-impedance formula `selected`.
The check compares the conductor axes with each conductive enclosure: a design whose
enclosure is eccentric or contains several cores has pipe topology. Because the built-in
formulas do not add a pipe term, they admit coaxial topology only. A backend or formula
that implements pipe topology adds its own method. Return `design`.

# Errors

- Throws `ArgumentError` for a pipe topology.
"""
function validate(design::CableDesign, selected::Formula, backend)
    pipe() = throw(ArgumentError("Pipe-type cable formulation is not yet implemented for " *
        "the $(description(backend)) backend. No pipe formulation is available."))
    # Compare conductor axes, not wire positions. A concentric sheath remains
    # ordinary coaxial geometry even when declared through pipe(...).
    for wall in design.geometry.regions
        wall.source.material.kind === :conductor || continue
        annular_wall = wall.primitive isa DataModel.Annulus && wall.primitive.ri > 0
        enclosed_wall = any(entry -> entry.pattern isa DataModel.EnclosureBoundary,
            wall.placement.patterns)
        annular_wall || enclosed_wall || continue
        wall.primitive isa DataModel.Annulus || pipe()
        axis = DataModel.radial_position(wall)
        for terminal in design.terminal_order
            terminal === wall.terminal && continue
            sources = filter(region -> region.terminal === terminal, design.geometry.regions)
            center = DataModel.conductor_zone_position(sources)
            distance = hypot(center[1] - axis[1], center[2] - axis[2])
            distance < wall.primitive.ri && !DataModel.same_radial_position(center, axis) &&
                pipe()
        end
    end
    return design
end

Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)
Formula(::Val{ID}; kwargs...) where {ID} = Formula{ID}(; kwargs...)
Formula(selected::PipeImpedanceFormulation) = selected

function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    isempty(selection.options.data) ||
        throw(ArgumentError("deferred pipe contribution has no numerical options"))
    selection.equivalent_earth === nothing ||
        throw(ArgumentError("pipe contribution cannot consume equivalent_earth"))
    return Formula{ID}(; parameters = selection.parameters)
end

formula_id(::Formula{ID}) where {ID} = ID
formula_id(::Type{<:Formula{ID}}) where {ID} = ID
# Identity-only dispatch also describes retained selections without constructors.
description(value::Formula; compact::Bool = false) = description(typeof(value); compact)
formulation_options(::Formula) = FormulationOptions()

"""Expose the selected pipe equation as a native record."""
Base.NamedTuple(value::Formula) = (identifier=formula_id(value),parameters=(;),options=(;))

"""Iterate the independently selectable child slots admitted by this formula family."""
Base.pairs(::Type{<:Formula}; quantity=nothing) = pairs((;))
