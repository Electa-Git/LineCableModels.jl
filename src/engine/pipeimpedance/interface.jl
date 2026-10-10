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

Check that `backend` can compute `design` with a built-in pipe-impedance formula. The
topology rule of `DataModel.flatten` rejects pipe-type cables, and the built-in formulas
do not add a pipe term, so they admit every design that flattens. A pipe formulation adds
its own method. Return `design`.
"""
validate(design::CableDesign, ::Formula, backend) = design

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
