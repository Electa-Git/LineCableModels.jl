"""
$(TYPEDSIGNATURES)

Map a finite space of completed cable designs to cable constants.

Each deterministic configuration or stochastic realization calls the scalar
`CableConstants(design; kwargs...)` surface exactly once.
"""
function DataModel.CableConstants(
        designs::Gridspace{<:DataModel.CableDesign};
        kwargs...
)
    caller = design -> DataModel.CableConstants(design; kwargs...)
    return Gridspace{DataModel.CableConstants}(caller, (designs,))
end
