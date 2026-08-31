"""
$(TYPEDSIGNATURES)

Map a finite space of completed cable designs to cable constants.

Each deterministic configuration or stochastic realization calls the scalar
`CableConstants(design; kwargs...)` surface exactly once.
"""
function Engine.CableConstants(
        designs::Gridspace{<:DataModel.CableDesign};
        kwargs...
)
    caller = design -> Engine.CableConstants(design; kwargs...)
    return Gridspace{Engine.CableConstants}(caller, (designs,))
end

"""
$(TYPEDSIGNATURES)

Map a finite space of completed cable designs to earth-free cable-constant
problems.
"""
function Engine.CableConstantsProblem(
        designs::Gridspace{<:DataModel.CableDesign};
        kwargs...
)
    caller = design -> Engine.CableConstantsProblem(design; kwargs...)
    return Gridspace{Engine.CableConstantsProblem}(caller, (designs,))
end
