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
        temperature = 20,
        frequency = 50,
        combine::Symbol = :product
)
    sources = (
        designs,
        temperature isa Union{AbstractGrid, Gridspace} ? temperature : Grid((temperature,)),
        frequency isa Union{AbstractGrid, Gridspace} ? frequency : Grid((frequency,))
    )
    return Gridspace{Engine.CableConstantsProblem}(
        _cable_constants_problem,
        sources;
        combine
    )
end

function _cable_constants_problem(design, temperature, frequency)
    return Engine.CableConstantsProblem(design; temperature, frequency)
end
