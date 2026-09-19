@testmodule UnifiedFormulaFixtures begin
    using LineCableModels
    const E=LineCableModels.Engine

    # Tests use the production allocation contract, not a second formula allocator.
    function buffers(geometry)
        T=eltype(geometry.radius)
        R=typeof(float(LineCableModels.nominal(one(T))))
        base=(quadrature = E.integration_workspace(R, Complex{T}; size = 0),
            observations = nothing)
        return E.initialize_buffers(
            E.EarthImpedance.Formula(:unified), T, (;), (; geometry), base)
    end
end
