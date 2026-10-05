@testmodule UnifiedFormulaFixtures begin
    using LineCableModels
    const E=LineCableModels.Engine
    const G=LineCableModels.Commons

    # Numerical controls use the production allocator. There is no second array layout.
    function buffers(geometry)
        T=eltype(geometry.radius)
        base=G.initialize_buffers(E.earth!, T, (; n_cables = length(geometry.radius)),
            (; geometry), (; observations = nothing))
        return G.initialize_buffers(
            E.EarthImpedance.Formula(:unified), T, (;), (; geometry), base)
    end

    # Build a real public problem and workspace. Independent equal-medium controls
    # replace completed material values after material evaluation.
    # Their artificial air conductivity is confined to this test fixture.
    function workspace(geometry,
            state,
            integration = E.formulation_options(E.SpectralIntegral, (
                method = :quad, options = (rtol = 1e-9,))))
        T=eltype(geometry.radius)
        metal=Material(:conductor, T(1.7e-8))
        origin=Pose2(zero(T), zero(T), zero(T))
        designs=[build(CableDesign,
                     "radial-$i",
                     Group(:core, origin,
                         Region(:core, Disk(radius), metal), nothing, nothing, nothing))
                 for (i, radius) in enumerate(geometry.radius)]
        system=build(LineCableSystem, designs,
            [Pose2(x, y) for (x, y) in zip(geometry.horizontal, geometry.height)];
            connections = [Dict(:core=>i) for i in eachindex(designs)], line_length = one(T))
        earth=convert(EarthModel{T}, homogeneous(rho = T(100)))
        problem=LineParametersProblem(system; temperature = T(20),
            earth_props = earth, frequencies = T[imag(state.jω) / (2pi)])
        selected=Formulation(
            earth_impedance = formula(:unified;
                options = (Γ = state.Γ,
                    integration = (method = :quad, options = integration.options))),
            earth_admittance = formula(:unified;
                options = (Γ = state.Γ,
                    integration = (method = :quad, options = integration.options))),
            options = (
                reduce_bundle = false, kron_reduction = false, ideal_transposition = false))
        execution=E.computation_options(LineCableModelsCoaxial, ComputationOptions())
        blueprints=only(E.flatten(LineCableModelsCoaxial(), problem.system.designs, T, [selected]))
        work=E.LineParametersWorkspace(problem, selected, execution, blueprints)
        E.materials!(work, selected)
        E.materials!(work, selected, 1)
        work.input.jω[1]=state.jω
        for materials in work.buffers.earth_materials
            for column in axes(materials.rho, 2), layer in 1:2

                materials.rho[layer, column]=iszero(state.sigma[layer]) ? T(Inf) :
                                             inv(state.sigma[layer])
                materials.epsilon[layer, column]=state.epsilon[layer]
                materials.mu[layer, column]=state.mu[layer]
            end
        end
        return work
    end

    function calculate(geometry, state, integration)
        work=workspace(geometry, state, integration)
        E.earth!(work, 1)
        return work.buffers
    end

    function field_state(geometry, state)
        work=workspace(geometry, state)
        binding=only(work.plan.earth_calculations)
        materials=only(work.buffers.earth_materials)
        calculation=binding.selection(materials, binding, work, 1)
        return merge(calculation.state, (radial_current = work.buffers.radial_current,))
    end
end
