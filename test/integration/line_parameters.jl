@testitem "Engine / solver / multicable reciprocity and modal transformation" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures,
    TestNumerics
] begin
    using LinearAlgebra

    problem=TestFixtures.line_parameters_problem(frequencies = [50.0, 500.0])
    formulation=Formulation(
        Val(:EMT);
        modal_transform = Fortescue(),
        options = (
            reduce_bundle = true,
            kron_reduction = true,
            ideal_transposition = false,
            temperature_correction = true
        )
    )
    workspace, parameters=LineCableModels.Engine._compute_with_workspace(
        problem,
        formulation,
        ComputeOptions(store_primitive_matrices = true)
    )

    @test domain(parameters) === ModalDomain
    @test size(parameters.Z) == (3, 3, 2)
    @test all(isfinite, parameters.Z)
    @test all(isfinite, parameters.Y)
    @test workspace.Zg !== nothing
    @test workspace.Pg !== nothing
    for frequency_index in eachindex(problem.frequencies)
        @test TestNumerics.isapprox_scaled(
            workspace.Zg[:, :, frequency_index],
            transpose(workspace.Zg[:, :, frequency_index])
        )
        @test TestNumerics.isapprox_scaled(
            workspace.Pg[:, :, frequency_index],
            transpose(workspace.Pg[:, :, frequency_index])
        )
        @test any(!iszero,
            workspace.Zg[:, :, frequency_index] .-
            Diagonal(diag(workspace.Zg[:, :, frequency_index])))
        @test any(!iszero,
            workspace.Pg[:, :, frequency_index] .-
            Diagonal(diag(workspace.Pg[:, :, frequency_index])))
    end
end

@testitem "Engine / solver / bundle-only and singleton reduction policies" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures
] begin
    design=TestFixtures.mv_cable_design()
    duplicate_mapping=Dict("core"=>1, "sheath"=>1, "jacket"=>0)
    duplicate_position=CablePosition(design, 0.0, -1.0, duplicate_mapping)
    duplicate_system=LineCableSystem("duplicate-bundle", 500.0, duplicate_position)
    frequencies=[50.0]
    earth=EarthModel(frequencies, 100.0, 10.0, 1.0)
    duplicate_problem=LineParametersProblem(
        duplicate_system;
        earth_props = earth,
        frequencies
    )
    bundle_only=Formulation(
        Val(:EMT);
        options = (
            reduce_bundle = true,
            kron_reduction = false,
            ideal_transposition = true,
            temperature_correction = false
        )
    )
    duplicate_result=compute!(
        duplicate_problem,
        bundle_only;
        options = ComputeOptions(store_primitive_matrices = false)
    )
    @test size(duplicate_result.Z) == (2, 2, 1)
    @test all(isfinite, duplicate_result.Z)
    @test all(isfinite, duplicate_result.Y)

    singleton_design=CableDesign("single-component", deepcopy(first(design.components)))
    singleton_position=CablePosition(
        singleton_design,
        0.0,
        -1.0,
        Dict("core"=>1)
    )
    singleton_system=LineCableSystem("singleton", 100.0, singleton_position)
    singleton_problem=LineParametersProblem(
        singleton_system;
        earth_props = earth,
        frequencies
    )
    unreduced=Formulation(
        Val(:EMT);
        options = (
            reduce_bundle = false,
            kron_reduction = false,
            ideal_transposition = true,
            temperature_correction = false
        )
    )
    singleton_result=compute!(singleton_problem, unreduced)
    @test size(singleton_result.Z) == (1, 1, 1)
    @test real(singleton_result.Z[1, 1, 1]) > 0
    @test imag(singleton_result.Y[1, 1, 1]) > 0
end
