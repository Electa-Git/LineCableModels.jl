@testitem "Engine / indexed reuse compares current material values and uncertainty" tags=[:unit] setup=[
    TestFixtures, EarthInteractionFixtures] begin
    using Measurements
    const E=LineCableModels.Engine
    const F=EarthInteractionFixtures
    # Equal plain floating-point infinities need no subtraction (Inf-Inf is NaN).
    @test E.same_physical_state(BigFloat(Inf), BigFloat(Inf))
    @test !E.same_physical_state(BigFloat(Inf), BigFloat(-Inf))
    problem=TestFixtures.three_bare_wires_problem(
        heights = (-1.0, -1.0, -1.0), frequencies = [50.0], rho = 100.0)
    selected=F.PairImpedance()
    work=F.workspace(problem, selected)
    binding=first(work.invariants.earth_calculations)
    @test count(iszero, binding.previous)==3
    E.earth!(work, 1)
    @test length(selected.calls)==3
    original=copy(work.buffers.Zearth)
    materials=first(work.buffers.earth_materials)
    # The final diagonal shares geometry with the first two, but its material changes.
    materials.rho[2, 9]=200.0
    empty!(selected.calls)
    E.earth!(work, 1)
    @test length(selected.calls)==4
    @test work.buffers.Zearth[3, 3]≈original[3, 3]+(1e-4+1e-3im)
    materials.rho[2, 9]=100.0
    empty!(selected.calls)
    E.earth!(work, 1)
    @test length(selected.calls)==3
    @test work.buffers.Zearth==original

    # A later pair can reuse the second material class in the same candidate chain.
    materials.rho[2, 5]=materials.rho[2, 9]=200.0
    empty!(selected.calls)
    E.earth!(work, 1)
    @test length(selected.calls)==4
    @test work.buffers.Zearth[2, 2]==work.buffers.Zearth[3, 3]
    @test work.buffers.Zearth[1, 1]==original[1, 1]

    uncertain_problem=TestFixtures.three_bare_wires_problem(
        heights = (-1.0, -1.0, -1.0), frequencies = [50.0], rho = measurement(100.0),
        radius = measurement(0.0425))
    uncertain=F.workspace(uncertain_problem)
    selected=first(uncertain.invariants.earth_calculations).selection
    materials=first(uncertain.buffers.earth_materials)
    rho=measurement(100.0, 1.0)
    independent=measurement(100.0, 1.0)
    materials.rho[2, :].=rho
    materials.rho[2, 9]=independent
    E.earth!(uncertain, 1)
    @test length(selected.calls)==4
    @test uncertainty(imag(uncertain.buffers.Zearth[1, 2]-uncertain.buffers.Zearth[2, 1]))==0
    @test uncertainty(imag(uncertain.buffers.Zearth[1, 1]-uncertain.buffers.Zearth[3, 3]))≈sqrt(2)*1e-5
    materials.rho[2, 9]=2rho-rho
    empty!(selected.calls)
    E.earth!(uncertain, 1)
    @test length(selected.calls)==3
    @test uncertainty(imag(uncertain.buffers.Zearth[1, 1]-uncertain.buffers.Zearth[3, 3]))==0

    height=measurement(-1.0, 0.001)
    independent_height=measurement(-1.0, 0.001)
    geometry_problem=TestFixtures.three_bare_wires_problem(
        heights = (height, independent_height, height), frequencies = [50.0], rho = 100.0,
        radius = measurement(0.0425))
    geometry_work=F.workspace(geometry_problem)
    geometry_binding=first(geometry_work.invariants.earth_calculations)
    @test geometry_binding.previous[5]==0
    @test geometry_binding.previous[9]==1
end

@testitem "Engine / unified shares indexed source coefficients before its complete solve" tags=[:unit] setup=[UnifiedFormulaFixtures] begin
    const E=LineCableModels.Engine
    positions=[(column*0.1905, -1.1335-row*0.19035) for row in 0:2 for column in -1:1]
    geometry=(horizontal = first.(positions),
        height = last.(positions), radius = fill(0.027358, 9))
    controls=E.formulation_options(E.SpectralIntegral, (
        method = :quad, options = (rtol = 1e-8,)))
    for (frequency, gamma) in ((50.0, 0im), (500.0, 1e-4+2e-4im))
        state=(jω = 2pi*im*frequency, Γ = gamma, sigma = (0.0, 0.01),
            epsilon = (8.8541878128e-12, 10*8.8541878128e-12), mu = (4pi*1e-7, 4pi*1e-7))
        reused=UnifiedFormulaFixtures.workspace(geometry, state, controls)
        separate=UnifiedFormulaFixtures.workspace(geometry, state, controls)
        binding=only(reused.invariants.earth_calculations)
        @test count(iszero, binding.previous)==27
        fill!(only(separate.invariants.earth_calculations).previous, 0)
        E.earth!(reused, 1)
        E.earth!(separate, 1)
        @test length(unique(reused.buffers.earth_interactions.representatives))==27
        @test length(unique(separate.buffers.earth_interactions.representatives))==81
        for name in (:axial_field, :source_potential, :current_map,
            :enclosed_impedance, :enclosed_potential)
            @test getproperty(reused.buffers, name)==getproperty(separate.buffers, name)
        end
    end
end

@testitem "Engine / ordinary earth formulations reuse through public compute" tags=[:unit] setup=[
    TestFixtures, EarthInteractionFixtures] begin
    const F=EarthInteractionFixtures
    problem=TestFixtures.three_bare_wires_problem(heights = (-1.0, -1.0, -1.0),
        frequencies = [50.0, 500.0], rho = 100.0)
    impedance, potential=F.PairImpedance(), F.PairPotential()
    selected=Formulation(earth_impedance = impedance, earth_admittance = potential;
        options = (
            reduce_bundle = false, kron_reduction = false, ideal_transposition = false))
    result=compute(problem, selected; options = (trace = true, verbosity = (default = 0,)))
    @test length(impedance.calls) == 6
    @test length(potential.calls) == 6
    trace=details(result).data.trace
    for frequency in eachindex(problem.frequencies), row in 1:3, column in 1:3
        coefficient=10(row==column)+3+abs(row-column)+4*0.0425+
                    1+problem.frequencies[frequency]/1000
        @test trace.Zg[row, column, frequency] ≈ coefficient*(1e-4+1e-3im)
        @test trace.Pg[row, column, frequency] ≈ coefficient*1e9
    end
    @test all(isfinite, Z(result)) && all(isfinite, Y(result))
end

@testitem "Engine / indexed integral reuse retains descriptions and logical positions" tags=[:unit] setup=[
    TestFixtures, EarthInteractionFixtures] begin
    using Logging
    const E=LineCableModels.Engine
    const F=EarthInteractionFixtures
    problem=TestFixtures.three_bare_wires_problem(
        heights = (-1.0, -1.0, -1.0), frequencies = [50.0], rho = 100.0)
    for description in (nothing, "earth integral"), trace in (false, true),
        positions in (false, true)
        selected=F.IntegralImpedance(description; positions)
        work=F.workspace(problem, selected; trace)
        logger=Test.TestLogger()
        with_logger(logger) do
            E.earth!(work, 1)
        end
        # Self and mutual are different indexed equations, each evaluated once.
        @test selected.calls[]==2
        @test all(==(work.buffers.Zearth[1]), work.buffers.Zearth)
        @test length(logger.logs)==9
        @test isempty(work.buffers.earth_interactions.warnings)
        if positions
            @test Set((record.kwargs[:context].receiver, record.kwargs[:context].source)
            for record in logger.logs)==Set(Tuple.(CartesianIndices(work.buffers.Zearth)))
        else
            @test all(record->record.kwargs[:context]===description, logger.logs)
        end
        trace||continue
        records=work.capture.integrals
        @test length(records)==9
        @test all(record->record.value==work.buffers.Zearth[1], records)
        if positions
            @test Set((record.context.receiver, record.context.source)
            for record in records) ==
                  Set(Tuple.(CartesianIndices(work.buffers.Zearth)))
        else
            @test all(record->record.context===description, records)
        end
    end
end
