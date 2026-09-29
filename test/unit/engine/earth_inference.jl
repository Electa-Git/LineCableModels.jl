@testitem "Engine / prepared earth tuples preserve equation types across layouts" tags=[:unit] setup=[FormulaFixtures] begin
    const E=LineCableModels.Engine
    design=build(CableDesign, "typed earth equations", terminal(:core,
        core(Material(:conductor, 1.72e-8); r=0.004),
        insulation(Material(:insulator, Inf, 2.3); t=0.002)))
    selections=(Formulation(),
        Formulation(earth_impedance=(air=:gary1976, earth=:saad1996, mixed=:lucca1994),
            earth_admittance=:ideal),
        Formulation(earth_impedance=FormulaFixtures.selection(E.EarthImpedance; layers=2:2),
            earth_admittance=FormulaFixtures.selection(E.EarthAdmittance; layers=2:2)))
    execution=E.computation_options(LineCableModelsCoaxial, ComputationOptions())
    for selected in selections, heights in ((1.0, 1.0, 1.0), (1.0, -1.0, -2.0))
        system=build(LineCableSystem, fill(design, 3),
            [Pose2(Float64(i-1), height) for (i, height) in enumerate(heights)];
            connections=[Dict(:core=>i) for i in 1:3])
        problem=LineParametersProblem(system; frequencies=[50.0], earth_props=homogeneous(rho=100.0))
        blueprints=only(E.flatten(LineCableModelsCoaxial(), problem.system.designs,
            Float64, [selected]))
        workspace=E.LineParametersWorkspace(problem, selected, execution, blueprints)
        calculations=workspace.invariants.earth_calculations
        materials=workspace.buffers.earth_materials
        @test calculations isa Tuple
        @test materials isa Tuple
        @test length(calculations)==length(materials)
        @test all(isconcretetype, fieldtypes(typeof(calculations)))
        for calculation in calculations
            @test calculation.equations isa Tuple
            @test all(isconcretetype, fieldtypes(typeof(calculation.equations)))
            for group in calculation.equations
                @test isconcretetype(fieldtype(typeof(group), :declaration))
                @test isconcretetype(fieldtype(typeof(group.declaration), :equation))
            end
        end
        @test (@inferred E._solve!(workspace, selected, calculations, materials)) === workspace
        @test all(isfinite, workspace.buffers.Zout)
        @test all(isfinite, workspace.buffers.Yout)
    end
end

@testitem "Engine / fixed public options retain phase and modal result types" tags=[:unit] setup=[TestFixtures] begin
    function fixed_options(problem, selected, ::Val{Basis}, ::Val{Trace}, ::Val{Timing}) where {Basis, Trace, Timing}
        compute(problem, selected;
            options=(output_basis=Basis, trace=Trace, timing=Timing, verbosity=(default=0,)))
    end
    function phase_and_modal(problem, selected)
        phase=compute(problem, selected; options=(output_basis=:pul, trace=false, timing=false))
        modal=compute(ModalAnalysisProblem(phase), ModalAnalysisFormulation(:default);
            options=(rotate=true,))
        return phase, modal
    end
    problem=TestFixtures.line_parameters_problem(TestFixtures.two_wire_system(); frequencies=[50.0, 500.0])
    selected=Formulation(options=(reduce_bundle=false, kron_reduction=false, ideal_transposition=false))
    reference=@inferred fixed_options(problem, selected, Val(:pul), Val(false), Val(false))
    for basis in (:pul, :total), trace in (false, true), timing in (false, true)
        result=@inferred fixed_options(problem, selected, Val(basis), Val(trace), Val(timing))
        multiplier=basis===:pul ? 1.0 : problem.system.line_length
        @test Z(result)==multiplier.*Z(reference)
        @test Y(result)==multiplier.*Y(reference)
        @test frequencies(result)==frequencies(reference)
        @test haskey(details(result).data, :trace)==trace
        @test haskey(details(result).data, :timing)==timing
    end
    phase, modal=@inferred phase_and_modal(problem, selected)
    @test Z(phase)==Z(reference)
    restored=transform(PhaseDomain, modal)
    @test Z(restored)≈Z(phase)
    @test Y(restored)≈Y(phase)
end
