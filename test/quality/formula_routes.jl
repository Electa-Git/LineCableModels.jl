@testitem "Quality / native equation bindings and closed built-in formula lists" tags=[:quality] setup=[FormulaFixtures] begin
    const E=LineCableModels.Engine
    const EP=LineCableModels.Earth
    const FM=LineCableModels.FormulaMethod
    owners=(E.InternalImpedance,E.InsulationImpedance,E.EarthImpedance,
        E.InsulationAdmittance,E.SemiconAdmittance,E.EarthAdmittance,
        EP.FrequencyDependent,EP.EquivalentHomogeneous,LineCableModels.ModalAnalysis,
        LineCableModels.Materials.TemperatureDependent)
    scalar_operations=(E.InsulationImpedance=>E.InsulationImpedance.insulation_impedance,
        E.InsulationAdmittance=>E.InsulationAdmittance.insulation_material,
        E.SemiconAdmittance=>E.SemiconAdmittance.semicon_material,
        EP.FrequencyDependent=>EP.FrequencyDependent.earth_material,
        LineCableModels.ModalAnalysis=>LineCableModels.ModalAnalysis.decompose!,
        LineCableModels.Materials.TemperatureDependent=>
            LineCableModels.Materials.TemperatureDependent.temperature_resistivity)
    for owner in owners, identifier in owner.formulas()
        selected=owner.Formula(identifier)
        @test fieldtype(typeof(selected), :options) <: FormulationOptions
        @test formulation_options(selected) === selected.options
        @test formula_id(selected) !== :default
        @test owner.Formula(selected) === selected
        routes = if owner in (E.EarthImpedance,E.EarthAdmittance)
            bindings=FM[]
            for kind in (:self,:mutual), source in 1:2, target in 1:2
                kind === :self && source != target && continue
                pair=E.EarthPair(1,kind === :self ? 1 : 2,
                    (source == 1 ? 1.0 : -1.0,target == 1 ? 1.0 : -1.0),
                    kind === :self ? 0.0 : 1.0,(source,target);
                    radius=kind === :self ? 0.01 : nothing)
                binding=FM(selected,pair)
                # Registration and indexed binding do not claim an implemented
                # equation. Actual supported, stub and unsupported calls are covered
                # by the execution tests, not a reflected coverage list.
                push!(bindings,binding)
            end
            bindings
        elseif owner === E.InternalImpedance
            Tuple(FM(selected,owner.internal_impedance,Val(kind)) for kind in (:inner,:outer,:transfer))
        elseif owner === EP.EquivalentHomogeneous
            (only(LineCableModels.Commons.bindings(selected,(E.EarthPair(1,1,(1.0,1.0),0.0,(1,1);radius=0.01),))).equation,)
        else
            operation=only(last(entry) for entry in scalar_operations if first(entry) === owner)
            (FM(selected,operation),)
        end
        for route in routes
            @test route.selection === selected
            @test typeof(route).parameters[1] === typeof(selected)
            @test parentmodule(route.method) === owner
            @test all(arg->arg isa Val,route.arguments)
            @test formulation_options(route) isa FormulationOptions
            @test_throws MethodError computation_options(route)
            @test_throws MethodError formulation_options(route, ComputationOptions())
        end
    end
    M=FormulaFixtures
    for (owner,custom) in ((E.InternalImpedance,M.SurfaceLaw()),
            (E.InsulationImpedance,M.InsulationReactance()),
            (E.InsulationAdmittance,M.InsulationLaw()),(E.SemiconAdmittance,M.SemiconLaw()),
            (E.EarthImpedance,M.selection(E.EarthImpedance)),
            (E.EarthAdmittance,M.selection(E.EarthAdmittance)),
            (EP.FrequencyDependent,M.DispersiveEarth()),
            (EP.EquivalentHomogeneous,M.MeanEarth()),
            (LineCableModels.Materials.TemperatureDependent,M.ConstantResistivity(1e-8)),
            (E.ShuntModel,M.UserCoaxialShunt()),(E.PipeImpedance,M.UserCoaxialPipe()))
        @test owner.Formula(custom) === custom
        @test formula_id(custom) ∉ owner.formulas()
    end
end

@testitem "Quality / user-owned shunt and pipe selections reach blueprint and compute" tags=[:quality] setup=[TestFixtures,FormulaFixtures] begin
    M=FormulaFixtures
    shunt=M.UserCoaxialShunt()
    selected=Formulation(shunt_model=shunt,pipe_impedance=M.UserCoaxialPipe())
    problem=TestFixtures.line_parameters_problem(frequencies=[50.0,500.0])
    actual=compute(problem,selected)
    expected=compute(problem,Formulation())
    @test Z(actual)==Z(expected) && Y(actual)==Y(expected)
    @test shunt.response_count[]==length(problem.system.designs)
    @test details(actual).data.formulations.methods.shunt_model.identifier === :UserCoaxialShunt
    @test details(actual).data.formulations.methods.pipe_impedance.identifier === :UserCoaxialPipe
    constants=CableConstantsProblem(first(problem.system.designs);frequency=50.0)
    @test compute(constants,CableConstantsFormulation(shunt_model=shunt,pipe_impedance=M.UserCoaxialPipe()))==
        compute(constants,CableConstantsFormulation())
    source=TestFixtures.two_conductor_results()
    # Modal selections are admitted through the modal action, not a callback bag.
    maps=operators(compute(ModalAnalysisProblem(source),ModalAnalysisFormulation()))
    custom=M.FixedModalMaps(maps.Tv,maps.Ti)
    action=ModalAnalysisFormulation(custom)
    @test action.formula === custom
    @test action.definition === custom
    @test all(isfinite,Z(compute(ModalAnalysisProblem(source),action)))
end
