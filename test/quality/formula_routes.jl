@testitem "Quality / native equation bindings and closed built-in catalogues" tags=[:quality] setup=[FormulaContractModels] begin
    const E=LineCableModels.Engine
    const EP=LineCableModels.Earth
    const FM=LineCableModels.FormulaMethod
    owners=(E.InternalImpedance,E.InsulationImpedance,E.EarthImpedance,
        E.InsulationAdmittance,E.SemiconAdmittance,E.EarthAdmittance,
        EP.FrequencyDependent,EP.EquivalentHomogeneous,LineCableModels.Transforms,
        LineCableModels.Materials.TemperatureDependent)
    scalar_operations=(E.InsulationImpedance=>E.InsulationImpedance.insulation_impedance,
        E.InsulationAdmittance=>E.InsulationAdmittance.insulation_material,
        E.SemiconAdmittance=>E.SemiconAdmittance.semicon_material,
        EP.FrequencyDependent=>EP.FrequencyDependent.earth_material,
        LineCableModels.Transforms=>LineCableModels.Transforms.modal_operators,
        LineCableModels.Materials.TemperatureDependent=>
            LineCableModels.Materials.TemperatureDependent.temperature_resistivity)
    for owner in owners, identifier in owner.formulas()
        selected=owner.Formula(identifier)
        @test !hasfield(typeof(selected), :hooks)
        @test !hasfield(typeof(selected), :binding)
        @test fieldtype(typeof(selected), :options) <: FormulationOptions
        @test formulation_options(selected) === selected.options
        @test formula_id(selected) !== :default
        @test owner.Formula(selected) === selected
        @test_throws MethodError owner.Formula(identifier; hooks=(;))
        routes = if owner in (E.EarthImpedance,E.EarthAdmittance)
            bindings=FM[]
            for kind in (:self,:mutual), source in 1:2, target in 1:2
                kind === :self && source != target && continue
                pair=E.EarthPair(1,kind === :self ? 1 : 2,
                    (source == 1 ? 1.0 : -1.0,target == 1 ? 1.0 : -1.0),
                    kind === :self ? 0.0 : 1.0,(source,target);
                    radius=kind === :self ? 0.01 : nothing)
                binding=FM(selected,pair)
                signature=Tuple{typeof(selected),typeof.(binding.arguments)...,Any,Any,Any}
                which(binding.method,signature) === owner.EQUATION_FALLBACK && continue
                push!(bindings,validate(selected,pair).equation)
            end
            @test !isempty(bindings)
            bindings
        elseif owner === E.InternalImpedance
            Tuple(FM(selected,owner.internal_impedance,Val(kind)) for kind in (:inner,:outer,:transfer))
        elseif owner === EP.EquivalentHomogeneous
            (validate(selected,E.EarthPair(1,1,(1.0,1.0),0.0,(1,1);radius=0.01)).equation,)
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
    M=FormulaContractModels
    for (owner,custom) in ((E.InternalImpedance,M.SurfaceLaw()),
            (E.InsulationImpedance,M.InsulationReactance()),
            (E.InsulationAdmittance,M.InsulationLaw()),(E.SemiconAdmittance,M.SemiconLaw()),
            (E.EarthImpedance,M.selection(E.EarthImpedance)),
            (E.EarthAdmittance,M.selection(E.EarthAdmittance)),
            (EP.FrequencyDependent,M.DispersiveEarth()),
            (EP.EquivalentHomogeneous,M.MeanEarth()),
            (LineCableModels.Materials.TemperatureDependent,M.ConstantResistivity(1e-8)),
            (E.ShuntModel,M.UserCoaxialShunt()),(E.PipeImpedance,M.CoaxialPipePolicy()))
        @test owner.Formula(custom) === custom
        @test formula_id(custom) ∉ owner.formulas()
    end
    @test_throws MethodError formula(:default;hooks=(;))
end

@testitem "Quality / user-owned shunt and pipe selections reach blueprint and compute" tags=[:quality] setup=[TestFixtures,FormulaContractModels] begin
    M=FormulaContractModels
    shunt=M.UserCoaxialShunt()
    selected=Formulation(shunt_model=shunt,pipe_impedance=M.CoaxialPipePolicy())
    problem=TestFixtures.line_parameters_problem(frequencies=[50.0,500.0])
    actual=compute(problem,selected)
    expected=compute(problem,Formulation())
    @test Z(actual)==Z(expected) && Y(actual)==Y(expected)
    @test shunt.preparations[]==length(problem.system.designs)
    @test details(actual).data.formulations.effective.shunt_model === :UserCoaxialShunt
    @test details(actual).data.formulations.effective.pipe_impedance === :CoaxialPipePolicy
    constants=CableConstantsProblem(first(problem.system.designs);frequency=50.0)
    @test compute(constants,CableConstantsFormulation(shunt_model=shunt,pipe_impedance=M.CoaxialPipePolicy()))==
        compute(constants,CableConstantsFormulation())
    source=TestFixtures.two_conductor_results()
    # Modal selections are admitted through the modal action, not a callback bag.
    maps=operators(compute(ModalTransformationProblem(source),ModalTransformationFormulation()))
    custom=M.FixedModalMaps(maps.voltage,maps.current)
    action=ModalTransformationFormulation(custom)
    @test action.formula === custom
    @test action.definition === custom
    @test all(isfinite,Z(compute(ModalTransformationProblem(source),action)))
end
