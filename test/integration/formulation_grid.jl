@testitem "Engine / formulation grids / exact batched calculations" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures
] begin
    function same_parameters(left, right)
        same_domain=if domain(left)===ModalDomain
            formula_id(left.domain.formula)==formula_id(right.domain.formula)&&
                operators(left).voltage==operators(right).voltage&&
                operators(left).current==operators(right).current
        else
            left.domain==right.domain
        end
        return same_domain&&
               left.Z.values==right.Z.values&&
               left.Y.values==right.Y.values&&
               left.f==right.f&&
               left.details==right.details
    end

    system=TestFixtures.three_phase_system()
    earth=Grid((
        EarthModel(10.0, 10.0, 1.0),
        EarthModel(100.0, 10.0, 1.0)
    ))
    problem_space=LineParametersProblem(
        system,
        earth;
        frequencies = [50.0]
    )
    formulation_space=Formulation(
        earth_impedance = Grid((:Pollaczek1926, :Papadopoulos2010)),
    )
    problems=collect(problem_space)
    formulations=collect(formulation_space)
    expected=[compute(problem, formulation)
              for formulation in formulations
              for problem in problems]

    run=compute(
        ParametricProblem(problem_space),
        Combinatorial(formulation_space; options = (retain_details = true,))
    )
    @test length(run) == length(problems) * length(formulations)
    @test length(details(run).points) == length(run)
    @test formula_id.(getproperty.(
        getproperty.(run.axes.formulations, :methods),
        :earth_impedance
    )) == [:Pollaczek1926, :Papadopoulos2010]
    for index in eachindex(expected)
        @test same_parameters(run[index], expected[index])
    end
    for formulation_index in eachindex(formulations)
        for problem_index in eachindex(problems)
            index=problem_index+
            (formulation_index-1)*length(problems)
            @test same_parameters(
                run[problem_index, formulation_index],
                expected[index]
            )
        end
    end

    direct=compute(first(problems), formulations)
    @test length(direct) == length(formulations)
    for formulation_index in eachindex(formulations)
        @test same_parameters(
            direct[formulation_index],
            compute(first(problems), formulations[formulation_index])
        )
    end

    design=TestFixtures.mv_cable_design()
    constants_problem=CableConstantsProblem(design; frequency = 50.0)
    constants_space=CableConstantsFormulation(
        insulation_admittance = Grid((:Ametani2004, :Gustavsen2013)),
    )
    constants_formulations=collect(constants_space)
    constants_batch=compute(constants_problem, constants_formulations)
    @test constants_batch == [compute(constants_problem, formulation)
           for formulation in constants_formulations]

    phase=first(expected)
    modal_problem=ModalTransformationProblem(phase)
    modal_space=ModalTransformationFormulation(
        Grid((:Fortescue, :Chrysochos2014)),
    )
    modal_formulations=collect(modal_space)
    modal_batch=compute(modal_problem, modal_formulations)
    modal_scalar=[compute(modal_problem, formulation)
                  for formulation in modal_formulations]
    @test length(modal_batch) == 2
    @test isconcretetype(eltype(modal_batch))
    @test typeof(modal_batch[1]) === typeof(modal_batch[2])
    @test fieldtype(typeof(modal_batch[1].domain), :formula) ===
          LineCableModels.Transforms.Formula
    for index in eachindex(modal_batch)
        @test same_parameters(modal_batch[index], modal_scalar[index])
    end

    transported=Gridspace{ModalTransformationProblem}(run)
    @test length(transported) == length(run)
    @test transported.grids === (run,)
    modal_run=compute(
        ParametricProblem(transported),
        Combinatorial(modal_space)
    )
    expected_modal=[compute(ModalTransformationProblem(parameters), formulation)
                    for formulation in modal_formulations
                    for parameters in run]
    @test length(modal_run) == length(run) * length(modal_formulations)
    @test all(
        same_parameters(modal_run[index], expected_modal[index])
    for index in eachindex(expected_modal)
    )
end

@testitem "Engine / formulation grids / one lowering per selected design" tags=[:integration] begin
    import LineCableModels.Engine as EN

    # Specialize only designs bearing a test-owned nominal-data type. invoke
    # then runs the unmodified production lowering; no global method is replaced.
    struct LoweringCounter
        calls::Base.RefValue{Int}
    end
    function EN.flatten(engine::LineCableModelsCoaxial,
            design::CableDesign{T, R, G, NamedTuple{(:counter,), Tuple{LoweringCounter}}},
            ::Type{S}) where {
            T<:Real, R<:AbstractCablePart,
            G<:LineCableModels.DataModel.CableGeometry, S<:Real
    }
        design.nominal_data.counter.calls[]+=1
        return invoke(EN.flatten,
            Tuple{LineCableModelsCoaxial, CableDesign, Type{S}}, engine, design, S)
    end

    counter=LoweringCounter(Ref(0))
    copper=Material(kind=:conductor, rho=1.7241e-8)
    dielectric=Material(kind=:insulator, rho=1e14, eps_r=2.3)
    design=build(CableDesign, "counted-lowering",
        terminal(:core, core(copper; r=5e-3), insulation(dielectric; t=2e-3)),
        terminal(:sheath, sheath(copper; t=0.5e-3));
        nominal_data=(; counter))
    system=build(LineCableSystem, design, Pose2(0.0, -1.0);
        connections=(core=1, sheath=0))
    problems=LineParametersProblem(system, homogeneous(rho=Grid((10.0, 100.0)));
        frequencies=[1.0, 50.0, 1000.0])
    formulations=Formulation(earth_impedance=Grid((:Pollaczek1926, :Papadopoulos2010)))
    @test counter.calls[] == 0

    phase=compute(ParametricProblem(problems), Combinatorial(formulations))
    @test length(phase) == 4
    @test counter.calls[] == 2
    @test all(parameters -> domain(parameters) === PhaseDomain, phase)
    @test all(parameters -> size(parameters.Z) == (1, 1, 3), phase)

    counter.calls[]=0
    scalar=compute(first(problems), first(formulations))
    @test counter.calls[] == 1
    @test scalar.Z.values == first(phase).Z.values
    @test scalar.Y.values == first(phase).Y.values

    counter.calls[]=0
    constants_problem=CableConstantsProblem(design; frequency=50.0)
    constants_formulations=collect(CableConstantsFormulation(
        insulation_admittance=Grid((:Ametani2004, :Gustavsen2013))))
    constants=compute(constants_problem, constants_formulations)
    @test length(constants) == 2
    @test counter.calls[] == 1

    counter.calls[]=0
    modal_problems=Gridspace{ModalTransformationProblem}(phase)
    modal_formulations=ModalTransformationFormulation(
        Grid((:Fortescue, :Chrysochos2014)))
    modal=compute(ParametricProblem(modal_problems), Combinatorial(modal_formulations))
    @test length(modal) == 8
    @test all(parameters -> domain(parameters) === ModalDomain, modal)
    @test counter.calls[] == 0
end
