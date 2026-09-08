@testitem "ParametricBuilder / formulation grids / construction and traversal" tags=[:unit] setup=[
    UseEngineSupport,
] begin
    inner=(f, w)->zero(f.state.jω)
    insulation=(r_in, r_ex, mu_r, s, values, options, workspace)->zero(s)
    @eval LineCableModels.computation_options(
        ::LineCableModels.FormulaMethod{
            :default, typeof(LineCableModels.Engine.InternalImpedance.internal_impedance),
            Tuple{Val{:inner}}},
        ::$(typeof(inner))) = (;)
    @eval LineCableModels.computation_options(
        ::LineCableModels.FormulaMethod{:default,
            typeof(LineCableModels.Engine.InsulationImpedance.insulation_impedance)},
        ::$(typeof(insulation))) = (;)
    selections=(
        internal_impedance = Grid((
            :default, formula(:default; hooks = (inner = inner,)))),
        insulation_impedance = Grid((:default,
            formula(:default; hooks = (contribution = insulation,)))),
        earth_impedance = Grid((:Pollaczek1926, :default)),
        insulation_admittance = Grid((:Ametani2004, :default)),
        semicon_admittance = Grid((:default, :Ametani2004)),
        earth_admittance = Grid((:IdealGround, :default)),
        earth_properties = Grid((nothing, :default)),
        pipe_impedance = Grid((:default, formula(:default))),
        options = Grid((
            (; ideal_transposition = false),
            (; ideal_transposition = true)
        ))
    )
    for name in keys(selections)
        keyword=NamedTuple{(name,)}((getproperty(selections, name),))
        space=Formulation(; keyword...)
        @test space isa Gridspace{LineParametersFormulation}
        @test length(space) == 2
        @test all(value -> value isa LineParametersFormulation, space)
        @test all(space) do value
            all(
                item -> !(item isa Union{
                    AbstractGrid,
                    Gridspace,
                    Symbol,
                    LineCableModels.FormulaDefinition
                }),
                values(value.methods)
            ) &&
                !(value.options isa Union{AbstractGrid, Gridspace})
        end
    end

    product=Formulation(
        earth_impedance = Grid((:Pollaczek1926, :default)),
        earth_admittance = Grid((:IdealGround, :default))
    )
    @test length(product) == 4
    @test Set((
                  formula_id(value.methods.earth_impedance),
                  formula_id(value.methods.earth_admittance)
              ) for value in product) == Set((
        (:Pollaczek1926, :IdealGround),
        (:default, :IdealGround),
        (:Pollaczek1926, :default),
        (:default, :default)
    ))

    zipped=Formulation(
        earth_impedance = Grid((:Pollaczek1926, :default)),
        earth_admittance = Grid((:IdealGround, :default));
        combine = :zip
    )
    @test length(zipped) == 2
    @test [(
               formula_id(value.methods.earth_impedance),
               formula_id(value.methods.earth_admittance)
           ) for value in zipped] == [
        (:Pollaczek1926, :IdealGround),
        (:default, :default)
    ]

    broadcast_zip=Formulation(
        earth_impedance = Grid((:Pollaczek1926, :default)),
        earth_admittance = Grid(:IdealGround);
        combine = :zip
    )
    @test length(broadcast_zip) == 2
    @test [formula_id(value.methods.earth_admittance)
           for value in broadcast_zip] == fill(:IdealGround, 2)

    reductions=Formulation(earth_impedance = Grid((
        formula(:default; equivalent_earth = formula(:default; order = :before)),
        formula(:default; equivalent_earth = formula(:default; order = :after)))))
    @test length(reductions) == 2
    @test first(reductions).methods.earth_impedance.equivalent_earth isa
          LineCableModels.Earth.EquivalentHomogeneous.BeforeFD
    @test last(collect(reductions)).methods.earth_impedance.equivalent_earth isa
          LineCableModels.Earth.EquivalentHomogeneous.AfterFD
    @test all(value -> value.methods.earth_admittance.equivalent_earth === nothing, reductions)

    constants=CableConstantsFormulation(
        insulation_admittance = Grid((:default, :Ametani2004)),
    )
    @test constants isa Gridspace{CableConstantsFormulation}
    @test length(constants) == 2

    modal=ModalTransformationFormulation(
        Grid((:default, :default)),
    )
    @test modal isa Gridspace{ModalTransformationFormulation}
    @test formula_id.(collect(modal)) == [:default, :default]

    modal_assumptions=ModalTransformationFormulation(Grid((
        formula(:default; options = (iteration = (convergence = 1e-4,),)),
        formula(:default; options = (iteration = (convergence = 1e-8,),))
    )))
    @test [value.formula.options.iteration.convergence for value in modal_assumptions] ==
          [1e-4, 1e-8]

    fem=LineCableModelsFEM(
        fem_options = Grid((
        (; mesh_policy = :reuse),
        (; mesh_policy = :remesh)
    )),
    )
    @test fem isa Gridspace{LineCableModelsFEM}
    @test length(fem) == 2

    for name in keys(selections)
        keyword=NamedTuple{(name,)}((getproperty(selections, name),))
        space=Formulation(:LineCableModelsFEM; keyword...)
        @test space isa Gridspace{LineCableModelsFEM}
        @test length(space) == 2
        @test all(value -> isconcretetype(typeof(value)), space)
        if name!==:options
            @test all(value -> haskey(value.definitions, name), space)
        end
    end
    fem_zipped=Formulation(:LineCableModelsFEM;
        insulation_admittance = Grid((:default, :Ametani2004)),
        semicon_admittance = Grid((:default, :Ametani2004)), combine = :zip)
    @test length(fem_zipped) == 2
    @test [(formula_id(value.methods.insulation_admittance),
               formula_id(value.methods.semicon_admittance)) for value in fem_zipped] ==
          [(:default, :default), (:Ametani2004, :Ametani2004)]

    struct CountedProblem<:AbstractProblemDefinition
        value::Int
    end
    LineCableModels.validate(problem::CountedProblem)=problem
    struct CountedFormulation{ID}<:AbstractFormulation end
    struct CountedResult<:AbstractCoreResult
        value::Int
    end
    builds=Ref(0)
    make_problem(value)=(builds[]+=1; CountedProblem(value))
    problem_space=Gridspace{CountedProblem}(make_problem, (Grid((1, 2)),))
    formulas=Grid((CountedFormulation{:first}(), CountedFormulation{:second}()))
    LineCableModels.compute(
        problem::CountedProblem,
        ::CountedFormulation{:first};
        options::NamedTuple = (;)
    )=CountedResult(problem.value)
    LineCableModels.compute(
        problem::CountedProblem,
        ::CountedFormulation{:second};
        options::NamedTuple = (;)
    )=CountedResult(10problem.value)

    run=compute(
        ParametricProblem(problem_space),
        Combinatorial(formulas)
    )
    @test builds[] == 2
    @test collect(run) == CountedResult.(Int[1, 2, 10, 20])
    @test run[1, 1] == CountedResult(1)
    @test run[2, 1] == CountedResult(2)
    @test run[1, 2] == CountedResult(10)
    @test run[2, 2] == CountedResult(20)
    @test run.axes.problems === problem_space
    @test run.axes.formulations isa Vector{<:AbstractFormulation}
    @test_throws BoundsError run[3, 1]

    # A completed scalar problem is an admitted singleton, not an iterable object.
    scalar=CountedProblem(7)
    scalar_problem=ParametricProblem(scalar)
    @test length(scalar_problem.space) == 1
    @test first(scalar_problem.space) === scalar
    scalar_run=compute(scalar_problem, Combinatorial(formulas))
    @test collect(scalar_run) == CountedResult.(Int[7, 70])
    @test scalar_run[1, 2] == CountedResult(70)

    linear=compute(
        ParametricProblem(problem_space),
        LinearError(CountedFormulation{:first}())
    )
    @test builds[] == 4
    @test collect(linear) == CountedResult.(Int[1, 2])

    legacy=ParametricResult(
        Combinatorial(CountedFormulation{:first}()),
        CountedResult[CountedResult(1)]
    )
    @test isempty(legacy.axes)
    @test_throws ArgumentError legacy[1, 1]
    @test_throws MethodError compute(CountedProblem(1), formulas)
    @test_throws ArgumentError Combinatorial(Grid((1, 2)))
    @test_throws ArgumentError Combinatorial(
        Gridspace{Int}(identity, (Grid((1, 2)),)),
    )

    struct UnsupportedFormulation<:AbstractFormulation end
    @test_throws MethodError compute(
        CountedProblem(1),
        UnsupportedFormulation[UnsupportedFormulation()]
    )
end
