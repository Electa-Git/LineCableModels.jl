@testitem "ParametricBuilder / formulation grids / construction and traversal" tags=[:unit] setup=[
    UseEngineSupport,
] begin
    selections=(
        internal_impedance = Grid((:default, :Schelkunoff1934)),
        insulation_impedance = Grid((:default, :Ametani1980)),
        earth_impedance = Grid((:Pollaczek1926, :Papadopoulos2010)),
        insulation_admittance = Grid((:Ametani2004, :Gustavsen2013)),
        semicon_admittance = Grid((:default, :Ametani2004)),
        earth_admittance = Grid((:IdealGround, :Papadopoulos2010)),
        earth_properties = Grid((nothing, :CIGRE2019)),
        equivalent_earth = Grid((
            formula(:Layer; layer = -1),
            formula(:Xue2021)
        )),
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
                    LineCableModels.FormulaSpec
                }),
                values(value.methods)
            ) &&
                !(value.options isa Union{AbstractGrid, Gridspace})
        end
    end

    product=Formulation(
        earth_impedance = Grid((:Pollaczek1926, :Papadopoulos2010)),
        earth_admittance = Grid((:IdealGround, :Papadopoulos2010))
    )
    @test length(product) == 4
    @test Set((
                  formula_id(value.methods.earth_impedance),
                  formula_id(value.methods.earth_admittance)
              ) for value in product) == Set((
        (:Pollaczek1926, :IdealGround),
        (:Papadopoulos2010, :IdealGround),
        (:Pollaczek1926, :Papadopoulos2010),
        (:Papadopoulos2010, :Papadopoulos2010)
    ))

    zipped=Formulation(
        earth_impedance = Grid((:Pollaczek1926, :Papadopoulos2010)),
        earth_admittance = Grid((:IdealGround, :Papadopoulos2010));
        combine = :zip
    )
    @test length(zipped) == 2
    @test [(
               formula_id(value.methods.earth_impedance),
               formula_id(value.methods.earth_admittance)
           ) for value in zipped] == [
        (:Pollaczek1926, :IdealGround),
        (:Papadopoulos2010, :Papadopoulos2010)
    ]

    broadcast_zip=Formulation(
        earth_impedance = Grid((:Pollaczek1926, :Papadopoulos2010)),
        earth_admittance = Grid(:IdealGround);
        combine = :zip
    )
    @test length(broadcast_zip) == 2
    @test [formula_id(value.methods.earth_admittance)
           for value in broadcast_zip] == fill(:IdealGround, 2)

    constants=CableConstantsFormulation(
        internal_impedance = Grid((:Schelkunoff1934, :Ametani2004)),
    )
    @test constants isa Gridspace{CableConstantsFormulation}
    @test length(constants) == 2

    modal=ModalTransformationFormulation(
        Grid((:Fortescue, :Chrysochos2014)),
    )
    @test modal isa Gridspace{ModalTransformationFormulation}
    @test formula_id.(collect(modal)) == [:Fortescue, :Chrysochos2014]

    modal_assumptions=ModalTransformationFormulation(Grid((
        formula(:Fortescue; tolerance = 1e-4),
        formula(:Fortescue; tolerance = 1e-8)
    )))
    @test LineCableModels.Transforms.assumptions.(collect(modal_assumptions)) == [
        (tolerance = 1e-4,),
        (tolerance = 1e-8,)
    ]

    fem=LineCableModelsFEM(
        fem_options = Grid((
        (; mesh_policy = :reuse),
        (; mesh_policy = :remesh)
    )),
    )
    @test fem isa Gridspace{LineCableModelsFEM}
    @test length(fem) == 2

    struct CountedProblem<:AbstractProblemDefinition
        value::Int
    end
    struct CountedFormulation{ID}<:AbstractFormulation end
    struct CountedResult<:AbstractCoreResult
        value::Int
    end
    builds=Ref(0)
    make_problem(value) = (builds[]+=1; CountedProblem(value))
    problem_space=Gridspace{CountedProblem}(make_problem, (Grid((1, 2)),))
    formulas=Grid((CountedFormulation{:first}(), CountedFormulation{:second}()))
    LineCableModels.compute(
        problem::CountedProblem,
        ::CountedFormulation{:first};
        options::NamedTuple = (;)
    ) = CountedResult(problem.value)
    LineCableModels.compute(
        problem::CountedProblem,
        ::CountedFormulation{:second};
        options::NamedTuple = (;)
    ) = CountedResult(10problem.value)

    run=compute(
        ParametricProblem(problem_space),
        Combinatorial(formulas)
    )
    @test builds[] == 2
    @test result(run) == CountedResult.(Int[1, 2, 10, 20])
    @test result(run, 1, 1) == CountedResult(1)
    @test result(run, 2, 1) == CountedResult(2)
    @test result(run, 1, 2) == CountedResult(10)
    @test result(run, 2, 2) == CountedResult(20)
    @test run.axes.problems === problem_space
    @test run.axes.formulations isa Vector{<:AbstractFormulation}
    @test_throws BoundsError result(run, 3, 1)

    linear=compute(
        ParametricProblem(problem_space),
        LinearError(CountedFormulation{:first}())
    )
    @test builds[] == 4
    @test result(linear) == CountedResult.(Int[1, 2])

    legacy=ParametricResult(
        Combinatorial(CountedFormulation{:first}()),
        CountedResult[CountedResult(1)]
    )
    @test isempty(legacy.axes)
    @test_throws ArgumentError result(legacy, 1, 1)
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
