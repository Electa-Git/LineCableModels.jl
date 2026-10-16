@testitem "LineModels / formula registry and formulation grammar" tags=[:unit, :engine] begin
    import LineCableModels.LineModels as LM

    @test LM.formulas(LM.Formula) == (:ulm, :fdcm, :marti, :default)
    @test formula_id(LineModelFormulation().formula) === :ulm
    for identifier in LM.formulas(LM.Formula)
        selected = LM.Formula(identifier)
        @test formula_id(selected) === (identifier === :default ? :ulm : identifier)
        @test !isempty(description(selected))
        @test LM.Formula(selected) === selected
        @test NamedTuple(selected).identifier === formula_id(selected)
        @test formula_id(LM.Formula(formula(identifier))) === formula_id(selected)
    end
    @test_throws ArgumentError LM.Formula(:unknown)
    @test_throws ArgumentError LineModelFormulation(:unknown)
    for identifier in (:ulm, :fdcm)
        controls = (yc=(npoles=8,), h=(npoles=6,), passivity=(enabled=true,), dc=(;))
        selected = LineModelFormulation(identifier; options=controls)
        @test formulation_options(selected.formula).data === controls
        @test isempty(formulation_options(LineModelFormulation(identifier).formula).data)
        @test_throws ArgumentError LineModelFormulation(identifier; options=(unknown=(;),))
        @test_throws ArgumentError LM.Formula(identifier; parameters=(unknown=1,))
    end
    marti = formulation_options(LineModelFormulation(:marti).formula).data
    @test marti == (transformation=(frequency=nothing,), bode=(tolerance_db=3.0,))
    supplied = LineModelFormulation(:marti; options=(transformation=(frequency=50.0,),))
    @test supplied.formula.options.data.transformation.frequency == 50.0
    @test_throws ArgumentError LineModelFormulation(:marti; options=(unknown=(;),))
    @test_throws ArgumentError LineModelFormulation(:marti; options=(bode=(unknown=1,),))
    selected = LineModelFormulation(:ulm)
    @test LineModelFormulation(selected) === selected
    @test_throws ArgumentError LineModelFormulation(selected; options=(yc=(;),))
    @test_throws ArgumentError LineModelFormulation(formula(:ulm); options=(yc=(;),))
    space = LineModelFormulation(Grid((:ulm, :fdcm, :marti)))
    @test space isa Gridspace{LineModelFormulation}
    @test [formula_id(value.formula) for value in space] == [:ulm, :fdcm, :marti]
end
