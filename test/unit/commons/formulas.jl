@testitem "Commons / formula options / one option projection for every formula family" tags=[:unit, :engine] begin
    using LineCableModels.Commons: bindings, formulation_options, Expression
    const E = LineCableModels.Engine
    const EH = LineCableModels.Earth.EquivalentHomogeneous
    soil = E.EarthPair(1, 2, (-1.0, -2.0), 1.0, (2, 2))
    self = E.EarthPair(1, 1, (-1.0, -1.0), 0.0, (2, 2); radius = 0.01)

    # One entry per distinct expression, in order of first appearance. Equal expressions
    # share their options.
    selected = E.EarthImpedance.Formula(:default;
        options = (integration = (method = :quad, options = (;)),))
    expressions = [Expression(selected, pair) for pair in (self, soil, self)]
    projected = formulation_options(selected, expressions)
    @test projected.expressions == unique(expressions)
    @test map(expression -> first(expression.arguments), projected.expressions) ==
          [Val(:self), Val(:mutual)]
    for (expression, options) in zip(projected.expressions, projected.options)
        @test options.data.integration.method === Val(:quad)
        @test keys(options.data) == keys(formulation_options(expression).data)
    end

    # An unused supplied section raises an error with the formula identifier.
    unused = E.EarthImpedance.Formula(:default; options = (unknown = 1,))
    @test_throws "unused formulation options (:unknown,) for :$(formula_id(unused))" formulation_options(
        unused, (Expression(unused, soil),))
    internal = E.InternalImpedance.Formula(:default)
    @test_throws "unused formulation options (:unknown,) for :$(formula_id(internal))" E.InternalImpedance.Formula(
        :default; options = (unknown = 1,))

    # Equivalent-earth reductions bind their equivalent_material equations.
    rule = EH.Formula(:bottommost)
    reduced = only(bindings(rule, (soil,)))
    @test reduced.expression.method === EH.equivalent_material
    @test reduced.options == formulation_options(reduced.expression)
end

@testitem "Commons / Functor / one evaluation point of a formula" tags=[:unit, :engine] begin
    using LineCableModels.Commons: Functor, Expression, formulation_options, AbstractFormulation
    const E = LineCableModels.Engine
    # Without a method of its own, a formula does not share values, and its state is empty.
    struct PlainFormula <: AbstractFormulation end
    plain = Functor(PlainFormula(), (; value = 1.0))
    @test plain.formula === PlainFormula() && plain.input.value == 1.0
    @test isempty(plain.state)
    selected = E.EarthImpedance.Formula(:wedepohl1973)
    rho, epsilon, mu = [Inf, 100.0], 8.8541878128e-12 .* [1, 10], fill(4pi*1e-7, 2)
    jω = complex(0.0, 100pi)
    pair = E.EarthPair(1, 2, (-1.0, -2.0), 0.75, (2, 2))
    expression = Expression(selected, pair)
    options = only(formulation_options(selected, (expression,)).options)
    # The earth formula's Functor of one pair, built as the standalone call builds it.
    functor = Functor(selected, (; jω, thickness = nothing, rho, epsilon, mu, options), (;))
    @test functor.formula === selected && functor.input.rho === rho
    # A conductor pair extends the input and keeps the formula and the state.
    point = Functor(functor, (; pair, physical = pair))
    @test point.formula === selected && point.state === functor.state
    @test point.input.pair === pair && point.input.rho === rho
    # The expression evaluates at the point, as the standalone call does.
    @test expression(point, nothing) == selected(rho, epsilon, mu, jω, pair)
    # The existence check matches the Functor and the workspace that the evaluation passes.
    @test validate(expression) === expression
    @test_throws ArgumentError validate(
        Expression(selected, E.EarthPair(1, 2, (10.0, 12.0), 0.75, (1, 1))))
end

@testitem "Commons / formula registry / a family lists its identifiers by its Formula type" tags=[:unit, :engine] begin
    using LineCableModels.Commons: formulas
    const EI = LineCableModels.Engine.EarthImpedance
    registered = formulas(EI.Formula)
    @test registered isa Tuple{Vararg{Symbol}}
    @test :default in registered && :unified in registered
    # Any member type of the family lists the same identifiers.
    @test formulas(typeof(EI.Formula(:unified))) === registered
    @test_throws MethodError formulas(LineParametersFormulation)
end
