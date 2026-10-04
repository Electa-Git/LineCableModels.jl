@testitem "Commons / formula bindings / one option projection for every formula family" tags=[:unit, :engine] begin
    using LineCableModels.Commons: bindings, formulation_options
    const E = LineCableModels.Engine
    const EH = LineCableModels.Earth.EquivalentHomogeneous
    soil = E.EarthPair(1, 2, (-1.0, -2.0), 1.0, (2, 2))
    self = E.EarthPair(1, 1, (-1.0, -1.0), 0.0, (2, 2); radius = 0.01)

    # One record per interaction, in order. Equal equations share their options.
    selected = E.EarthImpedance.Formula(:default;
        options = (integration = (method = :quad, options = (;)),))
    records = bindings(selected, [self, soil, self])
    @test length(records) == 3
    @test records[1] == records[3]
    @test (records[1].kind, records[2].kind) == (:self, :mutual)
    @test records[1].options === records[3].options
    for record in records
        @test record.options.data.integration.method === Val(:quad)
        @test keys(record.options.data) == keys(formulation_options(record.equation).data)
    end
    @test bindings(selected, (self, soil)) isa Tuple

    # An unused supplied section raises an error with the formula identifier.
    unused = E.EarthImpedance.Formula(:default; options = (unknown = 1,))
    @test_throws "unused formulation options (:unknown,) for :$(formula_id(unused))" bindings(
        unused, (soil,))
    internal = E.InternalImpedance.Formula(:default)
    @test_throws "unused formulation options (:unknown,) for :$(formula_id(internal))" E.InternalImpedance.Formula(
        :default; options = (unknown = 1,))

    # Equivalent-earth reductions bind their equivalent_material equations.
    rule = EH.Formula(:bottommost)
    reduced = only(bindings(rule, (soil,)))
    @test reduced.equation.method === EH.equivalent_material
    @test reduced.options == formulation_options(reduced.equation)
end
