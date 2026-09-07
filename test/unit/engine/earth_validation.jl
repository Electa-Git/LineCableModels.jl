@testitem "Engine / earth formulas / owner-defined geometry validation" tags=[:unit] begin
    using RequiredInterfaces
    const EN = LineCableModels.Engine
    const EI = EN.EarthImpedance
    const EA = EN.EarthAdmittance
    const FM = LineCableModels.FormulaMethod

    underground = EN.EarthPair(1, 2, (-1.0, -1.2), 1.0, (2, 2))
    overhead = EN.EarthPair(1, 2, (10.0, 12.0), 1.0, (1, 1))
    vertical = EN.EarthPair(1, 2, (-1.0, -1.2), 0.0, (2, 2))
    self = EN.EarthPair(1, 1, (-1.0, -1.0), 0.02, (2, 2))
    @test @inferred(validate(underground)) === underground
    @test validate(vertical) === vertical
    @test validate(self) === self
    @test_throws DomainError validate(EN.EarthPair(1, 1, (-1.0, -1.0), 0.0, (2, 2)))
    @test_throws DomainError validate(EN.EarthPair(1, 2, (-1.0, -1.0), 0.0, (2, 2)))
    @test_throws ArgumentError validate(EN.EarthPair(1, 2, (1.0, -1.0), 1.0, (2, 2)))

    for (owner, abstract_owner) in ((EI, EN.EarthImpedanceFormulation),
            (EA, EN.EarthAdmittanceFormulation))
        @test RequiredInterfaces.isInterface(abstract_owner)
        @test_throws ArgumentError validate(owner.Formula(:default), underground)
        for id in owner.formulas()
            formula = owner.Formula(id)
            for (name, route) in pairs(formula.routes)
                name === :Γ && continue
                method = which(validate, Tuple{typeof(underground), typeof(route), typeof(formula)})
                # Every native selection owns an explicit method. New catalogue
                # entries must not inherit the custom-callable geometry-only route.
                @test method.module === owner
                @test basename(String(method.file)) == lowercase(string(id)) * ".jl"
            end
        end
    end

    for (owner, good, bad, id) in ((EI, overhead, underground, :Wise1934),
            (EI, underground, overhead, :Xue2018),
            (EA, overhead, underground, :Wise1948),
            (EA, underground, overhead, :Xue2018))
        formula = owner.Formula(id)
        @test @inferred(validate(formula, good)) === formula
        @test_throws ArgumentError validate(formula, bad)
    end
    for id in (:Bridges1995, :Vance1978)
        formula = EI.Formula(id)
        @test validate(formula, self) === formula
        @test_throws ArgumentError validate(formula, underground)
    end
    for (owner, ids) in ((EI, (:Petrache2005, :Saad1996, :Theethayi2007, :WedepohlWilcox1973)),
            (EA, (:Theethayi2007, :Xue2021)))
        for id in ids
            formula = owner.Formula(id)
            @test validate(formula, self) === formula
            @test validate(formula, underground) === formula
            @test_throws DomainError validate(formula, vertical)
        end
    end
    for owner in (EI, EA)
        formula = owner.Formula(:Papadopoulos2011)
        @test validate(formula, underground) === formula
        @test_throws ArgumentError validate(formula,
            EN.EarthPair(1, 2, (-1.0, -2.0), 1.0, (2, 3)))
    end

    # Validation never calls the numerical override, and respects its identity.
    calls = Ref(0)
    custom = (functor, pair) -> (calls[] += 1; error("numerical sentinel"))
    supplied = EI.Formula(:Bridges1995; mutual=custom)
    @test validate(supplied, vertical) === supplied
    @test calls[] == 0
    replacement = EI.Formula(:Bridges1995; mutual=EI.Formula(:Xue2018).routes.mutual)
    @test validate(replacement, vertical) === replacement
    restricted = EI.Formula(:Xue2018; mutual=EI.Formula(:Bridges1995).routes.mutual)
    @test_throws ArgumentError validate(restricted, underground)
    nested = EI.Formula(:Ametani2009; underground=EI.Formula(:Bridges1995).routes.mutual)
    @test_throws ArgumentError validate(nested, underground)
    unrestricted = EI.Formula(:Ametani2009; underground=custom)
    @test validate(unrestricted, vertical) === unrestricted
    @test calls[] == 0

    # An unregistered native leaf has no silent validation fallback.
    missing = FM(Val(:MissingValidation), EI.earth_impedance, Val(:mutual))
    @test_throws ArgumentError validate(underground, missing, supplied)
end

@testitem "Engine / earth formulas / physical layers and material vectors" tags=[:unit] begin
    const EN = LineCableModels.Engine
    const EI = EN.EarthImpedance
    const EA = EN.EarthAdmittance
    const EP = LineCableModels.EarthProps
    halfspace = homogeneous(rho=100.0)
    layered = build(EP.EarthModel, (
        EP.EarthLayer(100.0, 10.0, 1.0, 4.0),
        EP.EarthLayer(200.0, 20.0, 1.0)))
    three_soils = build(EP.EarthModel, (
        EP.EarthLayer(100.0, 10.0, 1.0, 4.0),
        EP.EarthLayer(200.0, 20.0, 1.0, 6.0),
        EP.EarthLayer(300.0, 30.0, 1.0)))
    vertical = build(EP.EarthModel, (
        EP.EarthLayer(100.0, 10.0, 1.0),
        EP.EarthLayer(200.0, 20.0, 1.0)); vertical_layers=true)
    for (owner, ids) in ((EI, (:Ametani1974, :Papadopoulos2009, :Papadopoulos2011)),
            (EA, (:Papadopoulos2009, :Papadopoulos2011)))
        for id in ids
            selected = owner.Formula(id)
            @test @inferred(validate(selected, layered)) === selected
            @test_throws DimensionMismatch validate(selected, halfspace)
            @test_throws DimensionMismatch validate(selected, three_soils)
            @test_throws ArgumentError validate(selected, vertical)
        end
    end
    nakagawa = EI.Formula(:Nakagawa1973)
    @test validate(nakagawa, layered) === nakagawa
    @test validate(nakagawa, three_soils) === nakagawa
    @test_throws DimensionMismatch validate(nakagawa, halfspace)
    for id in (:Sunde1968, :Tsiamitros2008)
        selected = EI.Formula(id)
        @test validate(selected, halfspace) === selected
        @test validate(selected, layered) === selected
        @test validate(selected, three_soils) === selected
        @test_throws ArgumentError validate(selected, vertical)
    end
    # EHEM can supply a homogeneous formula from multiple physical soils.
    xue = EI.Formula(:Xue2018)
    @test validate(xue, layered) === xue
    @test validate(xue, vertical) === xue
    @test_throws DimensionMismatch validate(xue, [100.0])
    @test_throws DimensionMismatch validate(xue, [Inf, 100.0], [1.0], [1.0, 1.0], nothing)
    @test_throws DimensionMismatch validate(xue, [Inf, 100.0], [1.0, 1.0], [1.0, 1.0], [Inf])
    @test_throws DimensionMismatch validate(EI.Formula(:Sunde1968), [Inf, 100.0, 200.0], nothing)
    @test validate(EI.Formula(:Sunde1968), [Inf, 100.0], nothing) isa EI.Formula
    pair = EN.EarthPair(1, 2, (-1.0, -5.0), 1.0, (2, 3))
    @test @inferred(validate(pair, (Inf, 4.0, Inf))) === pair
    @test_throws ArgumentError validate(pair, (Inf, 6.0, Inf))
    @test_throws ArgumentError validate(pair, (Inf, Inf))
end

@testitem "Engine / earth formulas / preflight precedes numerical kernels" tags=[:unit] setup=[TestFixtures] begin
    using .TestFixtures
    const EN = LineCableModels.Engine
    problem = TestFixtures.line_parameters_problem(; frequencies=[50.0])
    calls = Ref(0)
    inner = (args...) -> (calls[] += 1; error("internal impedance should not execute"))
    formulation = Formulation(
        internal_impedance=formula(:Schelkunoff1934; inner),
        earth_impedance=:Bridges1995)
    @test_throws ArgumentError compute(problem, formulation)
    @test calls[] == 0
end
