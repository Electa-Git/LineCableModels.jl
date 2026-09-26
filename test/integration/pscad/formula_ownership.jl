@testitem "PSCAD / Engine owns scientific selections and explicit defaults" tags=[:integration] begin
    const E = LineCableModels.Engine
    const P = LineCableModels.PSCAD
    const FM = LineCableModels.FormulaMethod
    selected = Formulation(:pscad)
    @test !isdefined(P, :NativeFormula)
    @test selected.methods.internal_impedance isa E.InternalImpedance.Formula{:wedepohl1973}
    @test selected.methods.insulation_impedance isa E.InsulationImpedance.Formula{:ametani1980}
    @test selected.methods.earth_admittance isa E.EarthAdmittance.Formula{:ideal}
    choices = (air=:carson1926, earth=:pollaczek1926, mixed=:lucca1994)
    @test map(formula_id, selected.methods.earth_impedance) == choices
    @test all(value -> value isa E.EarthImpedance.Formula, selected.methods.earth_impedance)
    @test map(formula_id, Formulation(:pscad;
        earth_impedance=(air=:default, earth=:default, mixed=:default)).methods.earth_impedance) == choices
    @test keys(Formulation(:pscad; earth_impedance=(mixed=:default,)).methods.earth_impedance) == (:mixed,)
    @test Formulation(:pscad; earth_impedance=(mixed=nothing,)).methods.earth_impedance.mixed === nothing
    for (slot, obsolete) in ((:earth_impedance,:direct_lucca), (:earth_admittance,:coupled),
            (:internal_impedance,:cable_coax), (:insulation_impedance,:cable_coax))
        @test_throws ArgumentError Formulation(:pscad; NamedTuple{(slot,)}((obsolete,))...)
    end
    for owner in (E.InternalImpedance, E.InsulationImpedance, E.EarthImpedance, E.EarthAdmittance)
        @test isempty(intersect(owner.formulas(), (:direct_lucca, :coupled, :cable_coax)))
    end
    record = NamedTuple(selected)
    @test map(value -> value.identifier, record.methods.earth_impedance) == choices
    @test record.methods.earth_admittance.identifier === :ideal
    label = description(selected, Z)
    @test all(occursin(text, label) for text in ("air", "earth", "mixed", "Carson", "Pollaczek", "Lucca", "Wedepohl"))
    @test !any(occursin(text, label) for text in ("Direct/Lucca", "Cable_Coax", "not yet implemented"))
    @test occursin("Ideal", description(selected, Y))
    for kind in (:inner,:outer,:transfer)
        binding = FM(selected.methods.internal_impedance, E.InternalImpedance.internal_impedance, Val(kind))
        @test binding(Val(:pscad)) == (;)
        @test_throws r"not yet implemented" binding(nothing, nothing)
    end
    @test_throws r"not yet implemented" E.InternalImpedance.surface_impedances(
        selected.methods.internal_impedance, 0.003, 0.004, 1.72e-8, 1.0, 100pi*im)
    @test_throws ArgumentError Formulation(:pscad;
        internal_impedance=formula(:wedepohl1973; options=(integration=(method=:quad,),)))
    @test_throws ArgumentError Formulation(:pscad;
        earth_admittance=formula(:ideal; equivalent_earth=formula(:default)))
    @test_throws ArgumentError Formulation(:pscad;
        earth_admittance=formula(:ideal; options=(integration=(method=:quad,),)))
end

@testitem "PSCAD / ideal potential registration and native selection limits" tags=[:integration] begin
    const E = LineCableModels.Engine
    const P = LineCableModels.PSCAD
    metal = Material(:conductor, 1.72e-8, 1.0)
    dielectric = Material(:insulator, 1e14, 2.3)
    design = build(CableDesign, "ideal-potential", terminal(:core,
        core(metal; r=0.004), insulation(dielectric; t=0.002)))
    system = build(LineCableSystem, fill(design, 4),
        [Pose2(0,2), Pose2(1,3), Pose2(2,-1), Pose2(3,-2)];
        connections=[Dict(:core=>i) for i in 1:4])
    problem = LineParametersProblem(system; earth_props=homogeneous(rho=100.0), frequencies=[50.0])
    expected = Set((:ideal,kind,s,t) for (kind,s,t) in
        ((:self,1,1),(:mutual,1,1),(:self,2,2),(:mutual,2,2),(:mutual,1,2),(:mutual,2,1)))
    for choice in ((air=:carson1926,earth=:pollaczek1926,mixed=:lucca1994),
            (air=:gary1976,earth=:wedepohl1973,mixed=:ametani2009),
            (air=:gary1976,earth=:saad1996,mixed=:lucca1994))
        selected = Formulation(:pscad; earth_impedance=choice, earth_admittance=:ideal)
        settings = P.pscad_setting(selected, problem)
        @test Set((row.formula,row.kind,row.source,row.target) for row in settings.interactions.earth_admittance) == expected
        for row in settings.interactions.earth_admittance
            @test P.earth_potential_coefficient(selected.methods.earth_admittance,
                Val(row.kind),Val(row.source),Val(row.target),Val(:pscad)) == (;)
            @test_throws r"not yet implemented" E.EarthAdmittance.earth_potential_coefficient(
                selected.methods.earth_admittance,Val(row.kind),Val(row.source),Val(row.target),nothing,nothing,nothing)
        end
    end
    ideal = E.EarthAdmittance.Formula(:ideal)
    @test_throws ArgumentError P.earth_potential_coefficient(ideal,Val(:self),Val(1),Val(2),Val(:pscad))
    @test_throws ArgumentError P.earth_potential_coefficient(ideal,Val(:mutual),Val(2),Val(3),Val(:pscad))
    @test_throws ArgumentError P.pscad_setting(Formulation(:pscad; earth_admittance=:pollaczek1926),problem)
    @test_throws ArgumentError P.pscad_setting(Formulation(:pscad; internal_impedance=:schelkunoff1934),problem)
    for id in (:ametani2009,:lucca1994), (kind,s,t) in ((:self,1,1),(:mutual,1,1),(:self,2,2),(:mutual,2,2),(:self,1,2))
        @test_throws ArgumentError P.earth_impedance(E.EarthImpedance.Formula(id),Val(kind),Val(s),Val(t),Val(:pscad))
    end
end
