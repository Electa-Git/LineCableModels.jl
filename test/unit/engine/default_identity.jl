@testitem "Engine / retained defaults and PSCAD comparison inventory" tags=[:unit] begin
    @test (@inferred Formulation()) isa LineParametersFormulation
    const E=LineCableModels.Engine
    const EP=LineCableModels.Earth
    expected = (
        (E.InternalImpedance, (:default,)),
        (E.InsulationImpedance, (:default,)),
        (E.EarthImpedance,
            (:default, :Carson1926, :Pollaczek1926, :Gary1976,
                :WedepohlWilcox1973, :Saad1996, :Ametani2009, :Lucca1994)),
        (E.InsulationAdmittance, (:default, :Ametani2004)),
        (E.SemiconAdmittance, (:default, :Ametani2004)),
        (E.EarthAdmittance, (:default, :Pollaczek1926, :IdealGround)),
        (E.PipeImpedance, (:default,)), (EP.FrequencyDependent, (:default,)), (
            EP.EquivalentHomogeneous, (:default,)),
        (LineCableModels.Transforms, (:default,)))
    for (owner, ids) in expected
        @test Set(owner.formulas()) == Set(ids)
        for id in ids
            selected=owner.Formula(id)
            @test formula_id(selected) === id
            @test isconcretetype(typeof(selected))
        end
        @test_throws ArgumentError owner.Formula(:RemovedAuthor)
    end
    @test Formulation().methods.earth_impedance.equivalent_earth === nothing
    @test formula_id(Formulation(earth_impedance = formula(:default; equivalent_earth = formula(:default))).methods.earth_impedance.equivalent_earth.rule) ===
          :default
end
