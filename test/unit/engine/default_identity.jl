@testitem "Engine / retained defaults and PSCAD comparison inventory" tags=[:unit] begin
    @test Formulation() isa LineParametersFormulation
    const E=LineCableModels.Engine
    const EP=LineCableModels.Earth
    expected = (
        (E.InternalImpedance, (:default, :schelkunoff1934)),
        (E.InsulationImpedance, (:default, :ametani1980)),
        (E.EarthImpedance,
            (:default, :carson1926, :pollaczek1926, :gary1976,
                :wedepohl1973, :saad1996, :ametani2009, :lucca1994, :wise1934, :xue2018)),
        (E.InsulationAdmittance, (:default, :lossless, :lossy)),
        (E.SemiconAdmittance, (:default, :lossless, :lossy)),
        (E.EarthAdmittance, (:default, :pollaczek1926, :wise1948, :xue2018)),
        (E.PipeImpedance, (:default,)), (EP.FrequencyDependent, (:default,)), (
            EP.EquivalentHomogeneous, (:default,)),
        (LineCableModels.Transforms, (:default, :chrysochos2014)))
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
    @test_throws ArgumentError Formulation(earth_admittance = :IdealGround)
    @test formula_id(Formulation(earth_impedance = formula(:default; equivalent_earth = formula(:default))).methods.earth_impedance.equivalent_earth.rule) ===
          :default
end
