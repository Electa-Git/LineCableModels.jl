@testitem "Engine / retained defaults and PSCAD comparison inventory" tags=[:unit] begin
    @test Formulation() isa LineParametersFormulation
    const E=LineCableModels.Engine
    const EP=LineCableModels.Earth
    expected = (
        (E.InternalImpedance, (:default, :schelkunoff1934), :schelkunoff1934),
        (E.InsulationImpedance, (:default, :ametani1980), :ametani1980),
        (E.EarthImpedance,
            (:default, :unified, :carson1926, :pollaczek1926, :gary1976,
                :wedepohl1973, :saad1996, :ametani2009, :lucca1994, :wise1934, :xue2018), :unified),
        (E.InsulationAdmittance, (:default, :lossless, :lossy), :lossless),
        (E.SemiconAdmittance, (:default, :lossless, :lossy), :lossless),
        (E.EarthAdmittance, (:default, :unified, :pollaczek1926, :wise1948, :xue2018), :unified),
        (E.PipeImpedance, (:default, :none), :none), (EP.FrequencyDependent,
            (:default, :constant, :alipio2014, :cigre2019, :datsios2019,
                :longmire1975, :messier1985, :portela1999, :scott1967,
                :visacro1987, :visacro2012), :constant), (
            EP.EquivalentHomogeneous, (:default, :bottommost), :bottommost),
        (LineCableModels.Transforms, (:default, :chrysochos2014), :chrysochos2014),
        (LineCableModels.Materials.TemperatureDependent, (:default, :linear), :linear),
        (E.ShuntModel, (:default, :coaxial, :boundary), :coaxial))
    for (owner, ids, target) in expected
        @test Set(owner.formulas()) == Set(ids)
        for id in ids
            selected=owner.Formula(id)
            @test formula_id(selected) === (id === :default ? target : id)
            @test isconcretetype(typeof(selected))
        end
        @test_throws ArgumentError owner.Formula(:RemovedAuthor)
    end
    @test Formulation().methods.earth_impedance.equivalent_earth === nothing
    @test_throws ArgumentError Formulation(earth_admittance = :IdealGround)
    @test formula_id(Formulation(earth_impedance = formula(:default; equivalent_earth = formula(:default))).methods.earth_impedance.equivalent_earth.rule) ===
          :bottommost
end
