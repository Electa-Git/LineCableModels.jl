@testitem "Earth / static and user-owned frequency-dependent relations" tags=[:unit] setup=[FormulaContractModels] begin
    using Measurements: measurement,uncertainty
    const EP=LineCableModels.Earth
    relation=EP.FrequencyDependent.Formula(:default)
    for T in (Float32,Float64,BigFloat)
        material=EP.EarthMaterial(T(100),T(10),T(1))
        @test relation(material,T(50)) === material
    end
    material=EP.EarthMaterial(measurement(100.0,2.0),measurement(10.0,0.2),measurement(1.0,0.01))
    @test uncertainty(relation(material,50.0).rho)==2.0
    selected=FormulaContractModels.DispersiveEarth()
    @test EP.FrequencyDependent.Formula(selected) === selected
    @test selected(EP.EarthMaterial(100.0,10.0,1.0),100.0).rho==50.0
    @test only(selected.seen)[2]==100.0
    @test_throws DomainError selected(EP.EarthMaterial(100.0,10.0,1.0),0.0)
    @test_throws ArgumentError EP.FrequencyDependent.Formula(:default;parameters=(unknown=1,))
end

@testitem "Earth / soil coefficients are checked at construction" tags=[:unit] begin
    const FD=LineCableModels.Earth.FrequencyDependent
    for id in FD.formulas()
        selected=FD.Formula(id)
        for name in keys(selected.parameters), value in ("bad",true,NaN,Inf,1+im)
            @test_throws ArgumentError FD.Formula(id;parameters=NamedTuple{(name,)}((value,)))
        end
    end
    for (id,name,value) in ((:longmire1975,:corner_scale,0.0),
            (:visacro1987,:normalization_frequency,-1.0),
            (:visacro2012,:frequency_boundary,0.0),
            (:messier1985,:epsilon_infinity,-1.0),
            (:portela1999,:exponent,1.0),(:alipio2014,:exponent,-1.0))
        @test_throws ArgumentError FD.Formula(id;parameters=NamedTuple{(name,)}((value,)))
    end
end
