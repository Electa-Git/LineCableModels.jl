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

@testitem "Earth / published dispersive relations match numerical references" tags=[:unit] begin
    const EP=LineCableModels.Earth
    const FD=EP.FrequencyDependent
    material=EP.EarthMaterial(100.0,10.0,1.0)
    # Independently evaluated from the published expressions documented by
    # each formula, for rho_0 = 100 ohm m and f = 10 kHz.
    references=(
        alipio2014=(98.085845947624506,409.88833451812962),
        cigre2019=(98.078364240070741,408.02591429681854),
        datsios2019=(99.953087658615345,77.228148054469685),
        longmire1975=(91.26463633360153,548.99524080259334),
        messier1985=(97.1029453119354,544.28505232978239),
        portela1999=(97.229937848620864,1088.2832591726246),
        scott1967=(101.39113857366789,542.00089040162334),
        visacro1987=(71.779429127136169,815.10931718426775),
        visacro2012=(98.649670481737814,192.20336879472805),
    )
    for (identifier,(rho,eps_r)) in pairs(references)
        evaluated=FD.Formula(identifier)(material,1.0e4)
        @test evaluated.rho ≈ rho rtol=2e-14
        @test evaluated.eps_r ≈ eps_r rtol=2e-14
        @test evaluated.mu_r === material.mu_r
    end

    visacro=FD.Formula(:visacro2012)
    @test visacro(material,50.0) == visacro(material,100.0)
    datsios=FD.Formula(:datsios2019)
    @test datsios(material,50.0).eps_r == datsios(material,3000.0).eps_r
end
