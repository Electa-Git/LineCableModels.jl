@testitem "EarthProps / EHEM / discovery and reconstruction policy" tags=[:unit] begin
    const EP=LineCableModels.EarthProps
    const EH=EP.EHEM

    expected=(:MartinsBritto2020, :Xue2021)
    @test EH.formulas() == expected
    @test EH.DEFAULT === :Layer
    @test EH.Formula(:default) == EH.Layer(-1)
    @test EH.AfterFD(:default) == EH.AfterFD(EH.Layer(-1))
    @test EH.BeforeFD(:default) == EH.BeforeFD(EH.Layer(-1))
    files=sort(filter(endswith(".jl"), readdir(joinpath(
        pkgdir(LineCableModels), "src", "earthprops", "ehem", "formulas"
    ))))
    @test files == collect(lowercase.(string.(expected)) .* ".jl")

    for identifier in expected
        formula=EH.Formula(identifier)
        @test formula_id(formula) === identifier
        @test EH.assumptions(formula) == EH.assumptions(Val(identifier))
        @test !isempty(description(formula))
        @test isconcretetype(typeof(formula))
        @test EH.rule(EH.AfterFD(formula)) === formula
        @test EH.rule(EH.BeforeFD(formula)) === formula
    end
    @test_throws ArgumentError EH.Formula(:Unknown)
    @test_throws ArgumentError EH.Formula(:Xue2021; unknown = 1)
    @test_throws DomainError EH.Layer(0)
    @test_throws MethodError EH.Layer()

    route=(layout, rho, eps_r, mu_r, model, pair, frequency, values)->
        EP.EarthMaterial(rho[end] / values.scale, eps_r[end], mu_r[end])
    experimental=EH.Formula(:Experiment, route, (scale = 2.0,))
    @test_throws ArgumentError EH.Formula(:default; route)
    @test_throws ArgumentError EH.Formula(:default, route)

    model=build(EP.EarthModel, (
        EP.EarthLayer(100.0, 10.0, 1.0, 5.0),
        EP.EarthLayer(50.0, 5.0, 1.0)
    ))
    rho=collect(getfield.(model.layers, :rho))
    eps_r=collect(getfield.(model.layers, :eps_r))
    mu_r=collect(getfield.(model.layers, :mu_r))
    pair=LineCableModels.Engine.EarthPair(1, 1, (1.0, 1.0), 0.0, (1, 1))

    bottom=@inferred EH.Layer(-1)(Val(:overhead), rho, eps_r, mu_r, model, pair, 50.0)
    top=@inferred EH.Layer(2)(Val(:overhead), rho, eps_r, mu_r, model, pair, 50.0)
    @test bottom == EP.EarthMaterial(model.layers[3])
    @test top == EP.EarthMaterial(model.layers[2])
    custom=@inferred experimental(
        Val(:overhead), rho, eps_r, mu_r, model, pair, 50.0
    )
    @test custom.rho == model.layers[3].rho / 2
    @test_throws BoundsError EH.Layer(4)(
        Val(:overhead), rho, eps_r, mu_r, model, pair, 50.0
    )
end

@testitem "EarthProps / EHEM / legacy recurrence references and precision" tags=[:unit] begin
    using TOML
    const EP=LineCableModels.EarthProps
    const EH=EP.EHEM

    reference=TOML.parsefile(joinpath(
        pkgdir(LineCableModels), "test", "fixtures", "reference", "ehem.toml"
    ))
    rho_soil=reference["rho_ohm_m"]
    eps_soil=reference["eps_r"]
    mu_soil=reference["mu_r"]
    thickness=reference["thickness_m"]
    model=build(EP.EarthModel, (
        EP.EarthLayer(rho_soil[1], eps_soil[1], mu_soil[1], thickness[1]),
        EP.EarthLayer(rho_soil[2], eps_soil[2], mu_soil[2], thickness[2]),
        EP.EarthLayer(rho_soil[3], eps_soil[3], mu_soil[3])
    ))
    rho=collect(getfield.(model.layers, :rho))
    eps_r=collect(getfield.(model.layers, :eps_r))
    mu_r=collect(getfield.(model.layers, :mu_r))
    pair=LineCableModels.Engine.EarthPair(1, 1, (10.0, 10.0), 0.0, (1, 1))

    for identifier in EH.formulas()
        formula=EH.Formula(identifier)
        expected=reference["formula"][String(identifier)]
        for (index, frequency) in pairs(reference["frequencies_hz"])
            material=@inferred formula(
                Val(:overhead), rho, eps_r, mu_r, model, pair, frequency
            )
            @test material.rho ≈ expected["rho_ohm_m"][index] rtol=2e-13
            @test material.eps_r ≈ expected["eps_r"][index] rtol=2e-13
            @test material.mu_r === 1.0
        end
        @test_throws ArgumentError formula(
            Val(:underground), rho, eps_r, mu_r, model, pair, 50.0
        )
        @test_throws ArgumentError formula(
            Val(:mixed), rho, eps_r, mu_r, model, pair, 50.0
        )
    end

    martins=EH.Formula(:MartinsBritto2020)(
        Val(:overhead), rho, eps_r, mu_r, model, pair, 50.0
    )
    xue=EH.Formula(:Xue2021)(
        Val(:overhead), rho, eps_r, mu_r, model, pair, 50.0
    )
    @test xue.rho ≈ martins.rho rtol=2e-3

    model32=convert(EP.EarthModel{Float32}, model)
    rho32=collect(getfield.(model32.layers, :rho))
    eps32=collect(getfield.(model32.layers, :eps_r))
    mu32=collect(getfield.(model32.layers, :mu_r))
    pair32=LineCableModels.Engine.EarthPair(
        1, 1, (10.0f0, 10.0f0), 0.0f0, (1, 1)
    )
    for identifier in EH.formulas()
        material=@inferred EH.Formula(identifier)(
            Val(:overhead), rho32, eps32, mu32, model32, pair32, 1000.0f0
        )
        @test material isa EP.EarthMaterial{Float32}
    end
end

@testitem "EarthProps / EHEM / FD composition orders remain independent" tags=[:unit] setup=[
    UseEngineSupport,
] begin
    const EP=LineCableModels.EarthProps
    const EH=EP.EHEM
    const EN=LineCableModels.Engine

    model=build(EarthModel, (
        EP.EarthLayer(100.0, 10.0, 1.0, 5.0),
        EP.EarthLayer(500.0, 20.0, 1.0, 10.0),
        EP.EarthLayer(50.0, 5.0, 1.0)
    ))
    frequency=[1.0e6]
    pair=EN.EarthPair(1, 1, (10.0, 10.0), 0.0, (1, 1))

    after=Formulation(
        earth_properties = :CIGRE2019,
        equivalent_earth = EH.AfterFD(:Xue2021)
    )
    before=Formulation(
        earth_properties = :CIGRE2019,
        equivalent_earth = EH.BeforeFD(:Xue2021)
    )
    after_data=EN._earth_data(after, (earth = model, freq = frequency))
    before_data=EN._earth_data(before, (earth = model, freq = frequency))
    after_media=(
        rho = Matrix{Float64}(undef, 2, 1),
        epsilon = Matrix{Float64}(undef, 2, 1),
        mu = Matrix{Float64}(undef, 2, 1)
    )
    before_media=map(similar, after_media)
    @inferred EN.homogenize!(
        after_media, after.methods.equivalent_earth,
        after.methods.earth_properties, after_data, model, [pair], frequency[1], 1
    )
    EN.homogenize!(
        before_media, before.methods.equivalent_earth,
        before.methods.earth_properties, before_data, model, [pair], frequency[1], 1
    )
    @test after_media.rho[2, 1] != before_media.rho[2, 1]
    @test after_media.epsilon[2, 1] != before_media.epsilon[2, 1]
    after_allocation=@allocated EN.homogenize!(
        after_media, after.methods.equivalent_earth,
        after.methods.earth_properties, after_data, model, [pair], frequency[1], 1
    )
    @test after_allocation <= 1024

    layer_after=Formulation(
        earth_properties = :CIGRE2019,
        equivalent_earth = EH.AfterFD(EH.Layer(-1))
    )
    layer_before=Formulation(
        earth_properties = :CIGRE2019,
        equivalent_earth = EH.BeforeFD(EH.Layer(-1))
    )
    layer_after_data=EN._earth_data(layer_after, (earth = model, freq = frequency))
    layer_before_data=EN._earth_data(layer_before, (earth = model, freq = frequency))
    EN.homogenize!(
        after_media, layer_after.methods.equivalent_earth,
        layer_after.methods.earth_properties, layer_after_data,
        model, [pair], frequency[1], 1
    )
    EN.homogenize!(
        before_media, layer_before.methods.equivalent_earth,
        layer_before.methods.earth_properties, layer_before_data,
        model, [pair], frequency[1], 1
    )
    @test after_media.rho == before_media.rho
    @test after_media.epsilon == before_media.epsilon
    @test after_media.mu == before_media.mu
end
