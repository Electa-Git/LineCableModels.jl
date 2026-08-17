@testsnippet defs_earthprops begin
    const EP = LineCableModels.EarthProps
    const EEP = LineCableModels.Engine.EarthProperties
end

@testitem "EarthProps / static layer validation and promotion" tags=[:unit] setup=[
    defs_earthprops,
] begin
    using Measurements: Measurement, measurement, uncertainty, value

    layer=EP.EarthLayer(100, 10.0f0, 1, 20.0)
    @test layer isa EP.EarthLayer{Float64}
    @test (layer.rho, layer.eps_r, layer.mu_r, layer.thickness) ==
          (100.0, 10.0, 1.0, 20.0)

    layer32=EP.EarthLayer(100.0f0, 10.0f0, 1.0f0)
    @test layer32 isa EP.EarthLayer{Float32}
    @test isinf(layer32.thickness)
    @test convert(EP.EarthLayer{Float32}, layer32) === layer32

    layer_measurement=convert(EP.EarthLayer{Measurement{Float64}}, layer)
    @test value(layer_measurement.rho) == 100.0
    @test uncertainty(layer_measurement.rho) == 0.0
    @test value(layer_measurement.thickness) == 20.0

    promoted=EP.EarthLayer(
        measurement(100.0, 2.0), 10.0, 1.0, measurement(20.0, 0.1)
    )
    @test promoted isa EP.EarthLayer{Measurement{Float64}}
    @test uncertainty(promoted.rho) == 2.0
    @test uncertainty(promoted.thickness) == 0.1

    @test_throws DomainError EP.EarthLayer(0.0, 10.0, 1.0)
    @test_throws DomainError EP.EarthLayer(NaN, 10.0, 1.0)
    @test_throws DomainError EP.EarthLayer(100.0, 0.0, 1.0)
    @test_throws DomainError EP.EarthLayer(100.0, 10.0, Inf)
    @test_throws DomainError EP.EarthLayer(100.0, 10.0, 1.0, 0.0)
end

@testitem "EarthProps / static model ownership and atomic mutation" tags=[:unit] setup=[
    defs_earthprops,
] begin
    using Measurements: Measurement, measurement, uncertainty, value

    model=EP.EarthModel(100.0, 10.0, 1.0; thickness = 20.0)
    @test model isa EP.EarthModel{Float64}
    @test !model.vertical_layers
    @test length(model.layers) == 2
    @test isinf(first(model.layers).rho)
    @test isinf(first(model.layers).thickness)
    @test last(model.layers).thickness == 20.0
    @test !hasfield(typeof(model), :frequencies)
    @test !hasfield(typeof(model), :freq_dependence)

    add!(model, EP.EarthLayer(200.0, 15.0, 1.0, 40.0))
    add!(model, EP.EarthLayer(500.0, 20.0, 1.0))
    @test getproperty.(model.layers, :rho) == [Inf, 100.0, 200.0, 500.0]

    snapshot=copy(model.layers)
    @test_throws ArgumentError add!(model, EP.EarthLayer(800.0, 25.0, 1.0, 10.0))
    @test model.layers == snapshot
    @test_throws ArgumentError add!(model, EP.EarthLayer(800.0f0, 25.0f0, 1.0f0))
    @test model.layers == snapshot

    model32=EP.EarthModel(100.0f0, 10.0f0, 1.0f0; thickness = 20.0f0)
    add!(model32, 200.0f0, 15.0f0, 1.0f0; thickness = Float32(Inf))
    @test model32 isa EP.EarthModel{Float32}
    @test_throws ArgumentError add!(model32, 300.0, 15.0, 1.0)

    converted=convert(EP.EarthModel{Float32}, EP.EarthModel(100.0, 10.0, 1.0))
    @test converted isa EP.EarthModel{Float32}
    @test convert(EP.EarthModel{Float32}, converted) === converted

    uncertain=convert(EP.EarthModel{Measurement{Float64}}, model)
    @test uncertain isa EP.EarthModel{Measurement{Float64}}
    @test value(uncertain.layers[2].rho) == 100.0
    @test uncertainty(uncertain.layers[2].rho) == 0.0

    vertical=EP.EarthModel(100.0, 10.0, 1.0; vertical_layers = true)
    add!(vertical, EP.EarthLayer(200.0, 15.0, 1.0, 20.0))
    add!(vertical, EP.EarthLayer(500.0, 20.0, 1.0))
    @test vertical.vertical_layers
    @test length(vertical.layers) == 4
    @test_throws ArgumentError EP.EarthModel(
        100.0, 10.0, 1.0; thickness = 20.0, vertical_layers = true
    )
end

@testitem "EarthProps / tables and concise display" tags=[:unit] setup=[defs_earthprops] begin
    using DataFrames: DataFrame, nrow

    model=EP.EarthModel(100.0, 10.0, 1.0; thickness = 20.0)
    add!(model, EP.EarthLayer(500.0, 20.0, 1.0))

    table=DataFrame(model)
    @test names(table) == ["rho", "eps_r", "mu_r", "thickness"]
    @test nrow(table) == 3
    @test isinf(table.rho[1])
    @test table.rho[2:3] == [100.0, 500.0]

    shown=sprint(show, "text/plain", model)
    @test contains(shown, "2 horizontal earth layers (multilayer)")
    @test contains(shown, "Layer 1 (air)")
    @test contains(shown, "rho=100.0")
    @test !contains(shown, "frequency")
end

@testitem "Engine.EarthProperties / constant evaluation" tags=[:unit] setup=[
    defs_earthprops,
] begin
    using Measurements: Measurement, measurement, uncertainty, value

    formulation=EEP.CPEarth()
    @test LineCableModels.description(formulation) == "Constant properties (CP)"

    model=EP.EarthModel(100.0, 10.0, 1.0; thickness = 20.0)
    add!(model, EP.EarthLayer(500.0, 20.0, 2.0))
    frequencies=[50.0, 60.0, 1.0e3]
    properties=EEP.evaluate(formulation, model, frequencies)

    @test size(properties.rho) == (3, 3)
    @test size(properties.epsilon) == (3, 3)
    @test size(properties.mu) == (3, 3)
    @test properties.rho[:, 1] == [Inf, 100.0, 500.0]
    @test all(properties.rho[:, column] == properties.rho[:, 1]
    for column in axes(properties.rho, 2))
    @test properties.epsilon[1, 1] ≈ 8.8541878128e-12
    @test properties.epsilon[2, 1] ≈ 10 * 8.8541878128e-12
    @test properties.mu[1, 1] ≈ 4π * 1.0e-7
    @test properties.mu[3, 1] ≈ 2 * 4π * 1.0e-7

    model32=EP.EarthModel(100.0f0, 10.0f0, 1.0f0)
    properties32=EEP.evaluate(formulation, model32, Float32[50, 60])
    @test eltype(properties32.rho) === Float32
    @test eltype(properties32.epsilon) === Float32
    @test eltype(properties32.mu) === Float32

    uncertain_model=EP.EarthModel(
        measurement(100.0, 2.0), measurement(10.0, 0.5), measurement(1.0, 0.01)
    )
    uncertain_frequency=measurement.([50.0, 60.0], 0.0)
    uncertain=EEP.evaluate(formulation, uncertain_model, uncertain_frequency)
    @test eltype(uncertain.rho) === Measurement{Float64}
    @test value(uncertain.rho[2, 1]) == 100.0
    @test uncertainty(uncertain.rho[2, 1]) == 2.0
    @test uncertainty(uncertain.epsilon[2, 1]) > 0
    @test uncertainty(uncertain.mu[2, 1]) > 0
end
