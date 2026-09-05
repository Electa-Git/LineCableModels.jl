@testsnippet current_earthprops begin
    const EP = LineCableModels.EarthProps
end

@testitem "EarthProps / ephemeral material validation and promotion" tags=[:unit] setup=[
    current_earthprops,
] begin
    using Measurements: Measurement, measurement, uncertainty, value

    material=EP.EarthMaterial(100, 10.0f0, 1)
    @test material isa EP.EarthMaterial{Float64}
    @test material isa LineCableModels.AbstractMaterial
    @test (material.rho, material.eps_r, material.mu_r) == (100.0, 10.0, 1.0)
    @test convert(EP.EarthMaterial{Float64}, material) === material
    @test EP.EarthMaterial(100.0, -1.0, 1.0).eps_r == -1.0

    uncertain=EP.EarthMaterial(
        measurement(100.0, 2.0),
        measurement(10.0, 0.5),
        measurement(1.0, 0.01)
    )
    @test uncertain isa EP.EarthMaterial{Measurement{Float64}}
    @test value(uncertain.rho) == 100.0
    @test uncertainty(uncertain.eps_r) == 0.5

    @test_throws DomainError EP.EarthMaterial(0.0, 10.0, 1.0)
    @test_throws DomainError EP.EarthMaterial(NaN, 10.0, 1.0)
    @test_throws DomainError EP.EarthMaterial(100.0, 0.0, 1.0)
    @test_throws DomainError EP.EarthMaterial(100.0, 10.0, Inf)
end

@testitem "EarthProps / static layer validation and promotion" tags=[:unit] setup=[
    current_earthprops,
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

@testitem "EarthProps / immutable static model ownership" tags=[:unit] setup=[
    current_earthprops,
] begin
    using Measurements: Measurement, measurement, uncertainty, value

    model=@inferred build(EP.EarthModel,
        (
            EP.EarthLayer(100.0, 10.0, 1.0, 20.0),
            EP.EarthLayer(200.0, 15.0, 1.0, 40.0),
            EP.EarthLayer(500.0, 20.0, 1.0)
        ))
    @test model isa EP.EarthModel{Float64}
    @test !ismutabletype(typeof(model))
    @test !model.vertical_layers
    @test model.layers isa NTuple{4, EP.EarthLayer{Float64}}
    @test length(model.layers) == 4
    @test isinf(first(model.layers).rho)
    @test isinf(first(model.layers).thickness)
    @test isinf(last(model.layers).thickness)
    @test !hasfield(typeof(model), :frequencies)
    @test !hasfield(typeof(model), :freq_dependence)
    @test getproperty.(model.layers, :rho) == (Inf, 100.0, 200.0, 500.0)
    @test !applicable(add!, model, EP.EarthLayer(800.0, 25.0, 1.0))
    @test_throws MethodError setindex!(
        model.layers, EP.EarthLayer(800.0, 25.0, 1.0), 2
    )
    @test_throws ArgumentError build(
        EP.EarthModel, (
            EP.EarthLayer(100.0, 10.0, 1.0),
            EP.EarthLayer(800.0, 25.0, 1.0, 10.0)
        ))

    converted=@inferred convert(
        EP.EarthModel{Float32}, EP.EarthModel(100.0, 10.0, 1.0)
    )
    @test converted isa EP.EarthModel{Float32}
    @test convert(EP.EarthModel{Float32}, converted) === converted

    uncertain=convert(EP.EarthModel{Measurement{Float64}}, model)
    @test uncertain isa EP.EarthModel{Measurement{Float64}}
    @test value(uncertain.layers[2].rho) == 100.0
    @test uncertainty(uncertain.layers[2].rho) == 0.0

    vertical=build(EP.EarthModel,
        (
            EP.EarthLayer(100.0, 10.0, 1.0),
            EP.EarthLayer(200.0, 15.0, 1.0, 20.0),
            EP.EarthLayer(500.0, 20.0, 1.0)
        );
        vertical_layers = true)
    @test vertical.vertical_layers
    @test length(vertical.layers) == 4
    @test_throws ArgumentError EP.EarthModel(
        100.0, 10.0, 1.0; thickness = 20.0, vertical_layers = true
    )
end

@testitem "EarthProps / concise display" tags=[:unit] setup=[current_earthprops] begin
    model=build(EP.EarthModel, (
        EP.EarthLayer(100.0, 10.0, 1.0, 20.0),
        EP.EarthLayer(500.0, 20.0, 1.0)
    ))

    shown=sprint(show, "text/plain", model)
    @test contains(shown, "EarthModel · 2 earth layers")
    @test contains(shown, "air")
    @test contains(shown, "layer 1")
    @test contains(shown, "basement")
    @test contains(shown, "ρ=100 Ω·m")
    @test !contains(shown, "frequency")
end

@testitem "EarthProps / static constitutive pass-through" tags=[:unit] setup=[
    current_earthprops,
] begin
    using Measurements: Measurement, measurement, uncertainty, value

    model=build(EP.EarthModel, (
        EP.EarthLayer(100.0, 10.0, 1.0, 20.0),
        EP.EarthLayer(500.0, 20.0, 2.0)
    ))
    frequencies=[50.0, 60.0, 1.0e3]
    formulation=LineCableModels.Engine.Formulation()
    @test formulation.methods.earth_properties === nothing
    material=EP.EarthMaterial(model.layers[2])
    @test constitutive(nothing, material, first(frequencies)) === material
    properties=LineCableModels.Engine._earth_data(
        formulation,
        (earth = model, freq = frequencies)
    )
    static_rho=repeat(collect(getproperty.(model.layers, :rho)), 1, length(frequencies))
    static_eps_r=repeat(
        collect(getproperty.(model.layers, :eps_r)), 1, length(frequencies)
    )
    static_mu_r=repeat(
        collect(getproperty.(model.layers, :mu_r)), 1, length(frequencies)
    )

    @test size(properties.rho) == (3, 3)
    @test size(properties.eps_r) == (3, 3)
    @test size(properties.mu_r) == (3, 3)
    @test properties.rho[:, 1] == [Inf, 100.0, 500.0]
    @test all(properties.rho[:, column] == properties.rho[:, 1]
    for column in axes(properties.rho, 2))
    @test properties.rho == static_rho
    @test properties.eps_r == static_eps_r
    @test properties.mu_r == static_mu_r

    model32=EP.EarthModel(100.0f0, 10.0f0, 1.0f0)
    properties32=LineCableModels.Engine._earth_data(
        formulation,
        (earth = model32, freq = Float32[50, 60])
    )
    @test eltype(properties32.rho) === Float32
    @test eltype(properties32.eps_r) === Float32
    @test eltype(properties32.mu_r) === Float32

    uncertain_model=EP.EarthModel(
        measurement(100.0, 2.0), measurement(10.0, 0.5), measurement(1.0, 0.01)
    )
    uncertain_frequency=measurement.([50.0, 60.0], 0.0)
    uncertain=LineCableModels.Engine._earth_data(
        formulation,
        (earth = uncertain_model, freq = uncertain_frequency)
    )
    @test eltype(uncertain.rho) === Measurement{Float64}
    @test value(uncertain.rho[2, 1]) == 100.0
    @test uncertainty(uncertain.rho[2, 1]) == 2.0
    @test uncertainty(uncertain.eps_r[2, 1]) > 0
    @test uncertainty(uncertain.mu_r[2, 1]) > 0
end
