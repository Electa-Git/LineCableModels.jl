@testitem "BaseParams / geometry / helix and wire coordinates" tags=[:unit] setup=[
    BaseParamsTestSupport,
    UseBaseParamsSupport
] begin
    mean_diameter, pitch, overlength=helix(0.01, 0.015, 12.0)
    @test mean_diameter ≈ 0.025
    @test pitch ≈ 0.3
    @test overlength > 1
    @test helix(0.01, 0.015, 0.0) == (0.025, 0.0, 1.0)
    @test_throws DomainError helix(-0.01, 0.015, 12.0)
    @test_throws DomainError helix(0.02, 0.015, 12.0)
    @test_throws DomainError helix(0.01, 0.015, -1.0)

    coordinates=wire_coordinates(7, 0.001, 0.01)
    @test length(coordinates) == 7
    @test all(point->hypot(point...) ≈ 0.011, coordinates)
    @test wire_coordinates(1, 0.001, 0.0) == [(0.0, 0.0)]
    @test isempty(wire_coordinates(0, 0.001, 0.01))
    @test_throws DomainError wire_coordinates(7, -0.001, 0.01)

    for T in (Float32, Float64, BigFloat)
        values=helix(T(0.01), T(0.015), T(12))
        points=wire_coordinates(7, T(0.001), T(0.01))
        @test all(value->value isa T, values)
        @test all(point->all(value->value isa T, point), points)
    end

    uncertain=helix(
        measurement(0.01, 1e-4),
        measurement(0.015, 1e-4),
        measurement(12.0, 0.1)
    )
    @test all(value->value isa Measurement{Float64}, uncertain)
    @test uncertainty(last(uncertain)) > 0
end

@testitem "BaseParams / geometry / GMD and equivalent GMR" tags=[:unit] setup=[
    BaseParamsTestSupport,
    UseBaseParamsSupport
] begin
    material=Material(:conductor, 1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
    first_layer=CircStrands(0.01, 0.012, 0.001, 7, 10.0, material)
    second_layer=CircStrands(0.02, 0.022, 0.001, 12, 15.0, material)
    tube=Tubular(0.01, 0.02, material)

    @test gmd(first_layer, tube) ≈ gmd(tube, first_layer)
    @test gmd(tube, tube) ≈ tube.r_ex
    @test gmd(second_layer, tube) > gmd(first_layer, tube)
    @test equivalent_gmr(first_layer, second_layer) > 0
    @test equivalent_gmr(first_layer, second_layer) ≈
          equivalent_gmr(second_layer, first_layer)

    material32=convert(Material{Float32}, material)
    first32=CircStrands(
        Float32(0.01),
        Float32(0.012),
        Float32(0.001),
        7,
        Float32(10),
        material32
    )
    tube32=Tubular(Float32(0.01), Float32(0.02), material32)
    @test gmd(first32, tube32) isa Float32
    @test gmd(first32, tube32) ≈ Float32(gmd(first_layer, tube)) rtol=sqrt(eps(Float32))

    uncertain_radius=measurement(0.001, 1e-5)
    uncertain_layer=CircStrands(
        0.01,
        0.01+2uncertain_radius,
        uncertain_radius,
        7,
        10.0,
        material
    )
    @test gmd(uncertain_layer, tube) isa Measurement{Float64}
    @test equivalent_gmr(uncertain_layer, second_layer) isa Measurement{Float64}
    @test uncertainty(equivalent_gmr(uncertain_layer, second_layer)) > 0
end

@testitem "BaseParams / geometry / GMR kernels and solenoid factor" tags=[:unit] setup=[
    BaseParamsTestSupport,
    UseBaseParamsSupport
] begin
    strand=strand_gmr(0.011, 7, 0.001, 1.0)
    tube=tubular_gmr(0.02, 0.01, 1.0)
    @test strand > 0
    @test tube > 0
    @test strand_gmr(0.011, 7, 0.001, 2.0) < strand
    @test tubular_gmr(0.02, 0.01, 2.0) < tube
    @test tubular_gmr(0.02, 0.0, 1.0) ≈ 0.02*exp(-0.25)
    @test tubular_gmr(0.02, 0.02, 1.0) ≈ 0.02

    @test_throws DomainError strand_gmr(0.011, 0, 0.001, 1.0)
    @test_throws DomainError strand_gmr(0.011, 7, -0.001, 1.0)
    @test_throws DomainError tubular_gmr(0.01, 0.02, 1.0)
    @test_throws DomainError tubular_gmr(0.0, 0.0, 1.0)

    @test solenoid_factor(0.0, 0.01, 0.02) == 1.0
    @test solenoid_factor(10.0, 0.01, 0.02) > 1.0
    @test_throws DomainError solenoid_factor(-1.0, 0.01, 0.02)
    @test_throws DomainError solenoid_factor(1.0, 0.02, 0.01)

    for T in (Float32, Float64, BigFloat)
        @test strand_gmr(T(0.011), 7, T(0.001), one(T)) isa T
        @test tubular_gmr(T(0.02), T(0.01), one(T)) isa T
        @test solenoid_factor(T(10), T(0.01), T(0.02)) isa T
    end

    uncertain=tubular_gmr(
        measurement(0.02, 1e-4),
        measurement(0.01, 1e-4),
        measurement(1.0, 0.01)
    )
    @test uncertain isa Measurement{Float64}
    @test uncertainty(uncertain) > 0
end
