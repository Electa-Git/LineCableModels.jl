@testitem "DataModel / Tubular / constructor unit tests" tags=[:unit] setup=[
    DataModelTestSupport, UseDataModelSupport, TestNumerics, TestFixtures, MaterialFixtures] begin
    # Input Validation
    @testset "Input Validation" begin
        material = Material(1.7241e-8, 1.0, 1.0, 20.0, 0.00393)

        # Missing required arguments
        @test_throws ArgumentError Tubular()
        @test_throws ArgumentError Tubular(0.01)
        @test_throws ArgumentError Tubular(0.01, 0.02)
        # Invalid types
        @test_throws ArgumentError Tubular("0.01", 0.02, material)
        @test_throws ArgumentError Tubular(0.01, "0.02", material)
        @test_throws ArgumentError Tubular(0.01, 0.02, "material")
        @test_throws ArgumentError Tubular(0.01, 0.02, material, temperature = "25")
        @test_throws ArgumentError Tubular(-0.01, 0.02, material)
        @test_throws ArgumentError Tubular(0.01, -0.02, material)
        @test_throws ArgumentError Tubular(0.03, 0.02, material)
        # Invalid nothing/missing
        @test_throws ArgumentError Tubular(nothing, 0.02, material)
        @test_throws ArgumentError Tubular(0.01, nothing, material)
        @test_throws ArgumentError Tubular(0.01, 0.02, nothing)
        @test_throws ArgumentError Tubular(missing, 0.02, material)
        @test_throws ArgumentError Tubular(0.01, missing, material)
        @test_throws ArgumentError Tubular(0.01, 0.02, material, temperature = missing)
    end

    # Basic Functionality
    @testset "Basic Functionality" begin
        material = Material(1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
        t = Tubular(0.01, 0.02, material)
        @test t isa Tubular
        @test isapprox(t.r_in, 0.01, atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(t.r_ex, 0.02, atol = TestNumerics.absolute_floor(Float64))
        @test t.material_props == material
        @test isapprox(t.temperature, 20.0, atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(t.cross_section, π * (0.02^2 - 0.01^2), atol = TestNumerics.absolute_floor(Float64))
        t2 = Tubular(t.r_ex, material; thickness = 0.02)
        @test t2 isa Tubular
        @test isapprox(t2.r_in, t.r_ex, atol = TestNumerics.absolute_floor(Float64))
    end

    # Edge Cases
    @testset "Edge Cases" begin
        material = Material(1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
        # Very small but positive thickness
        eps = 1e-12
        t = Tubular(0.01, 0.01 + eps, material)
        @test t.r_ex > t.r_in
        @test t.cross_section > 0
        # Inf radii (should error)
        @test_throws DomainError Tubular(0.01, Inf, material)
    end

    # Physical Behavior
    @testset "Physical Behavior" begin
        material = Material(1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
        t1 = Tubular(0.01, 0.02, material)
        t2 = Tubular(0.01, 0.03, material)
        @test t2.cross_section > t1.cross_section
        @test t2.resistance < t1.resistance
    end

    # Type Stability & Promotion
    @testset "Type Stability & Promotion" begin
        material = Material(1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
        m = measurement(0.01, 0.001)
        # All Float64
        t1 = Tubular(0.01, 0.02, material)
        @test t1.r_in isa Float64
        # All Measurement
        mmat = Material(measurement(1.7241e-8, 1e-10), 1.0, 1.0, 20.0, 0.00393)
        t2 = Tubular(0.011, 0.021, mmat)
        @test t2.r_in isa Measurement
        # Mixed: r_in as Measurement
        t3 = Tubular(m, 0.02, material)
        @test t3.r_in isa Measurement
        # Mixed: r_ex as Measurement
        t4 = Tubular(0.001, m, material)
        @test t4.r_ex isa Measurement
        # Mixed: material_props as Measurement
        t5 = Tubular(0.01, 0.02, mmat)
        @test t5.material_props.rho isa Measurement
    end

    @testset "Strict numeric radii and named radial declarations" begin
        import LineCableModels.DataModel: _normalize_radii
        material = Material(1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
        @test _normalize_radii(Tubular, 0.01, 0.03) == (0.01, 0.03)
        @test_throws ArgumentError _normalize_radii(Tubular, :previous, 0.03)
        @test_throws ArgumentError _normalize_radii(Tubular, 0.01, :diameter)
        @test Tubular(0.01, material; radius = 0.03).r_ex == 0.03
        @test Tubular(0.01, material; thickness = 0.02).r_ex == 0.03
        @test_throws ArgumentError Tubular(0.01, material)
        @test_throws ArgumentError Tubular(
            0.01,
            material;
            radius = 0.03,
            thickness = 0.02
        )
    end
end
