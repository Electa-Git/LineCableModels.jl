@testitem "DataModel(Tubular): constructor unit tests" setup = [defaults, deps_datamodel, defs_materials] begin
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
        @test_throws ArgumentError Tubular(0.01, 0.02, material, temperature="25")
        @test_throws ArgumentError Tubular(-0.01, 0.02, material)
        @test_throws ArgumentError Tubular(0.01, -0.02, material)
        @test_throws ArgumentError Tubular(0.03, 0.02, material)
        # Invalid nothing/missing
        @test_throws ArgumentError Tubular(nothing, 0.02, material)
        @test_throws ArgumentError Tubular(0.01, nothing, material)
        @test_throws ArgumentError Tubular(0.01, 0.02, nothing)
        @test_throws ArgumentError Tubular(missing, 0.02, material)
        @test_throws ArgumentError Tubular(0.01, missing, material)
        @test_throws ArgumentError Tubular(0.01, 0.02, material, temperature=missing)
    end

    # Basic Functionality
    @testset "Basic Functionality" begin
        material = Material(1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
        t = Tubular(0.01, 0.02, material)
        @test t isa Tubular
        @test isapprox(t.radius_in, 0.01, atol=TEST_TOL)
        @test isapprox(t.radius_ext, 0.02, atol=TEST_TOL)
        @test t.material_props == material
        @test isapprox(t.temperature, 20.0, atol=TEST_TOL)
        @test isapprox(t.cross_section, π * (0.02^2 - 0.01^2), atol=TEST_TOL)
        t2 = Tubular(t, Thickness(0.02), material)
        @test t2 isa Tubular
        @test isapprox(t2.radius_in, t.radius_ext, atol=TEST_TOL)
    end

    # Edge Cases
    @testset "Edge Cases" begin
        material = Material(1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
        # Very small but positive thickness
        eps = 1e-12
        t = Tubular(0.01, 0.01 + eps, material)
        @test t.radius_ext > t.radius_in
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
        @test t1.radius_in isa Float64
        # All Measurement
        mmat = Material(measurement(1.7241e-8, 1e-10), 1.0, 1.0, 20.0, 0.00393)
        t2 = Tubular(0.011, 0.021, mmat)
        @test t2.radius_in isa Measurement
        # Mixed: radius_in as Measurement
        t3 = Tubular(m, 0.02, material)
        @test t3.radius_in isa Measurement
        # Mixed: radius_ext as Measurement
        t4 = Tubular(0.001, m, material)
        @test t4.radius_ext isa Measurement
        # Mixed: material_props as Measurement
        t5 = Tubular(0.01, 0.02, mmat)
        @test t5.material_props.rho isa Measurement
    end
end
