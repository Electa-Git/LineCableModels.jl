@testitem "DataModel(Strip): constructor unit tests" setup = [defaults, deps_datamodel, defs_materials] begin

    using Measurements

    @testset "Input Validation" begin
        # Missing required arguments
        @test_throws ArgumentError Strip()
        @test_throws ArgumentError Strip(radius_in=0.01)
        @test_throws ArgumentError Strip(radius_in=0.01, radius_ext=0.012)
        @test_throws ArgumentError Strip(radius_in=0.01, radius_ext=0.012, width=0.05)
        @test_throws ArgumentError Strip(radius_in=0.01, radius_ext=0.012, width=0.05, lay_ratio=10)

        # Invalid types
        @test_throws ArgumentError Strip("foo", 0.012, 0.05, 10, copper_props)
        @test_throws ArgumentError Strip(0.01, "bar", 0.05, 10, copper_props)
        @test_throws ArgumentError Strip(0.01, 0.012, "baz", 10, copper_props)
        @test_throws ArgumentError Strip(0.01, 0.012, 0.05, "qux", copper_props)
        @test_throws ArgumentError Strip(0.01, 0.012, 0.05, 10, "not_a_material")
        @test_throws ArgumentError Strip(0.01, 0.012, 0.05, 10, copper_props, "not_a_temp", 1)
        @test_throws ArgumentError Strip(0.01, 0.012, 0.05, 10, copper_props, temperature=20.0, "not_a_dir")

        # Out-of-range values
        @test_throws ArgumentError Strip(-0.01, 0.012, 0.05, 10, copper_props)
        @test_throws ArgumentError Strip(0.01, -0.012, 0.05, 10, copper_props)
        @test_throws ArgumentError Strip(0.01, 0.012, -0.05, 10, copper_props)
        @test_throws ArgumentError Strip(0.01, 0.012, 0.05, 10, copper_props, temperature=20.0, lay_direction=0)
        @test_throws ArgumentError Strip(0.01, 0.012, 0.05, 10, copper_props, temperature=20.0, lay_direction=2)
        @test_throws ArgumentError Strip(0.01, 0.012, 0.05, 10, copper_props, temperature=20.0, lay_direction=-2)

        # Geometrically impossible values
        @test_throws ArgumentError Strip(0.012, 0.01, 0.05, 10, copper_props)
        @test_throws ArgumentError Strip(0.01, 0.01, 0.05, 10, copper_props)

        # Invalid nothing/missing
        @test_throws ArgumentError Strip(nothing, 0.012, 0.05, 10, copper_props)
        @test_throws ArgumentError Strip(0.01, nothing, 0.05, 10, copper_props)
        @test_throws ArgumentError Strip(0.01, 0.012, nothing, 10, copper_props)
        @test_throws ArgumentError Strip(0.01, 0.012, 0.05, nothing, copper_props)
        @test_throws ArgumentError Strip(0.01, 0.012, 0.05, 10, nothing)
        @test_throws ArgumentError Strip(0.01, 0.012, 0.05, 10, copper_props, temperature=nothing, lay_direction=1)
        @test_throws ArgumentError Strip(0.01, 0.012, 0.05, 10, copper_props, temperature=20.0, lay_direction=nothing)
    end

    @testset "Basic Functionality" begin
        s = Strip(0.01, 0.012, 0.05, 10, copper_props)
        @test s isa Strip
        @test s.radius_in ≈ 0.01 atol = TEST_TOL
        @test s.radius_ext ≈ 0.012 atol = TEST_TOL
        @test s.width ≈ 0.05 atol = TEST_TOL
        @test s.lay_ratio ≈ 10 atol = TEST_TOL
        @test s.material_props === copper_props
        @test s.temperature ≈ 20.0 atol = TEST_TOL
        @test s.lay_direction == 1
        @test s.cross_section ≈ (0.012 - 0.01) * 0.05 atol = TEST_TOL
        # measurement type
        s2 = Strip(measurement(0.01, 1e-5), measurement(0.012, 1e-5), measurement(0.05, 1e-4), 10, copper_props, temperature=measurement(20.0, 0.1))
        @test s2 isa Strip
        @test value(s2.radius_in) ≈ 0.01 atol = TEST_TOL
        @test value(s2.radius_ext) ≈ 0.012 atol = TEST_TOL
        @test value(s2.width) ≈ 0.05 atol = TEST_TOL
        @test value(s2.temperature) ≈ 20.0 atol = TEST_TOL
    end

    @testset "Edge Cases" begin
        # radius_in very close to radius_ext
        s = Strip(1e-6, 1.0001e-6, 0.05, 10, copper_props)
        @test s.radius_in ≈ 1e-6 atol = TEST_TOL
        # width very small
        s2 = Strip(0.01, 0.012, 1e-6, 10, copper_props)
        @test s2.width ≈ 1e-6 atol = TEST_TOL
    end

    @testset "Physical Behavior" begin
        # Resistance should increase with temperature
        s20 = Strip(0.01, 0.012, 0.05, 10, copper_props)
        s80 = Strip(0.01, 0.012, 0.05, 10, copper_props, temperature=80.0)
        @test s80.resistance > s20.resistance
        # Cross-section increases with width
        s_small = Strip(0.01, 0.012, 0.01, 10, copper_props)
        s_large = Strip(0.01, 0.012, 0.1, 10, copper_props)
        @test s_large.cross_section > s_small.cross_section
    end

    @testset "Type Stability & Promotion" begin
        # All Float64
        s = Strip(0.01, 0.012, 0.05, 10.0, copper_props)
        @test typeof(s.radius_in) == Float64
        # All measurement
        sM = Strip(measurement(0.01, 1e-5), measurement(0.012, 1e-5), measurement(0.05, 1e-4), measurement(10.0, 0.1), copper_props, temperature=measurement(20.0, 0.1), lay_direction=1)
        @test typeof(sM.radius_in) <: Measurement
        # Mixed: radius_in as measurement
        sMix1 = Strip(measurement(0.01, 1e-5), 0.012, 0.05, 10.0, copper_props)
        @test typeof(sMix1.radius_in) <: Measurement
        # Mixed: width as measurement
        sMix2 = Strip(0.01, 0.012, measurement(0.05, 1e-4), 10.0, copper_props)
        @test typeof(sMix2.width) <: Measurement
        mmat = Material(measurement(1.7241e-8, 1e-10), 1.0, 1.0, 20.0, 0.00393)
        sMix3 = Strip(0.01, 0.012, 0.05, 10.0, mmat, temperature=20.0, lay_direction=1)
        @test typeof(sMix3.resistance) <: Measurement
    end

end
