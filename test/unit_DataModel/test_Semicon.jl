@testitem "DataModel(Semicon): constructor unit tests" setup = [defaults, deps_datamodel, defs_materials] begin

    using Measurements

    @testset "Input Validation" begin
        # Missing required arguments
        @test_throws ArgumentError Semicon()
        @test_throws ArgumentError Semicon(radius_in=0.01)
        @test_throws ArgumentError Semicon(radius_in=0.01, radius_ext=0.012)

        # Invalid types
        @test_throws ArgumentError Semicon("foo", 0.012, semicon_props, temperature=20.0)
        @test_throws ArgumentError Semicon(0.01, "bar", semicon_props, temperature=20.0)
        @test_throws ArgumentError Semicon(0.01, 0.012, "not_a_material", temperature=20.0)
        @test_throws ArgumentError Semicon(0.01, 0.012, semicon_props, temperature="not_a_temp")

        # Out-of-range values
        @test_throws ArgumentError Semicon(-0.01, 0.012, semicon_props, temperature=20.0)
        @test_throws ArgumentError Semicon(0.01, -0.012, semicon_props, temperature=20.0)

        # Geometrically impossible values
        @test_throws ArgumentError Semicon(0.012, 0.01, semicon_props, temperature=20.0)
        @test_throws ArgumentError Semicon(0.01, 0.01, semicon_props, temperature=20.0)

        # Invalid nothing/missing
        @test_throws ArgumentError Semicon(nothing, 0.012, semicon_props, temperature=20.0)
        @test_throws ArgumentError Semicon(0.01, nothing, semicon_props, temperature=20.0)
        @test_throws ArgumentError Semicon(0.01, 0.012, nothing, temperature=20.0)
        @test_throws ArgumentError Semicon(0.01, 0.012, semicon_props, nothing)
    end

    @testset "Basic Functionality" begin
        s = Semicon(0.01, 0.012, semicon_props, temperature=20.0)
        @test s isa Semicon
        @test s.radius_in ≈ 0.01 atol = TEST_TOL
        @test s.radius_ext ≈ 0.012 atol = TEST_TOL
        @test s.material_props === semicon_props
        @test s.temperature ≈ 20.0 atol = TEST_TOL
        @test s.cross_section ≈ π * (0.012^2 - 0.01^2) atol = TEST_TOL
        # Measurement type
        s2 = Semicon(measurement(0.01, 1e-5), measurement(0.012, 1e-5), semicon_props, temperature=measurement(20.0, 0.1))
        @test s2 isa Semicon
        @test value(s2.radius_in) ≈ 0.01 atol = TEST_TOL
        @test value(s2.radius_ext) ≈ 0.012 atol = TEST_TOL
        @test value(s2.temperature) ≈ 20.0 atol = TEST_TOL
    end

    @testset "Edge Cases" begin
        # radius_in very close to radius_ext
        s = Semicon(1e-6, 1.0001e-6, semicon_props, temperature=20.0)
        @test s.radius_in ≈ 1e-6 atol = TEST_TOL
    end

    @testset "Physical Behavior" begin
        # Cross-section increases with radius_ext
        s_small = Semicon(0.01, 0.011, semicon_props, temperature=20.0)
        s_large = Semicon(0.01, 0.013, semicon_props, temperature=20.0)
        @test s_large.cross_section > s_small.cross_section
    end

    @testset "Type Stability & Promotion" begin
        # All Float64
        s = Semicon(0.01, 0.012, semicon_props, temperature=20.0)
        @test typeof(s.radius_in) == Float64
        # All Measurement
        sM = Semicon(measurement(0.01, 1e-5), measurement(0.012, 1e-5), semicon_props, temperature=measurement(20.0, 0.1))
        @test typeof(sM.radius_in) <: Measurement
        # Mixed: radius_in as Measurement
        sMix1 = Semicon(measurement(0.01, 1e-5), 0.012, semicon_props, temperature=20.0)
        @test typeof(sMix1.radius_in) <: Measurement
        # Mixed: temperature as Measurement
        sMix2 = Semicon(0.01, 0.012, semicon_props, temperature=measurement(20.0, 0.1))
        @test typeof(sMix2.temperature) <: Measurement
        mmat = Material(1000.0, measurement(1000.0, 0.1), 1.0, 20.0, 0.0)
        sMix3 = Semicon(0.01, 0.012, mmat, temperature=20.0)
        @test typeof(sMix3.shunt_capacitance) <: Measurement
    end

end
