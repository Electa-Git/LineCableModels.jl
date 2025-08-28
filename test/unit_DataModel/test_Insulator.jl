@testitem "DataModel(Insulator): constructor unit tests" setup = [defaults, deps_datamodel, defs_materials] begin

    using Measurements

    @testset "Input Validation" begin
        # Missing required arguments
        @test_throws ArgumentError Insulator()
        @test_throws ArgumentError Insulator(radius_in=0.01)
        @test_throws ArgumentError Insulator(radius_in=0.01, radius_ext=0.015)

        # Invalid types
        @test_throws ArgumentError Insulator("foo", 0.015, insulator_props, temperature=20.0)
        @test_throws ArgumentError Insulator(0.01, "bar", insulator_props, temperature=20.0)
        @test_throws ArgumentError Insulator(0.01, 0.015, "not_a_material", temperature=20.0)
        @test_throws ArgumentError Insulator(0.01, 0.015, insulator_props, temperature="not_a_temp")

        # Out-of-range values
        @test_throws ArgumentError Insulator(-0.01, 0.015, insulator_props, temperature=20.0)
        @test_throws ArgumentError Insulator(0.01, -0.015, insulator_props, temperature=20.0)

        # Geometrically impossible values
        @test_throws ArgumentError Insulator(0.015, 0.01, insulator_props, temperature=20.0)
        @test_throws ArgumentError Insulator(0.01, 0.01, insulator_props, temperature=20.0)

        # Invalid nothing/missing
        @test_throws ArgumentError Insulator(nothing, 0.015, insulator_props, temperature=20.0)
        @test_throws ArgumentError Insulator(0.01, nothing, insulator_props, temperature=20.0)
        @test_throws ArgumentError Insulator(0.01, 0.015, nothing, temperature=20.0)
        @test_throws ArgumentError Insulator(0.01, 0.015, insulator_props, temperature=nothing)
    end

    @testset "Basic Functionality" begin
        i = Insulator(0.01, 0.015, insulator_props, temperature=20.0)
        @test i isa Insulator
        @test i.radius_in ≈ 0.01 atol = TEST_TOL
        @test i.radius_ext ≈ 0.015 atol = TEST_TOL
        @test i.material_props === insulator_props
        @test i.temperature ≈ 20.0 atol = TEST_TOL
        @test i.cross_section ≈ π * (0.015^2 - 0.01^2) atol = TEST_TOL
        # Measurement type
        i2 = Insulator(measurement(0.01, 1e-5), measurement(0.015, 1e-5), insulator_props, temperature=measurement(20.0, 0.1))
        @test i2 isa Insulator
        @test value(i2.radius_in) ≈ 0.01 atol = TEST_TOL
        @test value(i2.radius_ext) ≈ 0.015 atol = TEST_TOL
        @test value(i2.temperature) ≈ 20.0 atol = TEST_TOL
    end

    @testset "Edge Cases" begin
        # radius_in very close to radius_ext
        i = Insulator(1e-6, 1.0001e-6, insulator_props, temperature=20.0)
        @test i.radius_in ≈ 1e-6 atol = TEST_TOL
    end

    @testset "Physical Behavior" begin
        # Cross-section increases with radius_ext
        i_small = Insulator(0.01, 0.012, insulator_props, temperature=20.0)
        i_large = Insulator(0.01, 0.018, insulator_props, temperature=20.0)
        @test i_large.cross_section > i_small.cross_section
        # but the capacitance decreases (2 pi eps / log(rex/rin))
        @test i_large.shunt_capacitance < i_small.shunt_capacitance
    end

    @testset "Type Stability & Promotion" begin
        # All Float64
        i = Insulator(0.01, 0.015, insulator_props, temperature=20.0)
        @test typeof(i.radius_in) == Float64
        # All Measurement
        iM = Insulator(measurement(0.01, 1e-5), measurement(0.015, 1e-5), insulator_props, temperature=measurement(20.0, 0.1))
        @test typeof(iM.radius_in) <: Measurement
        # Mixed: radius_in as Measurement
        iMix1 = Insulator(measurement(0.01, 1e-5), 0.015, insulator_props, temperature=20.0)
        @test typeof(iMix1.radius_in) <: Measurement
        # Mixed: temperature as Measurement
        iMix2 = Insulator(0.01, 0.015, insulator_props, temperature=measurement(20.0, 0.1))
        @test typeof(iMix2.temperature) <: Measurement
        mmat = Material(1e14, measurement(5.0, 0.1), 1.0, 20.0, 0.0)
        iMix3 = Insulator(0.01, 0.015, mmat, temperature=20.0)
        @test typeof(iMix3.shunt_conductance) <: Measurement
    end

end
