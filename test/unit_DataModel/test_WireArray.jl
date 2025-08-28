@testitem "DataModel(WireArray): constructor unit tests" setup = [defaults, deps_datamodel, defs_materials] begin

    using Measurements

    @testset "Input Validation" begin
        # Missing required arguments
        @test_throws ArgumentError WireArray()
        @test_throws ArgumentError WireArray(radius_in=0.01)
        @test_throws ArgumentError WireArray(radius_in=0.01, radius_wire=0.002)
        @test_throws ArgumentError WireArray(radius_in=0.01, radius_wire=0.002, num_wires=7)
        @test_throws ArgumentError WireArray(radius_in=0.01, radius_wire=0.002, num_wires=7, lay_ratio=10)

        # Invalid types
        @test_throws ArgumentError WireArray("foo", 0.002, 7, 10, copper_props, temperature=20.0, lay_direction=1)
        @test_throws ArgumentError WireArray(0.01, "bar", 7, 10, copper_props, temperature=20.0, lay_direction=1)
        @test_throws ArgumentError WireArray(0.01, 0.002, "baz", 10, copper_props, temperature=20.0, lay_direction=1)
        @test_throws ArgumentError WireArray(0.01, 0.002, 7, "qux", copper_props, temperature=20.0, lay_direction=1)
        @test_throws ArgumentError WireArray(0.01, 0.002, 7, 10, "not_a_material", temperature=20.0, lay_direction=1)
        @test_throws ArgumentError WireArray(0.01, 0.002, 7, 10, copper_props, "not_a_temp", lay_direction=1)
        @test_throws ArgumentError WireArray(0.01, 0.002, 7, 10, copper_props, temperature=20.0, "not_a_dir")

        # Out-of-range values
        @test_throws ArgumentError WireArray(-0.01, 0.002, 7, 10, copper_props, temperature=20.0, lay_direction=1)
        @test_throws ArgumentError WireArray(0.01, -0.002, 7, 10, copper_props, temperature=20.0, lay_direction=1)
        @test_throws ArgumentError WireArray(0.01, 0.002, 0, 10, copper_props, temperature=20.0, lay_direction=1)
        @test_throws ArgumentError WireArray(0.01, 0.002, 7, 10, copper_props, temperature=20.0, lay_direction=0)
        @test_throws ArgumentError WireArray(0.01, 0.002, 7, 10, copper_props, temperature=20.0, lay_direction=2)
        @test_throws ArgumentError WireArray(0.01, 0.002, 7, 10, copper_props, temperature=20.0, lay_direction=-2)

        # Geometrically impossible values
        @test_throws ArgumentError WireArray(0.01, 0.0, 7, 10, copper_props, temperature=20.0, lay_direction=1)

        # Invalid nothing/missing
        @test_throws ArgumentError WireArray(nothing, 0.002, 7, 10, copper_props, temperature=20.0, lay_direction=1)
        @test_throws ArgumentError WireArray(0.01, nothing, 7, 10, copper_props, temperature=20.0, lay_direction=1)
        @test_throws ArgumentError WireArray(0.01, 0.002, nothing, 10, copper_props, temperature=20.0, lay_direction=1)
        @test_throws ArgumentError WireArray(0.01, 0.002, 7, nothing, copper_props, temperature=20.0, lay_direction=1)
        @test_throws ArgumentError WireArray(0.01, 0.002, 7, 10, nothing, temperature=20.0, lay_direction=1)
        @test_throws ArgumentError WireArray(0.01, 0.002, 7, 10, copper_props, temperature=nothing, lay_direction=1)
        @test_throws ArgumentError WireArray(0.01, 0.002, 7, 10, copper_props, temperature=20.0, lay_direction=nothing)
    end

    @testset "Basic Functionality" begin
        w = WireArray(0.01, 0.002, 7, 10, copper_props, temperature=20.0, lay_direction=1)
        @test w isa WireArray
        @test w.radius_in ≈ 0.01 atol = TEST_TOL
        @test w.radius_wire ≈ 0.002 atol = TEST_TOL
        @test w.num_wires == 7
        @test w.lay_ratio ≈ 10 atol = TEST_TOL
        @test w.material_props === copper_props
        @test w.temperature ≈ 20.0 atol = TEST_TOL
        @test w.lay_direction == 1
        @test w.cross_section ≈ 7 * π * 0.002^2 atol = TEST_TOL
        # Measurement type
        w2 = WireArray(measurement(0.01, 1e-5), measurement(0.002, 1e-6), 7, 10, copper_props, temperature=measurement(20.0, 0.1), lay_direction=1)
        @test w2 isa WireArray
        @test value(w2.radius_in) ≈ 0.01 atol = TEST_TOL
        @test value(w2.radius_wire) ≈ 0.002 atol = TEST_TOL
        @test value(w2.temperature) ≈ 20.0 atol = TEST_TOL
    end

    @testset "Edge Cases" begin
        # radius_in very close to radius_ext
        w = WireArray(1e-6, 0.002, 7, 10, copper_props, temperature=20.0, lay_direction=1)
        @test w.radius_in ≈ 1e-6 atol = TEST_TOL
        # num_wires = 1 (should set radius_ext = radius_wire)
        w1 = WireArray(0.0, 0.002, 1, 10, copper_props, temperature=20.0, lay_direction=1)
        @test w1.radius_ext ≈ 0.002 atol = TEST_TOL
    end

    @testset "Physical Behavior" begin
        # Resistance should increase with temperature
        w20 = WireArray(0.01, 0.002, 7, 10, copper_props, temperature=20.0, lay_direction=1)
        w80 = WireArray(0.01, 0.002, 7, 10, copper_props, temperature=80.0, lay_direction=1)
        @test w80.resistance > w20.resistance
        # Cross-section increases with wire radius
        w_small = WireArray(0.01, 0.001, 7, 10, copper_props, temperature=20.0, lay_direction=1)
        w_large = WireArray(0.01, 0.003, 7, 10, copper_props, temperature=20.0, lay_direction=1)
        @test w_large.cross_section > w_small.cross_section
    end

    @testset "Type Stability & Promotion" begin
        # All Float64
        w = WireArray(0.01, 0.002, 7, 10.0, copper_props, temperature=20.0, lay_direction=1)
        @test typeof(w.radius_in) == Float64
        # All Measurement
        wM = WireArray(measurement(0.01, 1e-5), measurement(0.002, 1e-6), 7, measurement(10.0, 0.1), copper_props, temperature=measurement(20.0, 0.1), lay_direction=1)
        @test typeof(wM.radius_in) <: Measurement
        # Mixed: radius_in as Measurement
        wMix1 = WireArray(measurement(0.01, 1e-5), 0.002, 7, 10.0, copper_props, temperature=20.0, lay_direction=1)
        @test typeof(wMix1.radius_in) <: Measurement
        # Mixed: lay_ratio as Measurement
        wMix2 = WireArray(0.01, 0.002, 7, measurement(10.0, 0.1), copper_props, temperature=20.0, lay_direction=1)
        @test typeof(wMix2.lay_ratio) <: Measurement
        # material as measurement
        mmat = Material(measurement(1.7241e-8, 1e-10), 1.0, 1.0, 20.0, 0.00393)
        wMix3 = WireArray(0.01, 0.002, 7, 10.0, mmat, temperature=20.0, lay_direction=1)
        @test typeof(wMix3.resistance) <: Measurement
    end

end
