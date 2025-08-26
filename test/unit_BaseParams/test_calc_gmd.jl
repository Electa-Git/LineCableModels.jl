@testitem "BaseParams: calc_gmd unit tests" setup = [defaults, deps_datamodel, defs_materials] begin
    @testset "Basic Functionality" begin
        material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
        wire_array = WireArray(0.01, Diameter(0.002), 7, 10, material_props, temperature=25)
        tubular = Tubular(0.01, 0.02, material_props, temperature=25)
        gmd = calc_gmd(wire_array, tubular)
        @test gmd > 0
        # Symmetry
        gmd2 = calc_gmd(tubular, wire_array)
        @test isapprox(gmd, gmd2, atol=TEST_TOL)
    end

    @testset "Edge Cases" begin
        material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
        # Identical objects (should return outer radius)
        tubular = Tubular(0.01, 0.02, material_props, temperature=25)
        gmd_same = calc_gmd(tubular, tubular)
        @test isapprox(gmd_same, 0.02, atol=TEST_TOL)
        # WireArray with itself
        wire_array = WireArray(0.01, Diameter(0.002), 7, 10, material_props, temperature=25)
        gmd_wa = calc_gmd(wire_array, wire_array)
        @test gmd_wa > 0
    end

    @testset "Numerical Consistency" begin
        material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
        wire_array_f32 = WireArray(Float32(0.01), Diameter(Float32(0.002)), 7, Float32(10), material_props, temperature=25)
        tubular_f32 = Tubular(Float32(0.01), Float32(0.02), material_props, temperature=25)
        gmd_f32 = calc_gmd(wire_array_f32, tubular_f32)
        wire_array_f64 = WireArray(0.01, Diameter(0.002), 7, 10, material_props, temperature=25)
        tubular_f64 = Tubular(0.01, 0.02, material_props, temperature=25)
        gmd_f64 = calc_gmd(wire_array_f64, tubular_f64)
        @test isapprox(gmd_f32, gmd_f64, atol=TEST_TOL)
    end

    @testset "Physical Behavior" begin
        material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
        wa1 = WireArray(0.01, Diameter(0.002), 7, 10, material_props, temperature=25)
        wa2 = WireArray(0.02, Diameter(0.002), 7, 10, material_props, temperature=25)
        tubular = Tubular(0.01, 0.02, material_props, temperature=25)
        gmd1 = calc_gmd(wa1, tubular)
        gmd2 = calc_gmd(wa2, tubular)
        @test gmd2 > gmd1
    end

    @testset "Type Stability & Promotion" begin
        material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
        wa = WireArray(0.01, Diameter(0.002), 7, 10, material_props, temperature=25)
        tub = Tubular(0.01, 0.02, material_props, temperature=25)
        mwa = WireArray(0.01, Diameter(measurement(0.002, 1e-4)), 7, 10, material_props, temperature=25)
        mtub = Tubular(0.01, measurement(0.02, 1e-4), material_props, temperature=25)
        # All Float64
        res1 = calc_gmd(wa, tub)
        @test typeof(res1) == Float64
        # All Measurement
        res2 = calc_gmd(mwa, mtub)
        @test res2 isa Measurement{Float64}
        # Mixed: first argument Measurement
        res3 = calc_gmd(mwa, tub)
        @test res3 isa Measurement{Float64}
        # Mixed: second argument Measurement
        res4 = calc_gmd(wa, mtub)
        @test res4 isa Measurement{Float64}
        mtub_temp = Tubular(0.01, 0.02, material_props, temperature=measurement(25, 1e-4))
        res5 = calc_gmd(wa, mtub_temp)
        @test res5 isa Measurement{Float64}
    end

    @testset "Uncertainty Quantification" begin
        material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
        wa = WireArray(0.01, Diameter(0.002), 7, 10, material_props, temperature=25)
        tub = Tubular(0.01, 0.02, material_props, temperature=25)
        mwa = WireArray(measurement(0.01, 1e-4), Diameter(0.002), 7, 10, material_props, temperature=25)
        mtub = Tubular(0.01, measurement(0.02, 1e-4), material_props, temperature=25)
        gmd = calc_gmd(mwa, mtub)
        @test gmd isa Measurement{Float64}
        @test uncertainty(gmd) > 0
    end
end
