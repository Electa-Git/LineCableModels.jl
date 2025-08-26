@testitem "BaseParams: calc_equivalent_gmr unit tests" setup = [defaults, deps_datamodel, defs_materials] begin
    @testset "Basic Functionality" begin
        # Example from docstring
        material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
        strip = Strip(0.01, Thickness(0.002), 0.05, 10, material_props)
        wirearray = WireArray(0.02, 0.002, 7, 15, material_props)
        gmr_eq = calc_equivalent_gmr(strip, wirearray)
        @test gmr_eq > 0
    end

    @testset "Edge Cases" begin
        # Identical layers (should reduce to geometric mean)
        material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
        part1 = WireArray(0.01, 0.002, 7, 10, material_props)
        part2 = WireArray(0.01, 0.002, 7, 10, material_props)
        gmr_eq = calc_equivalent_gmr(part1, part2)
        @test gmr_eq > 0
        # Very large cross-section for new_layer
        big_layer = WireArray(0.02, 0.002, 7, 1e6, material_props)
        gmr_eq2 = calc_equivalent_gmr(part1, big_layer)
        @test gmr_eq2 > 0
    end

    @testset "Numerical Consistency" begin
        # Float32 vs Float64
        material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
        part1f32 = WireArray(Float32(0.01), Float32(0.002), 7, Float32(10), material_props)
        part2f32 = WireArray(Float32(0.02), Float32(0.002), 7, Float32(15), material_props)
        gmr_eq_f32 = calc_equivalent_gmr(part1f32, part2f32)
        part1f64 = WireArray(0.01, 0.002, 7, 10, material_props)
        part2f64 = WireArray(0.02, 0.002, 7, 15, material_props)
        gmr_eq_f64 = calc_equivalent_gmr(part1f64, part2f64)
        @test isapprox(gmr_eq_f32, gmr_eq_f64, atol=TEST_TOL)
    end

    @testset "Physical Behavior" begin
        # Equivalent GMR increases as GMD increases
        material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
        part1 = WireArray(0.01, 0.002, 7, 10, material_props)
        part2 = WireArray(0.02, 0.002, 7, 15, material_props)
        part3 = WireArray(0.03, 0.002, 7, 15, material_props)
        gmr_eq1 = calc_equivalent_gmr(part1, part2)
        gmr_eq2 = calc_equivalent_gmr(part1, part3)
        @test gmr_eq2 > gmr_eq1
    end

    @testset "Type Stability & Promotion" begin
        using Measurements
        material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
        part1 = WireArray(0.01, 0.002, 7, 10, material_props)
        part2 = WireArray(0.02, 0.002, 7, 15, material_props)
        mpart1 = WireArray(measurement(0.01, 1e-4), 0.002, 7, 10, material_props)
        mpart2 = WireArray(0.02, measurement(0.002, 1e-4), 7, 15, material_props)
        # All Float64
        res1 = calc_equivalent_gmr(part1, part2)
        @test typeof(res1) == Float64
        # All Measurement
        res2 = calc_equivalent_gmr(mpart1, mpart2)
        @test res2 isa Measurement{Float64}
        # Mixed: first argument Measurement
        res3 = calc_equivalent_gmr(mpart1, part2)
        @test res3 isa Measurement{Float64}
        # Mixed: second argument Measurement
        res4 = calc_equivalent_gmr(part1, mpart2)
        @test res4 isa Measurement{Float64}
    end

    @testset "Uncertainty Quantification" begin
        using Measurements
        material_props = Material(1.7241e-8, 1.0, 0.999994, measurement(20, 10), 0.00393)
        part1 = WireArray(0.01, 0.002, 7, 10, material_props)
        part2 = WireArray(0.02, 0.002, 7, 15, material_props)
        mpart1 = WireArray(measurement(0.01, 1e-4), 0.002, 7, 10, material_props)
        mpart2 = WireArray(0.02, measurement(0.002, 1e-4), 7, 15, material_props)
        gmr_eq = calc_equivalent_gmr(mpart1, mpart2)
        @test gmr_eq isa Measurement{Float64}
        @test uncertainty(gmr_eq) > 0
    end
end
