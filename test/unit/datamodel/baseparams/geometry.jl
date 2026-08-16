@testitem "BaseParams / calc_gmd / contracts" tags=[:unit] setup=[
    BaseParamsTestSupport, UseBaseParamsSupport, TestNumerics] begin
    @testset "Basic Functionality" begin
        material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
        circ_strands = CircStrands(
            0.01, 0.001, 7, 10, material_props, temperature = 25)
        tubular = Tubular(0.01, 0.02, material_props, temperature = 25)
        gmd = calc_gmd(circ_strands, tubular)
        @test gmd > 0
        # Symmetry
        gmd2 = calc_gmd(tubular, circ_strands)
        @test isapprox(gmd, gmd2, atol = TestNumerics.absolute_floor(Float64))
    end

    @testset "Edge Cases" begin
        material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
        # Identical objects (should return outer radius)
        tubular = Tubular(0.01, 0.02, material_props, temperature = 25)
        gmd_same = calc_gmd(tubular, tubular)
        @test isapprox(gmd_same, 0.02, atol = TestNumerics.absolute_floor(Float64))
        # CircStrands with itself
        circ_strands = CircStrands(
            0.01, 0.001, 7, 10, material_props, temperature = 25)
        gmd_wa = calc_gmd(circ_strands, circ_strands)
        @test gmd_wa > 0
    end

    @testset "Numerical Consistency" begin
        material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
        circ_strands_f32 = CircStrands(
            Float32(0.01),
            Float32(0.001),
            7,
            Float32(10),
            material_props,
            temperature = 25
        )
        tubular_f32 = Tubular(Float32(0.01), Float32(0.02), material_props, temperature = 25)
        gmd_f32 = calc_gmd(circ_strands_f32, tubular_f32)
        circ_strands_f64 = CircStrands(
            0.01, 0.001, 7, 10, material_props, temperature = 25)
        tubular_f64 = Tubular(0.01, 0.02, material_props, temperature = 25)
        gmd_f64 = calc_gmd(circ_strands_f64, tubular_f64)
        @test isapprox(gmd_f32, gmd_f64, atol = TestNumerics.absolute_floor(Float64))
    end

    @testset "Physical Behavior" begin
        material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
        wa1 = CircStrands(0.01, 0.001, 7, 10, material_props, temperature = 25)
        wa2 = CircStrands(0.02, 0.001, 7, 10, material_props, temperature = 25)
        tubular = Tubular(0.01, 0.02, material_props, temperature = 25)
        gmd1 = calc_gmd(wa1, tubular)
        gmd2 = calc_gmd(wa2, tubular)
        @test gmd2 > gmd1
    end

    @testset "Type Stability & Promotion" begin
        material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
        wa = CircStrands(0.01, 0.001, 7, 10, material_props, temperature = 25)
        tub = Tubular(0.01, 0.02, material_props, temperature = 25)
        mwa = CircStrands(
            0.01,
            measurement(0.001, 5e-5),
            7,
            10,
            material_props,
            temperature = 25
        )
        mtub = Tubular(0.01, measurement(0.02, 1e-4), material_props, temperature = 25)
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
        mtub_temp = Tubular(0.01, 0.02, material_props, temperature = measurement(25, 1e-4))
        res5 = calc_gmd(wa, mtub_temp)
        @test res5 isa Measurement{Float64}
    end

    @testset "Uncertainty Quantification" begin
        material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
        wa = CircStrands(0.01, 0.001, 7, 10, material_props, temperature = 25)
        tub = Tubular(0.01, 0.02, material_props, temperature = 25)
        mwa = CircStrands(
            measurement(0.01, 1e-4),
            0.001,
            7,
            10,
            material_props,
            temperature = 25
        )
        mtub = Tubular(0.01, measurement(0.02, 1e-4), material_props, temperature = 25)
        gmd = calc_gmd(mwa, mtub)
        @test gmd isa Measurement{Float64}
        @test uncertainty(gmd) > 0
    end
end
@testitem "BaseParams / calc_helical_params / contracts" tags=[:unit] setup=[
    BaseParamsTestSupport, UseBaseParamsSupport, TestNumerics] begin
    using Measurements
    # --- Basic Functionality ---
    @testset "Basic Functionality" begin
        r_in = 0.01
        r_ex = 0.015
        lay_ratio = 12.0
        mean_diam, pitch, overlength = calc_helical_params(r_in, r_ex, lay_ratio)
        @test isapprox(mean_diam, 0.025, atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(pitch, 0.3, atol = TestNumerics.absolute_floor(Float64))
        @test overlength > 1.0
    end

    # --- Edge Cases ---
    @testset "Edge Cases" begin
        # Zero lay ratio (pitch_length = 0)
        m, p, o = calc_helical_params(0.01, 0.015, 0.0)
        @test isapprox(m, 0.025, atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(p, 0.0, atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(o, 1.0, atol = TestNumerics.absolute_floor(Float64))

        # Collapsing geometry (r_in == r_ex)
        m2, p2, o2 = calc_helical_params(0.02, 0.02, 10.0)
        @test isapprox(m2, 0.04, atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(p2, 0.4, atol = TestNumerics.absolute_floor(Float64))
        @test o2 > 1.0

        # Very large lay ratio
        m3, p3, o3 = calc_helical_params(0.01, 0.015, 1e6)
        @test isapprox(m3, 0.025, atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(p3, 25000.0, atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(o3, 1.0, atol = TestNumerics.absolute_floor(Float64))
    end

    # --- Numerical Consistency ---
    @testset "Numerical Consistency" begin
        # Float32
        m, p, o = calc_helical_params(Float32(0.01), Float32(0.015), Float32(12.0))
        @test isapprox(m, 0.025, atol = TestNumerics.absolute_floor(Float64))
        @test isapprox(p, 0.3; rtol = sqrt(eps(Float32)))
        @test o > 1.0
    end

    # --- Physical Behavior ---
    @testset "Physical Behavior" begin
        # Increasing lay_ratio increases pitch_length
        _, p1, _ = calc_helical_params(0.01, 0.015, 10.0)
        _, p2, _ = calc_helical_params(0.01, 0.015, 20.0)
        @test p2 > p1
        # Overlength approaches 1 as lay_ratio increases
        _, _, o1 = calc_helical_params(0.01, 0.015, 1e3)
        @test isapprox(o1, 1.0, atol = 1e-5)
    end

    # --- Type Stability & Promotion ---
    @testset "Type Stability & Promotion" begin
        # All Float64
        m, p, o = calc_helical_params(0.01, 0.015, 12.0)
        @test typeof(m) == Float64
        @test typeof(p) == Float64
        @test typeof(o) == Float64
        # All Measurement
        mM, pM, oM = calc_helical_params(
            measurement(0.01, 1e-5),
            measurement(0.015, 1e-5),
            measurement(12.0, 0.1)
        )
        @test mM isa Measurement{Float64}
        @test pM isa Measurement{Float64}
        @test oM isa Measurement{Float64}
        # Mixed: r_in as Measurement
        m1, p1, o1 = calc_helical_params(measurement(0.01, 1e-5), 0.015, 12.0)
        @test m1 isa Measurement{Float64}
        @test p1 isa Measurement{Float64}
        @test o1 isa Measurement{Float64}
        # Mixed: lay_ratio as Measurement
        m2, p2, o2 = calc_helical_params(0.01, 0.015, measurement(12.0, 0.1))
        @test m2 isa Measurement{Float64}
        @test p2 isa Measurement{Float64}
        @test o2 isa Measurement{Float64}
    end

    # --- Uncertainty Quantification ---
    @testset "Uncertainty Quantification" begin
        rin = measurement(0.01, 1e-5)
        rext = measurement(0.015, 1e-5)
        lrat = measurement(12.0, 0.1)
        m, p, o = calc_helical_params(rin, rext, lrat)
        # Check propagated uncertainties are nonzero
        @test uncertainty(m) > 0
        @test uncertainty(p) > 0
        @test uncertainty(o) > 0
    end
end
@testitem "BaseParams / calc_circstrands_coords / contracts" tags=[:unit] setup=[
    BaseParamsTestSupport, UseBaseParamsSupport, TestNumerics] begin
    @testset "Basic Functionality" begin
        @testset "Standard 6-wire array at origin" begin
            let num_wires = 6, radius_wire = 0.001, r_in = 0.01
                lay_radius = r_in + radius_wire # 0.011
                coords = calc_circstrands_coords(num_wires, radius_wire, r_in)

                @test length(coords) == num_wires
                @test coords isa Vector{Tuple{Float64, Float64}}

                # Expected coordinates for a 6-wire array (angle step = π/3)
                expected = [
                    (lay_radius, 0.0), # Angle 0
                    (lay_radius * cos(π / 3), lay_radius * sin(π / 3)),   # Angle π/3
                    (lay_radius * cos(2π / 3), lay_radius * sin(2π / 3)), # Angle 2π/3
                    (-lay_radius, 0.0), # Angle π
                    (lay_radius * cos(4π / 3), lay_radius * sin(4π / 3)), # Angle 4π/3
                    (lay_radius * cos(5π / 3), lay_radius * sin(5π / 3)) # Angle 5π/3
                ]

                @test length(coords) == length(expected)
                for (coord, exp_coord) in zip(coords, expected)
                    @test isapprox(coord[1], exp_coord[1]; atol = TestNumerics.absolute_floor(Float64))
                    @test isapprox(coord[2], exp_coord[2]; atol = TestNumerics.absolute_floor(Float64))
                end
            end
        end

        @testset "4-wire array with non-zero center" begin
            let num_wires = 4, radius_wire = 0.002, r_in = 0.02, C = (0.1, -0.2)
                lay_radius = r_in + radius_wire # 0.022
                coords = calc_circstrands_coords(num_wires, radius_wire, r_in, C)

                @test length(coords) == num_wires

                # Expected coordinates for a 4-wire array (angle step = π/2)
                expected = [
                    (C[1] + lay_radius, C[2]),           # Angle 0
                    (C[1], C[2] + lay_radius),           # Angle π/2
                    (C[1] - lay_radius, C[2]),           # Angle π
                    (C[1], C[2] - lay_radius)           # Angle 3π/2
                ]
                @test length(coords) == length(expected)
                for (coord, exp_coord) in zip(coords, expected)
                    @test isapprox(coord[1], exp_coord[1]; atol = TestNumerics.absolute_floor(Float64))
                    @test isapprox(coord[2], exp_coord[2]; atol = TestNumerics.absolute_floor(Float64))
                end
            end
        end
    end

    @testset "Edge Cases" begin
        @testset "Single wire is always at the center" begin
            # A single wire's lay radius is defined as 0.
            coords = calc_circstrands_coords(1, 0.001, 0.01)
            @test coords == [(0.0, 0.0)]

            C = (10.0, -20.0)
            coords_C = calc_circstrands_coords(1, 0.001, 0.01, C)
            @test coords_C == [C]
        end

        @testset "Zero wires returns an empty vector" begin
            coords = calc_circstrands_coords(0, 0.001, 0.01)
            @test isempty(coords)
            @test coords isa Vector
        end

        @testset "Zero radii places all wires at the center" begin
            # If lay radius is zero, all wires should be at the center C.
            num_wires = 7
            coords = calc_circstrands_coords(num_wires, 0.0, 0.0)
            @test length(coords) == num_wires
            @test all(c -> c == (0.0, 0.0), coords)

            C = (1.0, 1.0)
            coords_C = calc_circstrands_coords(num_wires, 0.0, 0.0, C)
            @test length(coords_C) == num_wires
            @test all(c -> c == C, coords_C)
        end
    end

    @testset "Type Stability and Promotion" begin
        @testset "Base case: Float64 inputs" begin
            coords = calc_circstrands_coords(6, 0.001, 0.01)
            @test coords isa Vector{Tuple{Float64, Float64}}
            @test eltype(first(coords)) == Float64
        end

        @testset "Fully promoted: All inputs are Measurement" begin
            num_wires = 3
            rw = 0.001 ± 0.0001
            ri = 0.01 ± 0.0002
            C = (0.1 ± 0.01, -0.2 ± 0.02)
            coords = calc_circstrands_coords(num_wires, rw, ri, C)

            @test coords isa Vector{Tuple{Measurement{Float64}, Measurement{Float64}}}
            @test eltype(first(coords)) == Measurement{Float64}

            # Check value and uncertainty propagation for the first wire (angle=0)
            lay_radius = rw + ri
            expected_x = C[1] + lay_radius
            expected_y = C[2] # sin(0) is 0, so lay_radius term is zero

            @test coords[1][1] ≈ expected_x
            @test coords[1][2] ≈ expected_y
        end

        @testset "Mixed types: radius_wire is Measurement" begin
            num_wires = 4
            rw = 0.001 ± 0.0001
            ri = 0.01 # Float64
            C = (0.1, -0.2) # Tuple{Float64, Float64}
            coords = calc_circstrands_coords(num_wires, rw, ri, C = C)

            @test coords isa Vector{Tuple{Measurement{Float64}, Measurement{Float64}}}
            lay_radius_val = Measurements.value(rw) + ri

            # Wire 1 (angle 0)
            @test Measurements.value(coords[1][1]) ≈ C[1] + lay_radius_val atol = TestNumerics.absolute_floor(Float64)
            @test Measurements.value(coords[1][2]) ≈ C[2] atol = TestNumerics.absolute_floor(Float64)
            @test Measurements.uncertainty(coords[1][1]) > 0
            @test Measurements.uncertainty(coords[1][2]) == 0 # sin(0) = 0, no uncertainty propagation

            # Wire 2 (angle π/2)
            @test Measurements.value(coords[2][1]) ≈ C[1] atol = TestNumerics.absolute_floor(Float64)
            @test Measurements.value(coords[2][2]) ≈ C[2] + lay_radius_val atol = TestNumerics.absolute_floor(Float64)
            @test isapprox(Measurements.uncertainty(coords[2][1]), 0,
                atol = TestNumerics.absolute_floor(Float64)) # cos(π/2) = 0, no uncertainty propagation
            @test Measurements.uncertainty(coords[2][2]) > 0
        end

        @testset "Mixed types: r_in is Measurement" begin
            num_wires = 4
            rw = 0.001 # Float64
            ri = 0.01 ± 0.0002
            coords = calc_circstrands_coords(num_wires, rw, ri)

            @test coords isa Vector{Tuple{Measurement{Float64}, Measurement{Float64}}}
            lay_radius_uncert = Measurements.uncertainty(ri)
            @test Measurements.uncertainty(coords[1][1]) ≈ lay_radius_uncert atol = TestNumerics.absolute_floor(Float64)
        end

        @testset "Mixed types: Center C is Measurement" begin
            num_wires = 4
            rw = 0.001 # Float64
            ri = 0.01 # Float64
            C = (0.1 ± 0.01, -0.2 ± 0.02)
            # Use keyword argument version to test the helper method
            coords = calc_circstrands_coords(num_wires, rw, ri; C = C)

            @test coords isa Vector{Tuple{Measurement{Float64}, Measurement{Float64}}}
            lay_radius = rw + ri

            # Check uncertainty propagation from center C
            @test Measurements.value(coords[1][1]) ≈ Measurements.value(C[1]) + lay_radius atol = TestNumerics.absolute_floor(Float64)
            @test Measurements.value(coords[1][2]) ≈ Measurements.value(C[2]) atol = TestNumerics.absolute_floor(Float64)
            @test Measurements.uncertainty(coords[1][1]) ≈ Measurements.uncertainty(C[1]) atol = TestNumerics.absolute_floor(Float64)
            @test Measurements.uncertainty(coords[1][2]) ≈ Measurements.uncertainty(C[2]) atol = TestNumerics.absolute_floor(Float64)
        end
    end
end
@testitem "BaseParams / calc_circstrands_gmr / contracts" tags=[:unit] setup=[
    BaseParamsTestSupport, UseBaseParamsSupport, TestNumerics] begin
    using Measurements: measurement, value, uncertainty

    @testset "Basic Functionality" begin
        # Example from docstring
        lay_rad = 0.05
        N = 7
        rad_wire = 0.002
        mu_r = 1.0
        gmr = calc_circstrands_gmr(lay_rad, N, rad_wire, mu_r)
        expected = exp((log(rad_wire * exp(-mu_r / 4) * N * lay_rad^(N - 1)) / N))
        @test isapprox(gmr, expected; atol = TestNumerics.absolute_floor(Float64))
        @test gmr > 0
    end

    @testset "Edge Cases" begin
        # N = 1 (single wire)
        lay_rad = 0.05
        N = 1
        rad_wire = 0.002
        mu_r = 1.0
        gmr = calc_circstrands_gmr(lay_rad, N, rad_wire, mu_r)
        expected = rad_wire * exp(-mu_r / 4)
        @test isapprox(gmr, expected; atol = TestNumerics.absolute_floor(Float64))

        # mu_r = 0 (non-magnetic)
        lay_rad = 0.05
        N = 7
        rad_wire = 0.002
        mu_r = 0.0
        gmr = calc_circstrands_gmr(lay_rad, N, rad_wire, mu_r)
        expected = exp((log(rad_wire * N * lay_rad^(N - 1)) / N))
        @test isapprox(gmr, expected; atol = TestNumerics.absolute_floor(Float64))

        # rad_wire = 0 (degenerate wire)
        lay_rad = 0.05
        N = 7
        rad_wire = 0.0
        mu_r = 1.0
        gmr = calc_circstrands_gmr(lay_rad, N, rad_wire, mu_r)
        @test gmr == 0.0

        # lay_rad = 0 (all wires at center)
        lay_rad = 0.0
        N = 7
        rad_wire = 0.002
        mu_r = 1.0
        gmr = calc_circstrands_gmr(lay_rad, N, rad_wire, mu_r)
        expected = exp((log(rad_wire * exp(-mu_r / 4) * N * 0.0^(N - 1)) / N))
        @test gmr == 0.0
    end

    @testset "Numerical Consistency" begin
        # Float64
        gmr1 = calc_circstrands_gmr(0.05, 7, 0.002, 1.0)
        # Measurement{Float64}
        gmr2 = calc_circstrands_gmr(measurement(0.05, 1e-4), 7, 0.002, 1.0)
        @test isapprox(value(gmr2), gmr1; atol = TestNumerics.absolute_floor(Float64))
        @test uncertainty(gmr2) > 0
    end

    @testset "Physical Behavior" begin
        # GMR increases with lay_rad
        gmr1 = calc_circstrands_gmr(0.01, 7, 0.002, 1.0)
        gmr2 = calc_circstrands_gmr(0.05, 7, 0.002, 1.0)
        @test gmr2 > gmr1
        # GMR decreases with mu_r
        gmr1 = calc_circstrands_gmr(0.05, 7, 0.002, 0.5)
        gmr2 = calc_circstrands_gmr(0.05, 7, 0.002, 2.0)
        @test gmr2 < gmr1
    end

    @testset "Type Stability & Promotion" begin
        # All Float64
        gmr = calc_circstrands_gmr(0.05, 7, 0.002, 1.0)
        @test typeof(gmr) == Float64
        # All Measurement
        gmr = calc_circstrands_gmr(measurement(0.05, 1e-4), 7, measurement(0.002, 1e-5), measurement(1.0, 0.01))
        @test gmr isa Measurement{Float64}
        # Mixed: lay_rad as Measurement
        gmr = calc_circstrands_gmr(measurement(0.05, 1e-4), 7, 0.002, 1.0)
        @test gmr isa Measurement{Float64}
        # Mixed: rad_wire as Measurement
        gmr = calc_circstrands_gmr(0.05, 7, measurement(0.002, 1e-5), 1.0)
        @test gmr isa Measurement{Float64}
        # Mixed: mu_r as Measurement
        gmr = calc_circstrands_gmr(0.05, 7, 0.002, measurement(1.0, 0.01))
        @test gmr isa Measurement{Float64}
    end

    @testset "Uncertainty Quantification" begin
        lay_rad = measurement(0.05, 1e-4)
        rad_wire = measurement(0.002, 1e-5)
        mu_r = measurement(1.0, 0.01)
        gmr = calc_circstrands_gmr(lay_rad, 7, rad_wire, mu_r)
        @test gmr isa Measurement{Float64}
        @test uncertainty(gmr) > 0
    end
end
@testitem "BaseParams / calc_equivalent_gmr / contracts" tags=[:unit] setup=[
    BaseParamsTestSupport, UseBaseParamsSupport, TestNumerics] begin
    @testset "Basic Functionality" begin
        # Example from docstring
        material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
        strip = Strip(0.01, 0.05, 10, material_props; thickness = 0.002)
        circstrands = CircStrands(0.02, 0.002, 7, 15, material_props)
        gmr_eq = calc_equivalent_gmr(strip, circstrands)
        @test gmr_eq > 0
    end

    @testset "Edge Cases" begin
        # Identical layers (should reduce to geometric mean)
        material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
        part1 = CircStrands(0.01, 0.002, 7, 10, material_props)
        part2 = CircStrands(0.01, 0.002, 7, 10, material_props)
        gmr_eq = calc_equivalent_gmr(part1, part2)
        @test gmr_eq > 0
        # Very large cross-section for new_layer
        big_layer = CircStrands(0.02, 0.002, 7, 1e6, material_props)
        gmr_eq2 = calc_equivalent_gmr(part1, big_layer)
        @test gmr_eq2 > 0
    end

    @testset "Numerical Consistency" begin
        # Float32 vs Float64
        material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
        part1f32 = CircStrands(
            Float32(0.01), Float32(0.002), 7, Float32(10), material_props)
        part2f32 = CircStrands(
            Float32(0.02), Float32(0.002), 7, Float32(15), material_props)
        gmr_eq_f32 = calc_equivalent_gmr(part1f32, part2f32)
        part1f64 = CircStrands(0.01, 0.002, 7, 10, material_props)
        part2f64 = CircStrands(0.02, 0.002, 7, 15, material_props)
        gmr_eq_f64 = calc_equivalent_gmr(part1f64, part2f64)
        @test isapprox(gmr_eq_f32, gmr_eq_f64, atol = TestNumerics.absolute_floor(Float64))
    end

    @testset "Physical Behavior" begin
        # Equivalent GMR increases as GMD increases
        material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
        part1 = CircStrands(0.01, 0.002, 7, 10, material_props)
        part2 = CircStrands(0.02, 0.002, 7, 15, material_props)
        part3 = CircStrands(0.03, 0.002, 7, 15, material_props)
        gmr_eq1 = calc_equivalent_gmr(part1, part2)
        gmr_eq2 = calc_equivalent_gmr(part1, part3)
        @test gmr_eq2 > gmr_eq1
    end

    @testset "Type Stability & Promotion" begin
        using Measurements
        material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
        part1 = CircStrands(0.01, 0.002, 7, 10, material_props)
        part2 = CircStrands(0.02, 0.002, 7, 15, material_props)
        mpart1 = CircStrands(measurement(0.01, 1e-4), 0.002, 7, 10, material_props)
        mpart2 = CircStrands(0.02, measurement(0.002, 1e-4), 7, 15, material_props)
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
        part1 = CircStrands(0.01, 0.002, 7, 10, material_props)
        part2 = CircStrands(0.02, 0.002, 7, 15, material_props)
        mpart1 = CircStrands(measurement(0.01, 1e-4), 0.002, 7, 10, material_props)
        mpart2 = CircStrands(0.02, measurement(0.002, 1e-4), 7, 15, material_props)
        gmr_eq = calc_equivalent_gmr(mpart1, mpart2)
        @test gmr_eq isa Measurement{Float64}
        @test uncertainty(gmr_eq) > 0
    end
end
@testitem "BaseParams / calc_tubular_gmr / contracts" tags=[:unit] setup=[
    BaseParamsTestSupport, UseBaseParamsSupport, TestNumerics] begin
    using Measurements: measurement, value, uncertainty

    @testset "Basic Functionality" begin
        # Example from docstring
        r_ex = 0.02
        r_in = 0.01
        mu_r = 1.0
        gmr = calc_tubular_gmr(r_ex, r_in, mu_r)
        # Manual calculation for expected value
        term1 = (r_in^4 / (r_ex^2 - r_in^2)^2) * log(r_ex / r_in)
        term2 = (3 * r_in^2 - r_ex^2) / (4 * (r_ex^2 - r_in^2))
        Lin = (μ₀ * mu_r / (2 * π)) * (term1 - term2)
        expected = exp(log(r_ex) - (2 * π / μ₀) * Lin)
        @test isapprox(gmr, expected; atol = TestNumerics.absolute_floor(Float64))
        @test gmr > 0
        @test_throws ArgumentError calc_tubular_gmr(r_in, r_ex, mu_r)
        @test_throws ArgumentError calc_tubular_gmr(0.0, r_in, mu_r)
    end

    @testset "Edge Cases" begin
        # Thin shell: r_ex ≈ r_in
        r_ex = 0.01
        r_in = 0.01
        mu_r = 1.0
        gmr = calc_tubular_gmr(r_ex, r_in, mu_r)
        @test isapprox(gmr, r_ex; atol = TestNumerics.absolute_floor(Float64))

        # Infinitely thick tube: r_in ≫ 0, r_in / r_ex ≈ 0
        r_ex = 1.0
        r_in = 1e-12
        mu_r = 1.0
        gmr = calc_tubular_gmr(r_ex, r_in, mu_r)
        @test isapprox(gmr, 0.7788; atol = 1e-4)

        # r_in = 0 (solid cylinder)
        r_ex = 0.02
        r_in = 0.0
        mu_r = 1.0
        gmr = calc_tubular_gmr(r_ex, r_in, mu_r)
        @test isapprox(gmr, 0.7788 * r_ex; atol = 1e-4)

        # r_ex < r_in (should throw)
        r_ex = 0.01
        r_in = 0.02
        mu_r = 1.0
        @test_throws ArgumentError calc_tubular_gmr(r_ex, r_in, mu_r)
    end

    @testset "Numerical Consistency" begin
        # Float64
        gmr1 = calc_tubular_gmr(0.02, 0.01, 1.0)
        # Measurement{Float64}
        gmr2 = calc_tubular_gmr(
            measurement(0.02, 1e-4),
            measurement(0.01, 1e-4),
            measurement(1.0, 0.01)
        )
        @test isapprox(value(gmr2), gmr1; atol = TestNumerics.absolute_floor(Float64))
        @test uncertainty(gmr2) > 0
    end

    @testset "Physical Behavior" begin
        # GMR increases with r_ex
        gmr1 = calc_tubular_gmr(0.01, 0.005, 1.0)
        gmr2 = calc_tubular_gmr(0.02, 0.005, 1.0)
        @test gmr2 > gmr1
        # GMR decreases with mu_r
        gmr1 = calc_tubular_gmr(0.02, 0.01, 0.5)
        gmr2 = calc_tubular_gmr(0.02, 0.01, 2.0)
        @test gmr2 < gmr1
    end

    @testset "Type Stability & Promotion" begin
        # All Float64
        gmr = calc_tubular_gmr(0.02, 0.01, 1.0)
        @test typeof(gmr) == Float64
        # All Measurement
        gmr = calc_tubular_gmr(
            measurement(0.02, 1e-4),
            measurement(0.01, 1e-4),
            measurement(1.0, 0.01)
        )
        @test gmr isa Measurement{Float64}
        # Mixed: r_ex as Measurement
        gmr = calc_tubular_gmr(measurement(0.02, 1e-4), 0.01, 1.0)
        @test gmr isa Measurement{Float64}
        # Mixed: r_in as Measurement
        gmr = calc_tubular_gmr(0.02, measurement(0.01, 1e-4), 1.0)
        @test gmr isa Measurement{Float64}
        # Mixed: mu_r as Measurement
        gmr = calc_tubular_gmr(0.02, 0.01, measurement(1.0, 0.01))
        @test gmr isa Measurement{Float64}
    end

    @testset "Uncertainty Quantification" begin
        r_ex = measurement(0.02, 1e-4)
        r_in = measurement(0.01, 1e-4)
        mu_r = measurement(1.0, 0.01)
        gmr = calc_tubular_gmr(r_ex, r_in, mu_r)
        @test gmr isa Measurement{Float64}
        @test uncertainty(gmr) > 0
    end
end
@testitem "BaseParams / calc_solenoid_correction / contracts" tags=[:unit] setup=[
    BaseParamsTestSupport, UseBaseParamsSupport, TestNumerics] begin
    @testset "Basic Functionality" begin
        # Example from docstring: 10 turns/m, conductor radius 5 mm, insulator radius 10 mm
        result = calc_solenoid_correction(10.0, 0.005, 0.01)
        expected = 1.0 + 2 * 10.0^2 * pi^2 * (0.01^2 - 0.005^2) / log(0.01 / 0.005)
        @test isapprox(result, expected; atol = TestNumerics.absolute_floor(Float64))
        @test result > 1.0

        # Non-helical cable (NaN turns)
        result = calc_solenoid_correction(NaN, 0.005, 0.01)
        @test result == 1.0
    end

    @testset "Edge Cases" begin
        # Zero turns (should be 1.0)
        result = calc_solenoid_correction(0.0, 0.005, 0.01)
        @test isapprox(result, 1.0; atol = TestNumerics.absolute_floor(Float64))

        # Collapsing geometry: radii nearly equal
        result = calc_solenoid_correction(10.0, 0.01, 0.010001)
        @test result > 1.0

        # Very large number of turns
        result = calc_solenoid_correction(1e6, 0.005, 0.01)
        @test result > 1e9

        # Inf/NaN radii
        @test isnan(calc_solenoid_correction(10.0, NaN, 0.01))
        @test isnan(calc_solenoid_correction(10.0, 0.005, NaN))
        @test isinf(calc_solenoid_correction(10.0, 0.0, 0.01)) == false  # log(0.01/0) = Inf, but numerator is finite
    end

    @testset "Numerical Consistency" begin
        # Float32 vs Float64
        r = calc_solenoid_correction(Float32(10.0), Float32(0.005), Float32(0.01))
        d = calc_solenoid_correction(10.0, 0.005, 0.01)
        @test isapprox(r, d; rtol = sqrt(eps(Float32)))
    end

    @testset "Physical Behavior" begin
        # Correction increases with more turns
        c1 = calc_solenoid_correction(5.0, 0.005, 0.01)
        c2 = calc_solenoid_correction(10.0, 0.005, 0.01)
        @test c2 > c1
        # Correction increases with larger insulator radius
        c3 = calc_solenoid_correction(10.0, 0.005, 0.02)
        @test c3 > c2
    end

    @testset "Type Stability & Promotion" begin
        using Measurements
        # All Float64
        r1 = calc_solenoid_correction(10.0, 0.005, 0.01)
        @test typeof(r1) == Float64
        # All Measurement
        r2 = calc_solenoid_correction(measurement(10.0, 0.1), measurement(0.005, 1e-5), measurement(0.01, 1e-5))
        @test r2 isa Measurement{Float64}
        # Mixed: num_turns as Measurement
        r3 = calc_solenoid_correction(measurement(10.0, 0.1), 0.005, 0.01)
        @test r3 isa Measurement{Float64}
        # Mixed: radius_ext_con as Measurement
        r4 = calc_solenoid_correction(10.0, measurement(0.005, 1e-5), 0.01)
        @test r4 isa Measurement{Float64}
        # Mixed: radius_ext_ins as Measurement
        r5 = calc_solenoid_correction(10.0, 0.005, measurement(0.01, 1e-5))
        @test r5 isa Measurement{Float64}
    end

    @testset "Uncertainty Quantification" begin
        using Measurements
        num_turns = measurement(10.0, 0.1)
        radius_ext_con = measurement(0.005, 1e-5)
        radius_ext_ins = measurement(0.01, 1e-5)
        result = calc_solenoid_correction(num_turns, radius_ext_con, radius_ext_ins)
        # Should propagate uncertainty
        @test result isa Measurement{Float64}
        # Uncertainty should be nonzero
        @test uncertainty(result) > 0
    end
end
