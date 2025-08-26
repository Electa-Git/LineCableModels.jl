@testitem "BaseParams module" setup = [defaults, defs_materials] begin
    @testset "Temperature correction" begin
        alpha = 0.004
        T0 = 20.0
        # Correction factor should be 1 at reference temperature
        @test calc_temperature_correction(alpha, T0, T0) ≈ 1.0 atol = TEST_TOL
        # Test T > T0 
        @test calc_temperature_correction(alpha, 30.0, T0) ≈ (1 + alpha * (30.0 - T0)) atol = TEST_TOL
        # Test T < T0
        @test calc_temperature_correction(alpha, 10.0, T0) ≈ (1 + alpha * (10.0 - T0)) atol = TEST_TOL
        # No correction if alpha is zero
        @test calc_temperature_correction(0.0, 50.0, T0) ≈ 1.0 atol = TEST_TOL
    end

    @testset "Parallel impedance calculations" begin
        # Parallel equivalent of two equal resistors
        @test calc_parallel_equivalent(10.0, 10.0) ≈ 5.0 atol = TEST_TOL
        # Parallel equivalent of two different resistors
        @test calc_parallel_equivalent(10.0, 5.0) ≈ (10.0 * 5.0) / (10.0 + 5.0) atol =
            TEST_TOL
        # Adding infinite resistance changes nothing
        @test calc_parallel_equivalent(10.0, Inf) ≈ 10.0 atol = TEST_TOL
        # Adding zero resistance results in zero (short circuit)
        @test calc_parallel_equivalent(10.0, 0.0) ≈ 0.0 atol = TEST_TOL

        # Complex numbers (impedances)
        Z1 = 3.0 + 4.0im
        Z2 = 8.0 - 6.0im
        Zeq_expected = (Z1 * Z2) / (Z1 + Z2)
        # Parallel equivalent of complex impedances
        @test calc_parallel_equivalent(Z1, Z2) ≈ Zeq_expected atol = TEST_TOL
        # Parallel equivalent of two equal complex impedances
        @test calc_parallel_equivalent(Z1, Z1) ≈ Z1 / 2 atol = TEST_TOL
    end

    @testset "Equivalent temperature coefficient" begin
        alpha1, R1 = 0.004, 10.0
        alpha2, R2 = 0.003, 5.0
        expected_alpha = (alpha1 * R2 + alpha2 * R1) / (R1 + R2)
        @test calc_equivalent_alpha(alpha1, R1, alpha2, R2) ≈ expected_alpha atol = TEST_TOL
        # Equivalent alpha of identical conductors
        @test calc_equivalent_alpha(alpha1, R1, alpha1, R1) ≈ alpha1 atol = TEST_TOL
        # Check symmetry
        @test calc_equivalent_alpha(alpha1, R1, alpha2, R2) ≈
              calc_equivalent_alpha(alpha2, R2, alpha1, R1) atol = TEST_TOL
    end

    @testset "Resistance calculations" begin
        # Using Copper properties
        rho = copper_props.rho
        alpha = copper_props.alpha
        T0 = copper_props.T0
        T = T0 # Test at reference temperature first

        # calc_tubular_resistance
        r_in, r_ext = 0.01, 0.02
        area_tube = π * (r_ext^2 - r_in^2)
        #  Tubular resistance at T0
        @test calc_tubular_resistance(r_in, r_ext, rho, alpha, T0, T) ≈ rho / area_tube atol =
            TEST_TOL
        # Solid conductor
        area_solid = π * r_ext^2
        # Solid conductor resistance (r_in=0)
        @test calc_tubular_resistance(0.0, r_ext, rho, alpha, T0, T) ≈ rho / area_solid atol =
            TEST_TOL
        # Temperature dependence
        T_hot = 70.0
        k = calc_temperature_correction(alpha, T_hot, T0)
        # Tubular resistance temperature dependence
        @test calc_tubular_resistance(r_in, r_ext, rho, alpha, T0, T_hot) ≈
              (rho / area_tube) * k atol = TEST_TOL
        # Thin tube limit (resistance should increase) - check relative magnitude
        r_in_thin = r_ext * 0.999
        # Thin tube has higher resistance
        @test calc_tubular_resistance(r_in_thin, r_ext, rho, alpha, T0, T) >
              calc_tubular_resistance(r_in, r_ext, rho, alpha, T0, T)

        # calc_strip_resistance
        thickness, width = 0.002, 0.05
        area_strip = thickness * width
        # Strip resistance at T0
        @test calc_strip_resistance(thickness, width, rho, alpha, T0, T) ≈ rho / area_strip atol =
            TEST_TOL
        # Strip resistance temperature dependence
        @test calc_strip_resistance(thickness, width, rho, alpha, T0, T_hot) ≈
              (rho / area_strip) * k atol = TEST_TOL
    end

    @testset "Helical parameters correction" begin
        r_in, r_ext = 0.01, 0.015
        mean_diam_expected = r_in + r_ext # 0.025
        lay_ratio = 12.0
        pitch_expected = lay_ratio * mean_diam_expected # 12.0 * 0.025 = 0.3

        mean_diam, pitch, overlength = calc_helical_params(r_in, r_ext, lay_ratio)
        @test mean_diam ≈ mean_diam_expected atol = TEST_TOL
        @test pitch ≈ pitch_expected atol = TEST_TOL
        @test overlength ≈ sqrt(1 + (π * mean_diam_expected / pitch_expected)^2) atol =
            TEST_TOL
        # Overlength factor must be > 1 for finite lay ratio
        @test overlength > 1.0

        # Edge case: No twist (infinite pitch length)
        mean_diam_no, pitch_no, overlength_no = calc_helical_params(r_in, r_ext, 0.0)
        # Note: lay_ratio=0 implies pitch=0 in the code, which makes overlength=1
        @test mean_diam_no ≈ mean_diam_expected atol = TEST_TOL
        @test pitch_no == 0.0
        # Overlength factor must be 1 for zero lay ratio (infinite pitch)
        @test overlength_no ≈ 1.0 atol = TEST_TOL
    end

    @testset "GMR calculations & consistency" begin
        r_in, r_ext = 0.01, 0.02
        mu_r = 1.0 # Non-magnetic

        # calc_tubular_gmr
        gmr_tube = calc_tubular_gmr(r_ext, r_in, mu_r)
        # GMR of tube should be less than outer radius
        @test gmr_tube < r_ext
        # GMR must be positive
        @test gmr_tube > 0

        # Solid conductor GMR
        gmr_solid_expected = r_ext * exp(-mu_r / 4.0)
        gmr_solid_calc = calc_tubular_gmr(r_ext, 0.0, mu_r)
        # Solid conductor GMR (analytical)
        @test gmr_solid_calc ≈ gmr_solid_expected atol = TEST_TOL

        # Thin shell GMR
        gmr_shell_calc = calc_tubular_gmr(r_ext, r_ext * (1 - 1e-12), mu_r) # Approx thin shell
        # Thin shell GMR approaches outer radius # Relax tolerance slightly
        @test gmr_shell_calc ≈ r_ext atol = 1e-5

        # Magnetic material
        mu_r_mag = 100.0
        gmr_solid_mag = calc_tubular_gmr(r_ext, 0.0, mu_r_mag)
        gmr_solid_mag_expected = r_ext * exp(-mu_r_mag / 4.0)
        # Solid conductor GMR with mu_r > 1
        @test gmr_solid_mag ≈ gmr_solid_mag_expected atol = TEST_TOL
        # Higher mu_r should decrease GMR
        @test gmr_solid_mag < gmr_solid_calc

        # Error handling
        @test_throws ArgumentError calc_tubular_gmr(r_in, r_ext, mu_r) # Should throw error if r_ext < r_in

        # calc_equivalent_mu (inverse consistency)
        @test calc_equivalent_mu(gmr_tube, r_ext, r_in) ≈ mu_r atol = TEST_TOL #  Inverse check: mu_r from tubular GMR
        @test calc_equivalent_mu(gmr_solid_calc, r_ext, 0.0) ≈ mu_r atol = TEST_TOL #  Inverse check: mu_r from solid GMR
        @test calc_equivalent_mu(gmr_solid_mag, r_ext, 0.0) ≈ mu_r_mag atol = TEST_TOL #  Inverse check: magnetic mu_r from solid GMR
        @test_throws ArgumentError calc_equivalent_mu(gmr_tube, r_in, r_ext) # Should throw error if r_ext < r_in

        # calc_wirearray_gmr
        wire_rad = 0.001
        num_wires = 7
        layout_rad = 0.005 # Center-to-center radius
        gmr_array = calc_wirearray_gmr(layout_rad, num_wires, wire_rad, mu_r)
        # Single wire case should match solid wire GMR
        @test gmr_array > 0

        gmr_single_wire_array = calc_wirearray_gmr(0.0, 1, wire_rad, mu_r) # Layout radius irrelevant for N=1
        gmr_single_wire_solid = calc_tubular_gmr(wire_rad, 0.0, mu_r)
        # GMR of 1-wire array matches solid wire GMR
        @test gmr_single_wire_array ≈ gmr_single_wire_solid atol = TEST_TOL

    end

    @testset "GMD and equivalent GMR" begin
        # Need some cable parts
        part1_solid = Tubular(0.0, 0.01, copper_props) # Solid conductor r=1cm
        part2_tubular = Tubular(0.015, 0.02, copper_props) # Tubular conductor, separate
        part3_wirearray = WireArray(0.03, 0.002, 7, 10.0, aluminum_props) # Wire array, separate

        # calc_gmd
        # Case 1: Two separate solid/tubular conductors (distance between centers)
        # Place part2 at (d, 0) relative to part1 at (0,0)
        d = 0.1 # 10 cm separation
        # GMD calculation for simple geometries relies on center-to-center distance if not wire arrays
        # This test might be trivial for Tubular/Tubular if code assumes center-to-center
        # Let's test Tubular vs WireArray where sub-elements exist
        gmd_1_3 = calc_gmd(part1_solid, part3_wirearray) # Should be approx layout_radius of part3 (0.03 + 0.002) if part1 is at center
        # GMD between central solid and wire array approx layout radius
        @test gmd_1_3 ≈ (0.03 + 0.002) atol = 1e-4

        # Case 2: Concentric Tubular Conductors (test based on comment in code)
        part_inner = Tubular(0.01, 0.02, copper_props)
        part_outer = Tubular(0.02, 0.03, copper_props) # Directly outside inner part
        # If truly concentric, d_ij = 0 for internal logic, should return max(r_ext1, r_ext2)
        gmd_concentric = calc_gmd(part_inner, part_outer)
        # GMD of concentric tubular conductors
        @test gmd_concentric ≈ part_outer.radius_ext atol = TEST_TOL

        # calc_equivalent_gmr
        # Create a conductor group to test adding layers
        core = ConductorGroup(part1_solid)
        layer2 = Tubular(core.radius_ext, 0.015, aluminum_props) # Add tubular layer outside
        beta = core.cross_section / (core.cross_section + layer2.cross_section)
        gmd_core_layer2 = calc_gmd(core.layers[end], layer2) # GMD between solid core and new layer

        gmr_eq_expected =
            (core.gmr^(beta^2)) * (layer2.gmr^((1 - beta)^2)) *
            (gmd_core_layer2^(2 * beta * (1 - beta)))
        gmr_eq_calc = calc_equivalent_gmr(core, layer2) # Test the function directly
        @test gmr_eq_calc ≈ gmr_eq_expected atol = TEST_TOL

        # Test adding a WireArray layer
        layer3_wa = WireArray(layer2.radius_ext, 0.001, 12, 15.0, copper_props)
        # Need to update core equivalent properties first before calculating next step
        core.gmr = gmr_eq_calc # Update core GMR based on previous step
        core.cross_section += layer2.cross_section # Update core area
        push!(core.layers, layer2) # Add layer for subsequent GMD calculation

        beta2 = core.cross_section / (core.cross_section + layer3_wa.cross_section)
        gmd_core_layer3 = calc_gmd(core.layers[end], layer3_wa) # GMD between tubular layer2 and wire array layer3

        gmr_eq2_expected =
            (core.gmr^(beta2^2)) * (layer3_wa.gmr^((1 - beta2)^2)) *
            (gmd_core_layer3^(2 * beta2 * (1 - beta2)))
        gmr_eq2_calc = calc_equivalent_gmr(core, layer3_wa)
        @test gmr_eq2_calc ≈ gmr_eq2_expected atol = TEST_TOL
    end


    @testset "Inductance calculations" begin
        # calc_tubular_inductance
        r_in, r_ext = 0.01, 0.02
        mu_r = 1.0
        L_expected = mu_r * μ₀ / (2 * π) * log(r_ext / r_in)
        @test calc_tubular_inductance(r_in, r_ext, mu_r) ≈ L_expected atol = TEST_TOL
        @test calc_tubular_inductance(r_in, r_ext, 2.0 * mu_r) ≈ 2.0 * L_expected atol =
            TEST_TOL #  Check mu_r scaling
        # Internal inductance of solid conductor is infinite in this simple model
        @test calc_tubular_inductance(0.0, r_ext, mu_r) == Inf

        # calc_inductance_trifoil - Requires benchmark data or simplified checks
        # This is complex. A simple check could be ensuring L > 0 for typical inputs.
        r_in_co, r_ext_co = 0.01, 0.015
        r_in_scr, r_ext_scr = 0.02, 0.022
        S = 0.1
        L_trifoil =
            calc_inductance_trifoil(r_in_co, r_ext_co, copper_props.rho, copper_props.mu_r,
                r_in_scr, r_ext_scr, copper_props.rho, copper_props.mu_r, S)
        # Trifoil inductance should be positive
        @test L_trifoil > 0
        # Could test sensitivity: increasing S should generally decrease L
        L_trifoil_S2 =
            calc_inductance_trifoil(r_in_co, r_ext_co, copper_props.rho, copper_props.mu_r,
                r_in_scr, r_ext_scr, copper_props.rho, copper_props.mu_r, 2 * S)
        # Increasing separation S should decrease L
        @test L_trifoil_S2 < L_trifoil

    end

    @testset "Capacitance & conductance" begin
        r_in, r_ext = 0.01, 0.02
        eps_r = insulator_props.eps_r
        rho_ins = insulator_props.rho

        # calc_shunt_capacitance
        C_expected = 2 * π * ε₀ * eps_r / log(r_ext / r_in)
        @test calc_shunt_capacitance(r_in, r_ext, eps_r) ≈ C_expected atol = TEST_TOL
        # Increasing r_ext decreases C
        @test calc_shunt_capacitance(r_in, r_ext * 10, eps_r) < C_expected
        # Decreasing r_in decreases C
        @test calc_shunt_capacitance(r_in / 10, r_ext, eps_r) < C_expected
        # If r_in -> r_ext, log -> 0, C -> Inf. Test approach?
        # Capacitance -> Inf as r_in approaches r_ext
        @test isinf(calc_shunt_capacitance(r_ext, r_ext, eps_r))

        # calc_shunt_conductance
        G_expected = 2 * π * (1 / rho_ins) / log(r_ext / r_in)
        @test calc_shunt_conductance(r_in, r_ext, rho_ins) ≈ G_expected atol = TEST_TOL
        # Lower rho increases G
        @test calc_shunt_conductance(r_in, r_ext, rho_ins / 10) ≈ 10 * G_expected atol =
            TEST_TOL
        # Infinite rho (perfect insulator) gives zero G
        @test calc_shunt_conductance(r_in, r_ext, Inf) ≈ 0.0 atol = TEST_TOL
        # Conductance -> Inf as r_in approaches r_ext
        @test isinf(calc_shunt_conductance(r_ext, r_ext, rho_ins))
    end

    @testset "Equivalent dielectric properties consistency" begin
        r_in, r_ext = 0.01, 0.02
        eps_r = insulator_props.eps_r
        rho_ins = insulator_props.rho
        C_eq = calc_shunt_capacitance(r_in, r_ext, eps_r)
        G_eq = calc_shunt_conductance(r_in, r_ext, rho_ins)

        # calc_equivalent_eps
        # Inverse check: eps_r from C_eq
        @test calc_equivalent_eps(C_eq, r_ext, r_in) ≈ eps_r atol = TEST_TOL

        # calc_sigma_lossfact & inverse check
        sigma_eq = calc_sigma_lossfact(G_eq, r_in, r_ext)
        # Check sigma_eq calculation
        @test sigma_eq ≈ 1 / rho_ins atol = TEST_TOL
        # Conductance from sigma
        G_from_sigma = 2 * π * sigma_eq / log(r_ext / r_in)
        # Inverse check: G_eq from sigma_eq
        @test G_from_sigma ≈ G_eq atol = TEST_TOL

        # calc_equivalent_lossfact
        f = 50.0
        ω = 2 * π * f
        tand_expected = G_eq / (ω * C_eq)
        @test calc_equivalent_lossfact(G_eq, C_eq, ω) ≈ tand_expected atol = TEST_TOL
    end

    @testset "Equivalent resistivity consistency" begin
        r_in, r_ext = 0.01, 0.02
        rho = copper_props.rho
        alpha = copper_props.alpha
        T0 = copper_props.T0

        R_tube = calc_tubular_resistance(r_in, r_ext, rho, alpha, T0, T0)
        rho_eq = calc_equivalent_rho(R_tube, r_ext, r_in)
        # Inverse check: rho from R_tube
        @test rho_eq ≈ rho atol = TEST_TOL

        R_solid = calc_tubular_resistance(0.0, r_ext, rho, alpha, T0, T0)
        rho_eq_solid = calc_equivalent_rho(R_solid, r_ext, 0.0)
        # Inverse check: rho from R_solid
        @test rho_eq_solid ≈ rho atol = TEST_TOL
    end

    @testset "Solenoid correction consistency" begin
        num_turns = 10.0 # turns/m
        r_con_ext = 0.01
        r_ins_ext = 0.015

        mu_r_corr = calc_solenoid_correction(num_turns, r_con_ext, r_ins_ext)
        #  Correction factor should be > 1 for non-zero turns
        @test mu_r_corr > 1.0

        # No twist (num_turns = NaN as per code comment)
        # Correction factor is 1 if num_turns is NaN
        @test calc_solenoid_correction(NaN, r_con_ext, r_ins_ext) ≈ 1.0 atol = TEST_TOL
        # Zero turns
        # Correction factor is 1 if num_turns is
        @test calc_solenoid_correction(0.0, r_con_ext, r_ins_ext) ≈ 1.0 atol = TEST_TOL

        # Edge case: r_con_ext == r_ins_ext (zero thickness insulator)
        # This leads to log(1) = 0 in denominator. Should return 1 or NaN/Inf?
        # Let's test the behavior. Assuming it might result in NaN due to 0/0 or X/0.
        # Correction factor is likely NaN if radii are equal (0/0 form)
        # Or maybe it should default to 1? Depends on desired behavior.
        # If the function should handle this, add a check inside it.
        @test isnan(calc_solenoid_correction(num_turns, r_con_ext, r_con_ext))

    end

    @testset "Basic uncertainty propagation" begin
        using Measurements
        r_in_m = (0.01 ± 0.001)
        r_ext_m = (0.02 ± 0.001)
        rho_m = (1.7241e-8 ± 0.001e-8)
        R_m = calc_tubular_resistance(r_in_m, r_ext_m, rho_m, (0.0 ± 0.0), (20.0 ± 0.0), (20.0 ± 0.0))
        @test Measurements.value(R_m) ≈ calc_tubular_resistance(
            Measurements.value(r_in_m),
            Measurements.value(r_ext_m),
            Measurements.value(rho_m),
            (0.0),
            (20.0),
            (20.0),
        ) atol =
            TEST_TOL
        @test Measurements.uncertainty(R_m) > 0
    end
end
