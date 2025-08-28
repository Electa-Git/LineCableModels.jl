@testitem "DataModel(CableComponent): unit tests" setup = [defaults, deps_datamodel, defs_materials] begin
    # Aliases (no `using` per project policy)
    const LM = LineCableModels
    const DM = LM.DataModel
    const MAT = LM.Materials

    # --- Canonical materials (fallbacks in case `defs_materials` lacks a key) ---
    copper = get(materials, "copper", MAT.Material(1.7241e-8, 1.0, 1.0, 20.0, 0.00393))
    aluminum = get(materials, "aluminum", MAT.Material(2.826e-8, 1.0, 1.0, 20.0, 0.00429))
    polyeth = get(materials, "polyethylene", MAT.Material(1e12, 2.3, 1.0, 20.0, 0.0))
    semimat = get(materials, "semicon", MAT.Material(1e3, 3.0, 1.0, 20.0, 0.0))

    # --- Helpers ---------------------------------------------------------------
    make_conductor_group_F = function ()
        # core: 1 wire at center (diameter d_w)
        d_w = 3e-3
        g = DM.ConductorGroup(DM.WireArray(0.0, DM.Diameter(d_w), 1, 0.0, aluminum))
        # add helical wire layer (defaults: radius_in = group.radius_ext)
        DM.add!(g, DM.WireArray, DM.Diameter(d_w), 6, 10.0, aluminum)
        # add thin strip layer
        DM.add!(g, DM.Strip, DM.Thickness(5e-4), 1.0e-2, 15.0, aluminum)
        # add tubular sheath (thickness)
        DM.add!(g, DM.Tubular, DM.Thickness(1e-3), aluminum)
        return g
    end

    make_insulator_group_F = function (rin::Real)
        # Build from inner radius `rin` outward
        g = DM.InsulatorGroup(DM.Insulator(rin, DM.Thickness(4e-3), polyeth))
        DM.add!(g, DM.Semicon, DM.Thickness(5e-4), semimat)
        return g
    end

    # Small helper for measurement creation
    m = x -> measurement(x, 0.1 * x + (x == 0 ? 1e-6 : 0))

    # --- Input Validation ------------------------------------------------------
    @testset "Input Validation" begin
        gC = make_conductor_group_F()
        # make insulator group that *does not* start exactly at gC.radius_ext
        bad_rin = gC.radius_ext + 1e-6
        gI_bad = DM.InsulatorGroup(DM.Insulator(bad_rin, DM.Thickness(4e-3), polyeth))
        DM.add!(gI_bad, DM.Semicon, DM.Thickness(5e-4), semimat)
        @test_throws ArgumentError DM.CableComponent("bad", gC, gI_bad)
    end

    # --- Basic Functionality (Float64 workflow) --------------------------------
    @testset "Basic Functionality (Float64)" begin
        gC = make_conductor_group_F()
        gI = make_insulator_group_F(gC.radius_ext)
        cc = DM.CableComponent("core", gC, gI)

        # Type & identity when no promotion needed
        @test eltype(cc) == Float64
        @test cc.id == "core"
        @test cc.conductor_group === gC
        @test cc.insulator_group === gI

        # Geometric continuity (nominal comparison)
        @test isapprox(cc.conductor_group.radius_ext, cc.insulator_group.radius_in; atol=TEST_TOL)

        # Physical sanity
        @test cc.conductor_props.rho > 0
        @test cc.conductor_props.mu_r > 0
        @test cc.insulator_props.eps_r > 0
        @test cc.insulator_props.rho > 0
    end

    # --- Edge Cases ------------------------------------------------------------
    @testset "Edge Cases" begin
        gC = DM.ConductorGroup(DM.WireArray(0.0, DM.Diameter(2e-3), 1, 0.0, copper))
        # hairline insulator: nearly zero thickness, but non-zero
        gI = DM.InsulatorGroup(DM.Insulator(gC.radius_ext, gC.radius_ext + 1e-6, polyeth))
        cc = DM.CableComponent("thin", gC, gI)
        @test cc.insulator_group.radius_ext > cc.conductor_group.radius_ext
        @test cc.insulator_props.eps_r > 0
    end

    # --- Physical Behavior (relationships that should hold) --------------------
    @testset "Physical Behavior" begin
        gC = make_conductor_group_F()
        gI = make_insulator_group_F(gC.radius_ext)
        cc = DM.CableComponent("phys", gC, gI)

        # Conductor alpha propagated
        @test isapprox(cc.conductor_props.alpha, gC.alpha; atol=TEST_TOL)

        gI2 = DM.InsulatorGroup(DM.Insulator(gC.radius_ext, DM.Thickness(6e-3), polyeth))
        DM.add!(gI2, DM.Semicon, DM.Thickness(5e-4), semimat)
        cc2 = DM.CableComponent("phys2", gC, gI2)
        @test cc2.insulator_group.shunt_capacitance < cc.insulator_group.shunt_capacitance
    end

    # --- Type Stability & Promotion -------------------------------------------
    @testset "Type Stability & Promotion" begin
        # Base: both Float64
        gC_F = make_conductor_group_F()
        gI_F = make_insulator_group_F(gC_F.radius_ext)
        cc_F = DM.CableComponent("F", gC_F, gI_F)
        @test eltype(cc_F) == Float64

        # Insulator as Measurement → component promotes
        gI_M = DM.coerce_to_T(gI_F, Measurement{Float64})
        cc_PM = DM.CableComponent("PM", gC_F, gI_M)
        @test eltype(cc_PM) <: Measurement
        @test eltype(cc_PM.conductor_group) <: Measurement
        @test eltype(cc_PM.insulator_group) <: Measurement
        # original groups untouched
        @test eltype(gC_F) == Float64
        @test eltype(gI_F) == Float64

        # Conductor as Measurement → component promotes
        gC_M = DM.coerce_to_T(gC_F, Measurement{Float64})
        cc_MP = DM.CableComponent("MP", gC_M, gI_F)
        @test eltype(cc_MP) <: Measurement

        # Both Measurement
        cc_MM = DM.CableComponent("MM", gC_M, gI_M)
        @test eltype(cc_MM) <: Measurement

        # Mixed raw creation using measurements inside groups
        gC_mix = DM.ConductorGroup(DM.WireArray(0.0, DM.Diameter(m(3e-3)), 1, 0.0, aluminum))
        DM.add!(gC_mix, DM.Tubular, DM.Thickness(m(5e-4)), copper)
        gI_mix = DM.InsulatorGroup(DM.Insulator(gC_mix.radius_ext, DM.Thickness(2e-3), polyeth))
        cc_mix = DM.CableComponent("mix", gC_mix, gI_mix)
        @test eltype(cc_mix) <: Measurement
    end

    # --- Combinatorial Type Testing (constructor path inside groups) -----------
    @testset "Combinatorial Type Testing" begin
        # Base floats
        gC = DM.ConductorGroup(DM.WireArray(0.0, DM.Diameter(2e-3), 1, 0.0, aluminum))
        gI = DM.InsulatorGroup(DM.Insulator(gC.radius_ext, DM.Thickness(3e-3), polyeth))

        # Case A: Float + Measurement (insulator group)
        gI_A = DM.coerce_to_T(gI, Measurement{Float64})
        cc_A = DM.CableComponent("A", gC, gI_A)
        @test eltype(cc_A) <: Measurement

        # Case B: Measurement + Float (conductor group)
        gC_B = DM.coerce_to_T(gC, Measurements.Measurement{Float64})
        cc_B = DM.CableComponent("B", gC_B, gI)
        @test eltype(cc_B) <: Measurement
    end
end
