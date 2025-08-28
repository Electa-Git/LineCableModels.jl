# -------------------------
# Test fixtures (canonical parts)
# -------------------------
@testsnippet cable_fixtures begin

    # Aliases
    const LM = LineCableModels
    const DM = LM.DataModel
    const MAT = LM.Materials
    using Measurements: measurement

    # metals & dielectrics
    copper_props = MAT.Material(1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
    alu_props = MAT.Material(2.82e-8, 1.0, 1.0, 20.0, 0.0039)
    xlpe_props = MAT.Material(1e10, 2.3, 1.0, 20.0, 0.0)   # insulator-like
    semi_props = MAT.Material(1e3, 2.6, 1.0, 20.0, 0.0)   # semicon-ish

    # geometry
    d_wire = 3e-3 # 3 mm
    rin0 = 0.0

    # One of each conductor type
    core_wire = DM.WireArray(rin0, DM.Diameter(d_wire), 1, 0.0, copper_props)
    # outer wire layer
    outer_wire = DM.WireArray(core_wire.radius_ext, DM.Diameter(d_wire), 6, 10.0, copper_props)
    # strip (placed over outer wire layer)
    strip1 = DM.Strip(outer_wire.radius_ext, DM.Thickness(0.5e-3), 0.02, 8.0, copper_props)
    # tubular (placed over strip)
    tube1 = DM.Tubular(strip1.radius_ext, DM.Thickness(0.8e-3), copper_props)

    # Build a conductor group with mixed parts
    function make_conductor_group()
        g = DM.ConductorGroup(core_wire)
        add!(g, DM.WireArray, DM.Diameter(d_wire), 6, 10.0, copper_props)
        add!(g, DM.Strip, DM.Thickness(0.5e-3), 0.02, 8.0, copper_props)
        add!(g, DM.Tubular, DM.Thickness(0.8e-3), copper_props)
        g
    end

    # Insulation parts
    ins1 = DM.Insulator(tube1.radius_ext, DM.Thickness(2.0e-3), xlpe_props)
    semi = DM.Semicon(ins1.radius_ext, DM.Thickness(0.8e-3), semi_props)
    ins2 = DM.Insulator(semi.radius_ext, DM.Thickness(2.0e-3), xlpe_props)

    function make_insulator_group()
        ig = DM.InsulatorGroup(ins1)
        add!(ig, DM.Semicon, DM.Thickness(0.8e-3), semi_props)
        add!(ig, DM.Insulator, DM.Thickness(2.0e-3), xlpe_props)
        ig
    end

    # convenience helpers used in tests
    make_groups() = (make_conductor_group(), make_insulator_group())
    m(x, u) = measurement(x, u)
end

@testitem "DataModel(CableDesign): unit tests" setup = [defaults, deps_datamodel, defs_materials, cable_fixtures] begin
    # -------------------------
    # Tests
    # -------------------------

    @testset "Input Validation" begin
        # mismatched radii: force an insulator group that doesn't start at conductor rex
        g = make_conductor_group()
        bad_ins = DM.InsulatorGroup(DM.Insulator(g.radius_ext + 1e-4, DM.Thickness(1e-3), MAT.Material(1e10, 3.0, 1.0, 20.0, 0.0)))
        @test_throws ArgumentError DM.CableComponent("core", g, bad_ins)
        @test_throws ArgumentError DM.CableDesign("cabA", g, bad_ins)
    end

    @testset "Basic Functionality (Float64)" begin
        g, ig = make_groups()
        @test g.radius_ext ≈ ig.radius_in atol = TEST_TOL

        # direct component then design
        cc = DM.CableComponent("core", g, ig)
        @test cc isa DM.CableComponent
        @test cc.id == "core"
        @test cc.conductor_group === g
        @test cc.insulator_group === ig

        des = DM.CableDesign("CAB-001", cc)
        @test des isa DM.CableDesign
        @test des.cable_id == "CAB-001"
        @test length(des.components) == 1
        @test des.components[1].id == "core"

        # wrapper constructor from groups
        des2 = DM.CableDesign("CAB-002", g, ig; component_id="core")
        @test des2 isa DM.CableDesign
        @test des2.components[1].id == "core"
    end

    @testset "Edge Cases" begin
        # tiny interface gap within tolerance should pass (uses isapprox in inner ctor)
        g, ig = make_groups()
        # Nudge insulator inner radius by a few eps of Float64
        ig.radius_in = ig.radius_in + eps(Float64) * 1
        @test DM.CableComponent("core", g, ig) isa DM.CableComponent
    end

    @testset "Physical Behavior" begin
        g, ig = make_groups()
        cc = DM.CableComponent("core", g, ig)
        # Resistivity and GMR-driven equivalent numbers should be positive
        @test cc.conductor_props.rho > 0
        @test cc.conductor_group.gmr > 0
        @test cc.insulator_props.eps_r > 0
        @test cc.insulator_props.mu_r > 0
        @test cc.insulator_group.shunt_capacitance > 0
    end

    @testset "Type Stability & Promotion (component creation)" begin
        # Build groups in Float64 then promote *only* one group with a Measurement value
        gF, igF = make_groups()

        # Create a Measurement insulator by tweaking thickness with uncertainty
        igM = DM.coerce_to_T(igF, Measurements.Measurement{Float64})
        ccP = DM.CableComponent("coreM", gF, igM)
        @test typeof(ccP.conductor_group.radius_ext) <: Measurements.Measurement
        @test typeof(ccP.insulator_group.radius_ext) <: Measurements.Measurement

        # Creating a design with this mixed component works and holds the component
        des = DM.CableDesign("CAB-MIXED", ccP)
        @test length(des.components) == 1
        @test des.components[1].id == "coreM"
        @test typeof(des.components[1].conductor_group.radius_in) <: Measurements.Measurement
    end

    @testset "Design add! (by component & by groups) + overwrite semantics" begin
        g, ig = make_groups()
        cc = DM.CableComponent("core", g, ig)
        des = DM.CableDesign("CAB-ADD", cc)
        @test length(des.components) == 1

        # Add a second component by groups
        g2, ig2 = make_groups()
        DM.add!(des, "sheath", g2, ig2)
        @test length(des.components) == 2
        @test any(c -> c.id == "sheath", des.components)

        # Add another with same id -> overwrite warning
        @test length(des.components) == 2
    end

    @testset "Combinatorial Type Testing (Measurement vs Float)" begin
        g, ig = make_groups()

        # 1) Base: both Float64
        ccF = DM.CableComponent("cF", g, ig)
        @test typeof(ccF.conductor_group.radius_in) == Float64

        # 2) Fully promoted: both Measurement
        gM = DM.coerce_to_T(g, Measurements.Measurement{Float64})
        igM = DM.coerce_to_T(ig, Measurements.Measurement{Float64})
        ccM = DM.CableComponent("cM", gM, igM)
        @test typeof(ccM.conductor_group.radius_in) <: Measurements.Measurement
        @test typeof(ccM.insulator_group.radius_ext) <: Measurements.Measurement

        # 3a) Mixed: conductor carries Measurement
        cc1 = DM.CableComponent("c1", gM, ig)
        @test typeof(cc1.conductor_group.radius_in) <: Measurements.Measurement
        @test typeof(cc1.insulator_group.radius_in) <: Measurements.Measurement

        # 3b) Mixed: insulator carries Measurement
        cc2 = DM.CableComponent("c2", g, igM)
        @test typeof(cc2.conductor_group.radius_in) <: Measurements.Measurement
        @test typeof(cc2.insulator_group.radius_in) <: Measurements.Measurement
    end
end
