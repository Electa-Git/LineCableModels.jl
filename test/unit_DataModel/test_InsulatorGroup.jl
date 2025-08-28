@testsnippet defs_ins_group begin
    using Measurements
    # Canonical dielectric material for tests
    const ins_props = Material(1e10, 3.0, 1.0, 20.0, 0.0)

    # Fresh inner insulator (Float64)
    make_ins_group() = InsulatorGroup(
        Insulator(0.02, 0.025, ins_props; temperature=20.0)
    )

    # Measurement helper
    m(x, u) = measurement(x, u)
end

@testitem "DataModel(InsulatorGroup.add!): unit tests" setup = [defaults, deps_datamodel, defs_materials, defs_ins_group] begin
    using Measurements

    @testset "Input Validation (wrapper triggers validate!)" begin
        g = make_ins_group()

        # Missing required args for Insulator: (radius_in provided by wrapper), need radius_ext, material
        @test_throws ArgumentError add!(g, Insulator)
        @test_throws ArgumentError add!(g, Insulator, 0.03)

        # Invalid types
        @test_throws ArgumentError add!(g, Insulator, "bad", ins_props)
        @test_throws ArgumentError add!(g, Insulator, 0.03, "not_a_material")

        # Geometry violations
        @test_throws ArgumentError add!(g, Insulator, 0.0, ins_props)   # outer cannot be 0 beyond rin
    end

    @testset "Basic Functionality (Float64)" begin
        g = make_ins_group()
        @test g isa InsulatorGroup
        @test length(g.layers) == 1
        @test g.radius_in ≈ 0.02 atol = TEST_TOL
        @test g.radius_ext ≈ 0.025 atol = TEST_TOL

        # Add a Semicon by thickness proxy (outer radius = rin + t). radius_in defaults to group.radius_ext
        t = 0.002
        rin_before = g.radius_ext
        g = add!(g, Semicon, Thickness(t), ins_props; f=60.0)
        @test g.layers[end] isa Semicon
        @test g.radius_ext ≈ rin_before + t atol = TEST_TOL
    end

    @testset "Edge Cases" begin
        g = make_ins_group()
        tsmall = 1e-6
        re0 = g.radius_ext
        g = add!(g, Semicon, Thickness(tsmall), ins_props; f=60.0)
        @test g.radius_ext ≈ re0 + tsmall atol = TEST_TOL
    end

    @testset "Physical Behavior (admittance parallel update)" begin
        g = make_ins_group()
        # Capture before
        C0 = g.shunt_capacitance
        G0 = g.shunt_conductance

        # Add another dielectric shell; admittances should combine → values typically decrease
        g = add!(g, Insulator, 0.03, ins_props; f=60.0)
        @test g.shunt_capacitance <= C0   # decreasing typical
        @test g.shunt_conductance <= G0   # decreasing typical
    end

    @testset "Type Stability & Promotion (group)" begin
        # Base Float64 group
        gF = make_ins_group()
        @test eltype(gF) == Float64
        @test typeof(gF.radius_ext) == Float64

        # Promote by Measurement temperature in part defaults
        gF_before = objectid(gF)
        gP = add!(gF, Insulator, 0.03, ins_props; f=60.0, temperature=m(20.0, 0.2))
        @test gP !== gF
        @test eltype(gP) <: Measurement
        @test typeof(gP.radius_ext) <: Measurement
        @test objectid(gF) == gF_before
        @test length(gF.layers) == 1

        # Already Measurement → in place
        gM = LineCableModels.DataModel.coerce_to_T(make_ins_group(), Measurement{Float64})
        id0 = objectid(gM)
        gM2 = add!(gM, Semicon, Thickness(m(0.001, 1e-6)), ins_props; f=60.0)
        @test gM2 === gM
        @test objectid(gM) == id0
        @test eltype(gM) <: Measurement
    end

    @testset "Combinatorial Type Testing" begin
        # All Float64
        g = make_ins_group()
        g1 = add!(g, Insulator, 0.03, ins_props; f=60.0)
        @test eltype(g1) == Float64

        # All Measurement
        g = make_ins_group()
        g2 = add!(g, Insulator, m(0.03, 1e-6), ins_props; f=60.0, temperature=m(20.0, 0.1))
        @test eltype(g2) <: Measurement

        # Mixed A: radius_ext is Measurement
        g = make_ins_group()
        g3 = add!(g, Insulator, m(0.03, 1e-6), ins_props; f=60.0)
        @test eltype(g3) <: Measurement

        # Mixed B: pass Measurement frequency (promotes group by wrapper decision)
        g = make_ins_group()
        g4 = add!(g, Insulator, 0.03, ins_props; f=m(60.0, 0.5))
        @test eltype(g4) <: Measurement
    end
end
