@testsnippet defs_con_group begin
    # Canonical geometry and helpers reused across tests
    # Materials come from `defs_materials` (e.g., `copper_props`, `materials`)
    using Measurements

    const rad_wire0 = 0.0015           # 1.5 mm wire radius
    const lay_ratio0 = 12.0            # arbitrary, > 0

    # A fresh single‑wire core (center conductor)
    make_core_group() = ConductorGroup(
        WireArray(0.0, rad_wire0, 1, 0.0, copper_props; temperature=20.0, lay_direction=1)
    )

    # Convenience: a Float64 Tubular sleeve over the given inner radius
    make_tubular_over(rin, t, mat) = Tubular(rin, Thickness(t), mat; temperature=20.0)

    # Measurement helpers
    m(x, u) = measurement(x, u)
end

@testitem "DataModel(ConductorGroup.add!): unit tests" setup = [defaults, deps_datamodel, defs_materials, defs_con_group] begin
    using Measurements

    @testset "Input Validation (wrapper triggers validate!)" begin
        g = make_core_group()

        # Missing required args for WireArray (radius_wire, num_wires, lay_ratio, material)
        @test_throws ArgumentError add!(g, WireArray)
        @test_throws ArgumentError add!(g, WireArray, rad_wire0)
        @test_throws ArgumentError add!(g, WireArray, rad_wire0, 6)
        @test_throws ArgumentError add!(g, WireArray, rad_wire0, 6, lay_ratio0)

        # Invalid types forwarded to part validator
        @test_throws ArgumentError add!(g, WireArray, "bad", 6, lay_ratio0, copper_props)
        @test_throws ArgumentError add!(g, WireArray, rad_wire0, "bad", lay_ratio0, copper_props)
        @test_throws ArgumentError add!(g, WireArray, rad_wire0, 6, "bad", copper_props)
        @test_throws ArgumentError add!(g, WireArray, rad_wire0, 6, lay_ratio0, "not_a_material")

        # Out of range / geometry violations caught by rules
        @test_throws ArgumentError add!(g, WireArray, -rad_wire0, 6, lay_ratio0, copper_props)
        @test_throws ArgumentError add!(g, WireArray, 0.0, 6, lay_ratio0, copper_props)  # radius_wire > 0
        @test_throws ArgumentError add!(g, WireArray, rad_wire0, 0, lay_ratio0, copper_props) # num_wires > 0

        # Unknown keyword should be rejected by sanitize/keyword_fields policy (if enforced upstream)
        # NOTE: If this currently passes, add rejection in `sanitize` for unknown keywords.
        @test_throws ArgumentError add!(g, WireArray, rad_wire0, 6, lay_ratio0, copper_props; not_a_kw=1)
    end

    @testset "Basic Functionality (Float64)" begin
        g = make_core_group()
        @test g isa ConductorGroup
        @test length(g.layers) == 1
        @test g.radius_in == 0.0
        @test g.radius_ext ≈ rad_wire0 atol = TEST_TOL

        # Add another wire layer using Diameter convenience + defaults (radius_in auto = g.radius_ext)
        d_w = 2 * rad_wire0
        g = add!(g, WireArray, Diameter(d_w), 6, 15.0, copper_props)  # lay_direction defaults to 1
        @test g isa ConductorGroup
        @test length(g.layers) == 2
        @test g.layers[end] isa WireArray
        @test g.layers[end].lay_direction == 1  # came from keyword_defaults for WireArray

        # Geometry stacks outward and resistance decreases (parallel)
        @test g.radius_ext > rad_wire0
        @test g.resistance < g.layers[1].resistance

        # Add an outer tubular sleeve by thickness proxy
        outer_before = g.radius_ext
        g = add!(g, Tubular, Thickness(0.002), copper_props)  # temperature default from keyword_defaults(Tubular)
        @test g.layers[end] isa Tubular
        @test g.radius_ext ≈ outer_before + 0.002 atol = TEST_TOL
    end

    @testset "Edge Cases" begin
        # Very thin sleeve
        g = make_core_group()
        outer0 = g.radius_ext
        g = add!(g, Tubular, Thickness(1e-6), copper_props)
        @test g.radius_ext ≈ outer0 + 1e-6 atol = TEST_TOL

        # WireArray with lay_ratio very small but positive
        g = make_core_group()
        g = add!(g, WireArray, rad_wire0, 3, 1e-6, copper_props)
        @test g.layers[end] isa WireArray
    end

    @testset "Physical Behavior" begin
        g = make_core_group()
        # Add two conductor layers: resistance should drop further
        R0 = g.resistance
        g = add!(g, WireArray, rad_wire0, 6, 10.0, copper_props)
        R1 = g.resistance
        g = add!(g, WireArray, rad_wire0, 12, 10.0, copper_props)
        R2 = g.resistance
        @test R1 < R0
        @test R2 < R1

        # Cross-section should be monotone increasing
        @test g.cross_section > 0
        cs = [p.cross_section for p in g.layers if p isa LineCableModels.DataModel.AbstractConductorPart]
        @test all(>(0), cs)
    end

    @testset "Type Stability & Promotion (group)" begin
        # Base: purely Float64 group
        gF = make_core_group()
        @test eltype(gF) == Float64
        @test typeof(gF.radius_ext) == Float64

        # Promote by adding a Measurement argument (e.g., temperature)
        gF_before_id = objectid(gF)
        gP = add!(gF, WireArray, rad_wire0, 6, 10.0, copper_props; temperature=m(20.0, 0.1))
        @test gP !== gF                     # returned a promoted group
        @test eltype(gP) <: Measurement
        @test typeof(gP.radius_ext) <: Measurement
        # Original left intact
        @test objectid(gF) == gF_before_id
        @test length(gF.layers) == 1

        # In‑place when already Measurement
        gM = LineCableModels.DataModel.coerce_to_T(make_core_group(), Measurement{Float64})
        id_before = objectid(gM)
        gM2 = add!(gM, WireArray, m(rad_wire0, 1e-6), 6, 10.0, copper_props)
        @test gM2 === gM                    # mutated in place
        @test objectid(gM) == id_before
        @test eltype(gM) <: Measurement
    end

    @testset "Combinatorial Type Testing (constructor path of added part)" begin
        # All Float64
        g = make_core_group()
        g1 = add!(g, WireArray, rad_wire0, 6, 10.0, copper_props)
        @test eltype(g1) == Float64

        # All Measurement (radius_wire, lay_ratio, temperature)
        g = make_core_group()
        g2 = add!(g, WireArray, m(rad_wire0, 1e-6), 6, m(10.0, 0.1), copper_props; temperature=m(20.0, 0.1))
        @test eltype(g2) <: Measurement

        # Mixed case A: first numeric arg is Measurement
        g = make_core_group()
        g3 = add!(g, WireArray, m(rad_wire0, 1e-6), 6, 10.0, copper_props)
        @test eltype(g3) <: Measurement

        # Mixed case B: middle arg (lay_ratio) is Measurement
        g = make_core_group()
        g4 = add!(g, WireArray, rad_wire0, 6, m(10.0, 0.1), copper_props)
        @test eltype(g4) <: Measurement
    end
end
