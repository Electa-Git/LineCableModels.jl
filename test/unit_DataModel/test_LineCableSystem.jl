@testsnippet defs_linesys begin
    using Test
    using LineCableModels
    const DM = LineCableModels.DataModel
    const MAT = LineCableModels.Materials
    using Measurements

    # --- helpers ----------------------------------------------------------------

    # Minimal Float64 design with matching interface radii
    function _make_design_F64(; id="CAB")
        mC = MAT.Material(1e-8, 1.0, 1.0, 20.0, 0.0)
        mI = MAT.Material(1e12, 2.5, 1.0, 20.0, 0.0)

        cg = DM.ConductorGroup(DM.Tubular(0.010, 0.012, mC))
        ig = DM.InsulatorGroup(DM.Insulator(0.012, 0.016, mI))

        cc = DM.CableComponent("core", cg, ig)
        return DM.CableDesign(id, cc)
    end

    # Promote a design to Measurement{Float64}
    function _make_design_M(; id="CABM")
        des = _make_design_F64(; id)
        return DM.coerce_to_T(des, Measurement{Float64})
    end

    # Outermost radius of the last component (for placement checks)
    _out_radius(des) = max(
        des.components[end].conductor_group.radius_ext,
        des.components[end].insulator_group.radius_ext,
    )

    # Position with explicit mapping (phase 1 by default)
    function _make_position_F64(des; phase::Int=1)
        rmax = _out_radius(des)
        conn = Dict(des.components[1].id => phase)
        return DM.CablePosition(des, 0.0, rmax + 0.20, conn)
    end

    # Measurement position (promotes design through the constructor)
    function _make_position_M(desF; phase::Int=1)
        rmax = _out_radius(desF)
        vertM = measurement(rmax + 0.25, 1e-6)
        conn = Dict(desF.components[1].id => phase)
        return DM.CablePosition(desF, 0.0, vertM, conn)
    end
end

@testitem "DataModel(LineCableSystem): constructor & add! unit tests" setup = [defaults, defs_linesys] begin
    # 1) Basic construction from CablePosition (Float64)
    @testset "Basic construction (from CablePosition)" begin
        des = _make_design_F64()
        posF = _make_position_F64(des; phase=1)
        sys = DM.LineCableSystem("SYS", 1000.0, posF)

        @test sys isa DM.LineCableSystem
        @test DM.eltype(sys) == Float64
        @test sys.line_length == 1000.0
        @test sys.num_cables == 1
        @test sys.num_phases == 1
        @test sys.cables[1] === posF                      # identity preserved
        @test DM.coerce_to_T(sys, Float64) === sys        # no-op coercion
    end

    # 2) Loose constructor from CableDesign + coordinates
    @testset "Loose constructor (from CableDesign + coords)" begin
        des = _make_design_F64()
        rmax = _out_radius(des)
        conn = Dict(des.components[1].id => 1)
        sys2 = DM.LineCableSystem("SYS2", 500.0, des, 0.10, rmax + 0.30, conn)

        @test DM.eltype(sys2) == Float64
        @test sys2.num_cables == 1
        @test sys2.num_phases == 1
        @test sys2.cables[1].design_data === des          # built via position constructor
        @test typeof(sys2.cables[1].horz) == Float64
        @test typeof(sys2.cables[1].vert) == Float64
    end

    # 3) Phase counting across multiple positions
    @testset "Phase accounting" begin
        des = _make_design_F64()
        pos1 = _make_position_F64(des; phase=1)
        sys = DM.LineCableSystem("SYS-PH", 100.0, pos1)

        # Add a second cable mapped to phase 2 at a non-overlapping position
        des2 = _make_design_F64(; id="CAB2")
        r1 = _out_radius(des)
        r2 = _out_radius(des2)
        dx = r1 + r2 + 0.05                  # strictly beyond contact
        y = max(r1, r2) + 0.20
        conn2 = Dict(des2.components[1].id => 2)

        pos2 = DM.CablePosition(des2, dx, y, conn2)
        sys = DM.add!(sys, pos2)

        @test sys.num_cables == 2
        @test sys.num_phases == 2
        @test all(p -> any(x -> x > 0, p.conn), sys.cables)
    end

    # 4) Promotion on add! (Float64 system + Measurement position → promoted system)
    @testset "Promotion on add! (Float64 → Measurement)" begin
        desF = _make_design_F64()
        posF = _make_position_F64(desF)
        sysF = DM.LineCableSystem("SYS-PR", 200.0, posF)

        # Measurement position created from a Float64 design (promotes inside)
        posM = _make_position_M(desF)
        sysP = DM.add!(sysF, posM)                       # returns promoted system

        @test DM.eltype(sysF) == Float64
        @test DM.eltype(sysP) <: Measurement
        @test sysP !== sysF
        @test sysP.num_cables == 2
        @test typeof(sysP.cables[1].horz) <: Measurement  # existing coerced during promotion
        @test typeof(sysP.cables[end].vert) <: Measurement
    end

    # 5) No-op add! when already Measurement (system mutates in place)
    @testset "No-op add! when types match (Measurement system)" begin
        desF = _make_design_F64()
        posM0 = _make_position_M(desF)                    # Measurement position
        sysM0 = DM.LineCableSystem("SYS-M", measurement(1000.0, 1e-6), posM0)

        # Add a Float64 position → coerced to Measurement; system should mutate in place
        des2 = _make_design_F64(; id="CAB-F2")
        posF2 = _make_position_F64(des2)
        sysM1 = DM.add!(sysM0, posF2)

        @test sysM1 === sysM0
        @test DM.eltype(sysM0) <: Measurement
        @test sysM0.num_cables == 2
        @test typeof(sysM0.cables[end].horz) <: Measurement
        @test typeof(sysM0.line_length) <: Measurement
    end

    # 6) Combinatorial type testing (length × cable position)
    @testset "Combinatorial type testing (constructors)" begin
        desF = _make_design_F64()
        posF = _make_position_F64(desF)
        posM = _make_position_M(desF)

        lengths = (
            250.0,
            measurement(250.0, 1e-6),
        )

        positions = (
            posF,
            posM,
        )

        for L in lengths, p in positions
            sys = DM.LineCableSystem("SYS-COMB", L, p)
            if (L isa Measurement) || (DM.eltype(p) <: Measurement)
                @test DM.eltype(sys) <: Measurement
            else
                @test DM.eltype(sys) == Float64
            end
            # Round-trip no-op coercion at current T
            @test DM.coerce_to_T(sys, DM.eltype(sys)) === sys
        end
    end
end

@testitem "DataModel(LineCableSystem): promotion safety (intern-proof)" setup = [defaults, defs_linesys] begin
    using Test
    using Measurements

    # local helper: make N Float64 positions spaced along x
    function _many_positions(des, N::Int)
        rmax = _out_radius(des)
        conn = Dict(des.components[1].id => 1)
        [DM.CablePosition(des, 0.1 * i, rmax + 0.20, conn) for i in 1:N]
    end

    @testset "Promote whole system when a single Measurement cable is added" begin
        # Build a big Float64 system (N deterministic cables)
        N = 200
        desF = _make_design_F64()
        possF = _many_positions(desF, N)

        # Build system by adding positions incrementally
        sysF = DM.LineCableSystem("SYS-BIG", 1000.0, possF[1])
        for i in 2:N
            sysF = DM.add!(sysF, possF[i])  # no promotion; mutates in place
        end

        @test DM.eltype(sysF) == Float64
        @test sysF.num_cables == N
        @test all(p -> DM.eltype(p) == Float64, sysF.cables)

        # Now the intern adds ONE measurement-typed cable position
        posM = _make_position_M(desF)  # promotes design inside position ctor

        # Expect: add! returns a promoted system, original unchanged
        local sysP
        @test_logs (:warn, r"promoted") (sysP = DM.add!(sysF, posM))

        @test sysP !== sysF
        @test DM.eltype(sysP) <: Measurement

        # Original remains Float64 and unchanged
        @test DM.eltype(sysF) == Float64
        @test sysF.num_cables == N

        # New system is Measurement everywhere
        @test DM.eltype(sysP) <: Measurement
        @test sysP !== sysF
        @test sysP.num_cables == N + 1
        @test all(p -> DM.eltype(p) <: Measurement, sysP.cables)
        @test typeof(sysP.line_length) <: Measurement

        # The new position inside the promoted system is exactly the object we added
        @test sysP.cables[end] === posM
        @test DM.eltype(sysP.cables[end].design_data) <: Measurement

        # Existing positions were coerced during promotion
        @test typeof(sysP.cables[1].vert) <: Measurement
        @test DM.eltype(sysP.cables[1].design_data) <: Measurement
    end

    @testset "No-op add! when already Measurement" begin
        # Start with a Measurement system
        desF = _make_design_F64()
        posM0 = _make_position_M(desF)
        sysM = DM.LineCableSystem("SYS-M", measurement(500.0, 1e-6), posM0)

        # Add a Float64 position → it should be coerced to Measurement and mutate in place
        posF1 = _make_position_F64(_make_design_F64())
        sysM2 = DM.add!(sysM, posF1)

        @test sysM2 === sysM
        @test DM.eltype(sysM) <: Measurement
        @test sysM.num_cables == 2
        @test typeof(sysM.cables[end].horz) <: Measurement
    end
end

