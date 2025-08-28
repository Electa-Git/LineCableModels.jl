@testsnippet defs_cablepos begin
    using Test
    using LineCableModels
    const DM = LineCableModels.DataModel
    const MAT = LineCableModels.Materials
    using Measurements

    # ---- helpers ----------------------------------------------------------

    # Minimal Float64 design with matching interface radii
    function _make_design_F64()
        mC = MAT.Material(1e-8, 1.0, 1.0, 20.0, 0.0)
        mI = MAT.Material(1e12, 2.5, 1.0, 20.0, 0.0)

        cg = DM.ConductorGroup(DM.Tubular(0.010, 0.012, mC))
        ig = DM.InsulatorGroup(DM.Insulator(0.012, 0.016, mI))

        cc = DM.CableComponent("core", cg, ig)
        return DM.CableDesign("CAB", cc)
    end

    # Outermost radius of the last component (for placement checks)
    _out_radius(des) = max(
        des.components[end].conductor_group.radius_ext,
        des.components[end].insulator_group.radius_ext,
    )
end

@testitem "DataModel(CablePosition): constructor unit tests" setup = [defaults, defs_cablepos] begin
    @testset "Basic construction (Float64)" begin
        des = _make_design_F64()
        rmax = _out_radius(des)
        pos = DM.CablePosition(des, 1.0, rmax + 0.10)   # default mapping

        @test pos isa DM.CablePosition
        @test DM.eltype(pos) == Float64
        @test pos.design_data === des                    # no promotion → same object
        @test pos.horz == 1.0
        @test pos.vert == rmax + 0.10
        @test length(pos.conn) == length(des.components)
        @test any(!iszero, pos.conn)                     # at least one non-grounded
    end

    @testset "Phase mapping (Dict-based)" begin
        des = _make_design_F64()
        rmax = _out_radius(des)

        # map by component id (unknown ids are rejected, missing ids default to 0)
        conn = Dict(des.components[1].id => 1)
        pos = DM.CablePosition(des, 0.0, rmax + 0.05, conn)

        @test pos.conn[1] == 1
        @test all(i == 1 ? pos.conn[i] == 1 : pos.conn[i] == 0 for i in 1:length(pos.conn))

        bad = Dict("does-not-exist" => 1)
        @test_throws ArgumentError DM.CablePosition(des, 0.0, rmax + 0.05, bad)
    end

    @testset "Geometry validation" begin
        des = _make_design_F64()
        rmax = _out_radius(des)

        # exactly at z=0 is forbidden
        @test_throws ArgumentError DM.CablePosition(des, 0.0, 0.0)

        # inside outer radius (crossing interface) is forbidden
        @test_throws ArgumentError DM.CablePosition(des, 0.0, rmax * 0.5)
    end

    @testset "Type stability & promotion" begin
        desF = _make_design_F64()
        rmax = _out_radius(desF)
        vertM = measurement(rmax + 0.10, 1e-6)

        posM = DM.CablePosition(desF, 0.0, vertM)

        @test DM.eltype(posM) <: Measurement
        @test DM.eltype(posM.design_data) <: Measurement   # design promoted with position
        @test posM.vert === vertM                          # identity preserved
        @test typeof(posM.horz) <: Measurement             # coerced to same scalar type
        @test DM.coerce_to_T(posM, DM.eltype(posM)) === posM
    end
end
