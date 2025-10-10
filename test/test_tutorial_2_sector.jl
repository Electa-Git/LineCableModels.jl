@testitem "examples/tutorial2_sector.jl tests" setup = [defaults] begin
    # Replicate the setup from the tutorial

    # === Materials ===
    materials = MaterialsLibrary(add_defaults=true)
    pvc = Material(Inf, 8.0, 1.0, 20.0, 0.1)
    add!(materials, "pvc", pvc)
    copper = get(materials, "copper")
    aluminum = get(materials, "aluminum")

    @testset "Material setup" begin
        @test get(materials, "pvc") isa LineCableModels.Materials.Material
        @test get(materials, "aluminum") isa LineCableModels.Materials.Material
        @test get(materials, "copper") isa LineCableModels.Materials.Material
    end

    # === Sector (core) geometry ===
    @testset "Sector core construction" begin
        n_sectors = 3
        r_back_mm = 10.24
        d_sector_mm = 9.14
        r_corner_mm = 1.02
        theta_cond_deg = 119.0
        ins_thick = 1.1e-3

        sector_params = SectorParams(
            n_sectors,
            r_back_mm / 1000,
            d_sector_mm / 1000,
            r_corner_mm / 1000,
            theta_cond_deg,
            ins_thick
        )

        rot_angles = (0.0, 120.0, 240.0)
        sectors = [Sector(sector_params, ang, aluminum) for ang in rot_angles]
        insulators = [SectorInsulator(sectors[i], ins_thick, pvc) for i in 1:3]

        @test length(sectors) == 3
        @test all(s -> s isa Sector, sectors)
        @test length(insulators) == 3
        @test all(i -> i isa SectorInsulator, insulators)

        components = [
            CableComponent("core1", ConductorGroup(sectors[1]), InsulatorGroup(insulators[1])),
            CableComponent("core2", ConductorGroup(sectors[2]), InsulatorGroup(insulators[2])),
            CableComponent("core3", ConductorGroup(sectors[3]), InsulatorGroup(insulators[3]))
        ]
        @test length(components) == 3
        @test components[1].id == "core1"
    end

    # === Concentric neutral ===
    @testset "Concentric neutral construction" begin
        n_neutral = 30
        r_strand = 0.79e-3
        R_N = 14.36e-3
        R_O = 17.25e-3

        inner_radius_neutral = R_N - r_strand
        outer_jacket_thickness = R_O - (R_N + r_strand)

        neutral_wires = WireArray(
            inner_radius_neutral,
            Diameter(2*r_strand),
            n_neutral,
            0.0,
            copper
        )
        @test neutral_wires isa WireArray

        neutral_jacket = Insulator(neutral_wires, Thickness(outer_jacket_thickness), pvc)
        @test neutral_jacket isa Insulator

        neutral_component = CableComponent("neutral", ConductorGroup(neutral_wires), InsulatorGroup(neutral_jacket))
        @test neutral_component.id == "neutral"
    end

    # === Assemble cable design ===
    @testset "Full cable design assembly" begin
        # Re-create components for this testset to be self-contained
        n_sectors = 3
        r_back_mm = 10.24
        d_sector_mm = 9.14
        r_corner_mm = 1.02
        theta_cond_deg = 119.0
        ins_thick = 1.1e-3
        sector_params = SectorParams(n_sectors, r_back_mm/1000, d_sector_mm/1000, r_corner_mm/1000, theta_cond_deg, ins_thick)
        rot_angles = (0.0, 120.0, 240.0)
        sectors = [Sector(sector_params, ang, aluminum) for ang in rot_angles]
        insulators = [SectorInsulator(sectors[i], ins_thick, pvc) for i in 1:3]
        components = [
            CableComponent("core1", ConductorGroup(sectors[1]), InsulatorGroup(insulators[1])),
            CableComponent("core2", ConductorGroup(sectors[2]), InsulatorGroup(insulators[2])),
            CableComponent("core3", ConductorGroup(sectors[3]), InsulatorGroup(insulators[3]))
        ]

        n_neutral = 30
        r_strand = 0.79e-3
        R_N = 14.36e-3
        R_O = 17.25e-3
        inner_radius_neutral = R_N - r_strand
        outer_jacket_thickness = R_O - (R_N + r_strand)
        neutral_wires = WireArray(inner_radius_neutral, Diameter(2*r_strand), n_neutral, 0.0, copper)
        neutral_jacket = Insulator(neutral_wires, Thickness(outer_jacket_thickness), pvc)
        neutral_component = CableComponent("neutral", ConductorGroup(neutral_wires), InsulatorGroup(neutral_jacket))

        design = CableDesign("NAYCWY_O_3x95_30x2_5", components[1])
        add!(design, components[2])
        add!(design, components[3])
        add!(design, neutral_component)

        @test length(design.components) == 4
        @test design.cable_id == "NAYCWY_O_3x95_30x2_5"
        @test design.components[1].id == "core1"
        @test design.components[2].id == "core2"
        @test design.components[3].id == "core3"
        @test design.components[4].id == "neutral"
    end

    @testset "DataFrame and preview" begin
        # Re-create the full design
        n_sectors = 3
        r_back_mm = 10.24
        d_sector_mm = 9.14
        r_corner_mm = 1.02
        theta_cond_deg = 119.0
        ins_thick = 1.1e-3
        sector_params = SectorParams(n_sectors, r_back_mm/1000, d_sector_mm/1000, r_corner_mm/1000, theta_cond_deg, ins_thick)
        rot_angles = (0.0, 120.0, 240.0)
        sectors = [Sector(sector_params, ang, aluminum) for ang in rot_angles]
        insulators = [SectorInsulator(sectors[i], ins_thick, pvc) for i in 1:3]
        components = [
            CableComponent("core1", ConductorGroup(sectors[1]), InsulatorGroup(insulators[1])),
            CableComponent("core2", ConductorGroup(sectors[2]), InsulatorGroup(insulators[2])),
            CableComponent("core3", ConductorGroup(sectors[3]), InsulatorGroup(insulators[3]))
        ]
        n_neutral = 30
        r_strand = 0.79e-3
        R_N = 14.36e-3
        R_O = 17.25e-3
        inner_radius_neutral = R_N - r_strand
        outer_jacket_thickness = R_O - (R_N + r_strand)
        neutral_wires = WireArray(inner_radius_neutral, Diameter(2*r_strand), n_neutral, 0.0, copper)
        neutral_jacket = Insulator(neutral_wires, Thickness(outer_jacket_thickness), pvc)
        neutral_component = CableComponent("neutral", ConductorGroup(neutral_wires), InsulatorGroup(neutral_jacket))
        design = CableDesign("NAYCWY_O_3x95_30x2_5", components[1])
        add!(design, components[2])
        add!(design, components[3])
        add!(design, neutral_component)

        # Test that DataFrame constructors do not throw errors
        @test DataFrame(design, :detailed) isa DataFrame
        @test DataFrame(design, :components) isa DataFrame
        @test DataFrame(design, :baseparams) isa DataFrame

        # Test that preview functions execute without error
        @test preview(design, display_plot=false) isa Any
    end

    @testset "Error handling" begin
        # Test invalid geometric parameters for Sector
        @test_throws ArgumentError SectorParams(3, -10.24/1000, 9.14/1000, 1.02/1000, 119.0, 1.1e-3)
        @test_throws ArgumentError SectorParams(3, 10.24/1000, -9.14/1000, 1.02/1000, 119.0, 1.1e-3)
        @test_throws ArgumentError SectorParams(3, 10.24/1000, 9.14/1000, -1.02/1000, 119.0, 1.1e-3)
        @test_throws ArgumentError SectorParams(3, 10.24/1000, 9.14/1000, 1.02/1000, 119.0, -1.1e-3)
    end
end
