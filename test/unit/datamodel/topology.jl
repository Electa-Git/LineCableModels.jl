# The topology rule of `flatten`: coaxial groups, pipe-type cables and containers, read
# from which conductors enclose which. Several rejected designs come from the gallery of
# the data model documentation.

@testsnippet TopologyGallery begin
    copper = Material(kind = :conductor, rho = 1.7241e-8, mu_r = 1.0)
    xlpe = Material(kind = :insulator, rho = 1.0e14, eps_r = 2.3)
    air = Material(kind = :insulator, rho = Inf, eps_r = 1.0)
    jacket_material = Material(kind = :insulator, rho = 1.0e13, eps_r = 2.8)
    pe = Material(kind = :insulator, rho = 1.0e15, eps_r = 2.3)
    coaxial_member = @terminal :phase begin
        core(copper; r = 1.6e-3)
        insulation(xlpe; t = 1.2e-3)
    end
    coaxial_trefoil = cores(coaxial_member; n = 3, r = 3.5e-3, names = (:a, :b, :c),
        φ0 = -pi / 2)
    # A screened cable with an 8.5 mm outer radius.
    screened_cable() = Stack(
        terminal(:core, core(copper; r = 4.0e-3), insulation(xlpe; t = 3.4e-3)),
        terminal(:screen, sheath(copper; t = 0.1e-3)),
        jacket(pe; t = 1.0e-3))
end

@testitem "DataModel / topology / groups in a shared insulating container are rejected" tags=[:unit, :parametric] setup=[TestFixtures, TopologyGallery] begin
    const DM = LineCableModels.DataModel
    shared = "coaxial groups off the center of their container have no coaxial reduction yet"
    # flatten rejects groups that share a container, screened or not.
    screened = TestFixtures.screened_three_core_design()
    @test_throws shared DM.flatten(screened)
    @test_throws "(:a_core, :a_screen), (:b_core, :b_screen), (:c_core, :c_screen) lie " *
                 "in :pipe_fill" DM.flatten(screened)
    @test_throws shared homogenize(screened)

    coaxial_pipe = @cable "coaxial-trefoil-pipe" begin
        @pipe shape=Disk(7.0e-3) fill=air wall=jacket(jacket_material; t = 1.0e-3) begin
            coaxial_trefoil
        end
    end
    @test_throws "(:a), (:b), (:c) lie in :pipe_fill" DM.flatten(coaxial_pipe)

    sector_member = @terminal :segment begin
        solid(copper, Sector(span = 2pi / 3, r_base = 0.60e-3, r_back = 5.0e-3,
            fillet = 0.20e-3); tag = :solid_sector)
        insulation(xlpe; t = 0.25e-3)
    end
    sector_cable = @cable "solid-three-core-sector" begin
        @pipe shape=Disk(5.5e-3) fill=air wall=insulation(xlpe; t = 0.8e-3) begin
            cores(sector_member; n = 3, r = 0, names = (:phase_1, :phase_2, :phase_3))
        end
    end
    @test_throws shared DM.flatten(sector_cable)
end

@testitem "DataModel / topology / a conductor around groups off its center is pipe-type" tags=[:unit, :parametric] setup=[TopologyGallery] begin
    const DM = LineCableModels.DataModel
    metal_wall = @terminal :pipe begin
        sheath(copper; t = 0.3e-3)
        jacket(jacket_material; t = 1.0e-3)
    end
    metal_pipe = @cable "coaxial-trefoil-metal-pipe" begin
        @pipe shape=Disk(7.0e-3) fill=air wall=metal_wall begin
            coaxial_trefoil
        end
    end
    @test_throws "pipe-type cables are not supported: the conductor :pipe encloses the " *
                 "coaxial groups (:a), (:b), (:c)" DM.flatten(metal_pipe)
    @test_throws "pipe-type cables are not supported" homogenize(metal_pipe)

    # In the gallery's umbilical subcable, a copper sheath surrounds two cores and two
    # fillers.
    function umbilical_core(name)
        @terminal name begin
            core(copper; r = 1.5e-3)
            insulation(pe; t = 0.65e-3)
        end
    end
    pair_wall = @terminal :pair_sheath begin
        sheath(copper; t = 0.10e-3)
        jacket(pe; t = 0.35e-3)
    end
    pair_filler = filler(pe, Disk(1.45e-3); tag = :pe_filler)
    umbilical = @cable "umbilical-subcable" begin
        @pipe shape=Disk(4.5e-3) fill=jacket_material wall=pair_wall begin
            @at umbilical_core(:pair_a) (-2.25e-3, 0.0)
            @at umbilical_core(:pair_b) (2.25e-3, 0.0)
            @at pair_filler (0.0, -2.85e-3)
            @at pair_filler (0.0, 2.85e-3)
        end
    end
    @test_throws "the conductor :pair_sheath encloses the coaxial groups (:pair_a), " *
                 "(:pair_b)" DM.flatten(umbilical)
end

@testitem "DataModel / topology / one container per group, centered and circular" tags=[:unit, :parametric] setup=[TestFixtures, TopologyGallery] begin
    const DM = LineCableModels.DataModel
    # A conduit centered on its cable adds two annuli. The air fill runs from the cable's
    # 8.5 mm jacket to the 12 mm bore, and the conduit's wall follows.
    centered = @cable "centered-conduit" begin
        @pipe shape=Disk(12.0e-3) fill=air wall=jacket(pe; t = 2.0e-3) begin
            @at screened_cable() (0.0, 0.0)
        end
    end
    blueprint = DM.flatten(centered)
    @test blueprint.assembly_ranges == [1:2]
    layers = blueprint[2].layers
    @test all(zip(layers, ((7.5e-3, 8.5e-3), (8.5e-3, 12.0e-3), (12.0e-3, 14.0e-3)))) do (
            layer, radii)
        layer.r_in ≈ radii[1] && layer.r_ex ≈ radii[2]
    end
    @test layers[2].material == air

    # The open wire screen has no jacket, and the fill keeps its wires as holes. The
    # annular part of the fill still joins the group.
    wire_screened = Stack(
        terminal(:core, core(copper; r = 4.0e-3), insulation(xlpe; t = 3.4e-3)),
        terminal(:screen, wires(copper; shape = Disk(0.5e-3), n = 16, r = 7.9e-3)))
    open_screen = @cable "wire-screen-conduit" begin
        @pipe shape=Disk(12.0e-3) fill=air wall=jacket(pe; t = 2.0e-3) begin
            wire_screened
        end
    end
    @test any(region -> region.primitive isa DM.DifferenceShape,
        open_screen.geometry.regions)
    layers = DM.flatten(open_screen).dielectrics
    @test all(zip(layers, ((4.0e-3, 7.4e-3), (8.4e-3, 12.0e-3), (12.0e-3, 14.0e-3)))) do (
            layer, radii)
        layer.r_in ≈ radii[1] && layer.r_ex ≈ radii[2]
    end

    # A cable on the bottom of its conduit leaves a crescent of fill, off its center.
    offset = TestFixtures.offset_conduit_design()
    crescent = only(filter(region -> region.primitive isa DM.DifferenceShape,
        offset.geometry.regions))
    @test area(crescent.primitive) ≈ pi * (12.0e-3^2 - 8.5e-3^2)
    @test_throws "coaxial groups off the center of their container have no coaxial " *
                 "reduction yet: (:core, :screen) lie in :pipe_fill" DM.flatten(offset)

    elliptical = @cable "elliptical-duct" begin
        @duct shape=Ellipse(14.0e-3, 10.0e-3) fill=air wall=jacket(pe; t = 2.0e-3) begin
            @at screened_cable() (0.0, 0.0)
        end
    end
    @test_throws "a non-circular container has no coaxial reduction: the coaxial groups " *
                 "(:core, :screen) lie in :duct_fill" DM.flatten(elliptical)
end
