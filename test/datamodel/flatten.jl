using LineCableModels.DataModel: gmr, EmptyBoundary
using LineCableModels.Commons: tolerance

@testset "DataModel flatten" begin
    flatten = LineCableModels.flatten
    copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
    steel = Material(:conductor, 1.7475e-7, 1.0, 36.6326, 20.0, 0.0045)
    polymer = Material(:insulator, 1.0e14, 2.3)
    air = material(kind = :insulator, rho = Inf)

    @testset "a terminal over one circle or annulus is its own equivalent" begin
        design = build(CableDesign, "thin_tube",
            Region(:bore, Circle(0.03 * (1 - 1e-8)), polymer),
            Terminal(:core, Region(:tube, Annulus(0.03 * (1 - 1e-8), 0.03), steel)),
            Region(:jacket, Offset(1e-3), polymer))
        core = flatten(design.tree.items[2])
        @test core.name == :core
        @test core.part.tag == :tube
        @test core.part.primitive === Annulus(0.03 * (1 - 1e-8), 0.03)
        @test core.part.material === steel
        wire = build(CableDesign, "wire",
            Terminal(:core, Region(:wire, Circle(2e-3, Placement(5e-3, 0)), copper)))
        @test flatten(only(wire.tree.items)).part.primitive ===
              Circle(2e-3, Placement(5e-3, 0))
    end

    @testset "a round outline gives the tube up to the farthest metal" begin
        wires = [Circle(1e-3, Placement(3e-3cos(θ), 3e-3sin(θ))) for θ in (0:5) .* (π / 3)]
        tube = flatten(Circle(5e-3), Circle(2e-3), wires)
        @test [r_in(tube), r_ex(tube)] ≈ [2e-3, 4e-3]
        solid = flatten(Circle(5e-3), EmptyBoundary(), wires)
        @test r_in(solid) == 0
        @test r_ex(solid) ≈ 4e-3
        # The tube lies about the center of its outline.
        at = Placement(20e-3, -5e-3)
        moved = flatten(Circle(5e-3, at), Circle(2e-3, at),
            [Circle(1e-3, at * wire.at) for wire in wires])
        @test [r_in(moved), r_ex(moved)] ≈ [2e-3, 4e-3]
        @test moved.at == at
    end

    @testset "any other outline is the cross-section of its merged region" begin
        bar = Rectangle(10e-3, 5e-3)
        wires = [Circle(1e-3, Placement(x, 0)) for x in (-3e-3, -1e-3, 1e-3, 3e-3)]
        @test flatten(bar, EmptyBoundary(), wires) === bar
    end

    @testset "a terminal over one region that is not round is its declaration" begin
        bar = Region(:bar, Rectangle(10e-3, 5e-3), copper)
        screen = Region(:screen, Offset(0.2e-3), copper)
        design = build(CableDesign, "bar", Terminal(:core, bar),
            Region(:insulation, Offset(1e-3), polymer), Terminal(:screen, screen))
        core, insulation, tube = design.tree.items
        @test flatten(core) == Terminal(:core, bar)
        # A layer that an offset declares has no cross-section of its own: it stays resolved.
        @test flatten(tube) === tube
        @test flatten(Terminal(:core, bar)) == Terminal(:core, bar)
    end

    @testset "a conductor in a fill that is not round merges into its outline" begin
        square = Rectangle(20.0e-3, 20.0e-3)
        design = build(CableDesign, "square",
            Terminal(:core,
                Group(Region(:fill, square, polymer), Region(:wire, Circle(5.0e-3), copper))),
            Region(:insulation, Offset(2.0e-3), polymer))
        core, insulation = flatten(design.tree).items
        @test core.part.tag == :core
        @test core.part.primitive === square
        # The square has the resistance and the geometric mean radius of the wire.
        @test core.part.material.rho / area(square) ≈ copper.rho / (π * 5.0e-3^2) rtol = 1e-12
        @test gmr(square, core.part.material.mu_r) ≈ gmr(Circle(5.0e-3), 1.0) rtol = 1e-12
        @test core.part.material.alpha ≈ copper.alpha
        @test insulation == Region(:insulation, Offset(2.0e-3), polymer)
    end

    @testset "a conductor with a tube around it merges into the outline of the tube" begin
        bar = Region(:bar, Rectangle(10e-3, 5e-3), copper)
        design = build(CableDesign, "plated",
            Terminal(:core, Stack(bar, Region(:plating, Offset(0.5e-3), steel))),
            Region(:insulation, Offset(1e-3), polymer))
        core, insulation = flatten(design.tree).items
        metal = design.regions[1:2]
        outline = core.part.primitive
        @test area(outline) ≈ sum(placed -> area(placed.shape), metal) rtol = 1e-12
        @test area(outline) / core.part.material.rho ≈
              sum(placed -> area(placed.shape) / placed.source.material.rho, metal) rtol = 1e-12
        @test insulation == Region(:insulation, Offset(1e-3), polymer)
        # The merged cable is a cable, and it is its own flattened form.
        cable = build(CableDesign, "merged", core, insulation)
        @test flatten(cable.tree) == Stack(core, insulation)
        @test area(cable.outer) ≈ area(design.outer) rtol = 1e-12
    end

    @testset "the strands of a terminal associate in parallel" begin
        courses = [Lay(Helix(LayRatio(12)),
                       Repeat(CircularPattern(n; r = 3.0e-3k),
                           Region(Symbol(:wires_, k), Circle(1.5e-3), copper)))
                   for (k, n) in enumerate((6, 12))]
        design = build(CableDesign, "stranded",
            Terminal(:core,
                Group(Region(:core_fill, Circle(7.5e-3), air),
                    Stack(Region(:core_centre, Circle(1.5e-3), copper), courses...))),
            Region(:insulation, Offset(3.0e-3), polymer))
        core = flatten(EmptyBoundary(), first(design.tree.items))
        region = core.part.part
        @test core.name == :core
        @test region.tag == :core
        @test r_in(region.primitive) == 0
        @test r_ex(region.primitive) ≈ 7.5e-3
        # The centre wire lies straight, and each of the 18 laid wires is longer by its
        # overlength.
        conductance = (1 + 18 / sqrt(1 + (π / 12)^2)) * π * 1.5e-3^2 / copper.rho
        @test region.material.rho / (π * 7.5e-3^2) ≈ inv(conductance) rtol = 1e-12
        @test region.material.alpha ≈ copper.alpha
        # Equal wires of unit permeability have the geometric mean radius of their
        # distances, with the weights of their conductances.
        wires = [placed.shape for placed in design.regions[1:19]]
        w = [1; fill(inv(sqrt(1 + (π / 12)^2)), 18)]
        w ./= sum(w)
        log_g = sum(w[i] * w[j] *
                    log(i == j ? gmr(wires[i], 1.0) :
                        hypot((centroid(wires[i]) .- centroid(wires[j]))...))
        for i in eachindex(wires), j in eachindex(wires))
        @test log(7.5e-3) - region.material.mu_r / 4 ≈ log_g
    end

    @testset "a terminal holds its conductor alone" begin
        bar = Region(:bar, Rectangle(10e-3, 5e-3), copper)
        design = build(CableDesign, "wrapped",
            Terminal(:core, Stack(bar, Region(:insulation, Offset(1e-3), polymer))))
        @test_throws ArgumentError flatten(EmptyBoundary(), only(design.tree.items))
        @test_throws ArgumentError flatten(design.tree)
    end

    @testset "a stack declares the conductor of a terminal in one part" begin
        tape = Region(:sheath_tape, Annulus(2e-3, 2.2e-3), copper)
        wires = Region(:sheath_tube, Offset(1e-3), copper)
        core = Region(:bore, Circle(2e-3), polymer)
        twice = build(CableDesign, "twice", core, Terminal(:sheath, tape),
            Terminal(:sheath, wires))
        @test_throws ArgumentError flatten(twice.tree)
        once = build(CableDesign, "once", core, Terminal(:sheath, Stack(tape, wires)))
        bore, sheath = flatten(once.tree).items
        @test sheath.name == :sheath
        @test [r_in(sheath.part.primitive), r_ex(sheath.part.primitive)] ≈ [2e-3, 3.2e-3]
    end

    @testset "a stack flattens to the tubes of its terminals and their layers" begin
        helix = Helix(LayRatio(12))
        courses = [Lay(helix,
                       Repeat(CircularPattern(n; r = 3.0e-3k),
                           Region(Symbol(:wires_, k), Circle(1.5e-3), copper)))
                   for (k, n) in enumerate((6, 12))]
        strands = Stack(Region(:core_centre, Circle(1.5e-3), copper), courses...)
        design = build(CableDesign, "loose",
            Terminal(:core, Group(Region(:core_fill, Circle(8.0e-3), air), strands)),
            Region(:insulation, Offset(3.0e-3), polymer),
            Terminal(:sheath, Region(:sheath, Offset(1.0e-3), copper)),
            Region(:jacket, Offset(2.0e-3), polymer))
        flat = flatten(design.tree)
        core, fill_layer, insulation, sheath, jacket = flat.items
        @test (core.name, fill_layer.tag, insulation.tag, sheath.name, jacket.tag) ==
              (:core, :core_fill, :insulation, :sheath, :jacket)
        # The fill beyond the wires is a layer of air, and the insulation is its annulus.
        tube = core.part.part.primitive
        @test [r_in(tube), r_ex(tube)] ≈ [0, 7.5e-3]
        @test [r_in(fill_layer.primitive), r_ex(fill_layer.primitive)] ≈ [7.5e-3, 8.0e-3]
        @test [r_in(insulation.primitive), r_ex(insulation.primitive)] ≈ [8.0e-3, 11.0e-3]
        # The lay over the region of the core gives the turns of its current. Its helix has
        # the pitch 1 / N at the outer radius of the region. The layers of the interval have
        # the materials of their declarations.
        radii = [3.0e-3, 6.0e-3]
        G = [1; fill(inv(overlength(helix, radii[1])), 6);
             fill(inv(overlength(helix, radii[2])), 12)]
        turns = [0; fill(inv(helix.dir * pitch(helix, radii[1])), 6);
                 fill(inv(helix.dir * pitch(helix, radii[2])), 12)]
        N = sum(G .* turns) / sum(G)
        @test core.part.helix.dir / pitch(core.part.helix, r_ex(tube)) ≈ N rtol = 1e-12
        @test fill_layer.material === air
        @test insulation.material === polymer
        # A tube of one annulus has no turns.
        @test sheath.part.material === copper
        @test jacket.material === polymer

        reached = build(CableDesign, "reached",
            Terminal(:core, Group(Region(:core_fill, Circle(7.5e-3), air), strands)),
            Region(:insulation, Offset(3.0e-3), polymer))
        # Wires that reach the container of their group leave its fill as a layer of no
        # width.
        core, flush, insulation = flatten(reached.tree).items
        @test flush.tag == :core_fill && flush.material === air
        @test r_in(flush.primitive) ≈ r_ex(flush.primitive) ≈ 7.5e-3

        nested = build(CableDesign, "nested",
            Stack(Terminal(:core, Region(:wire, Circle(2e-3), copper)),
                Region(:insulation, Offset(1e-3), polymer)),
            Region(:jacket, Offset(1e-3), polymer))
        @test length(flatten(nested.tree).items) == 3
    end

    @testset "a tape and its fill are one layer of the composite material" begin
        polyester = material(kind = :insulator, rho = 1.0e15, eps_r = 3.3, tan_delta = 2.0e-3)
        compound = material(kind = :insulator, rho = 1.0e13, eps_r = 2.2, tan_delta = 5.0e-4)
        design = build(CableDesign, "open_tape",
            Terminal(:core, Region(:core, Circle(10.0e-3), copper)),
            tape(polyester; w = 10.0e-3, t = 0.2e-3, lay = LayRatio(10), fill = compound),
            Region(:insulation, Offset(1.0e-3), polymer))
        core, layer, insulation = flatten(design.tree).items
        f = 10.0e-3 * 0.2e-3 / (π * (10.2e-3^2 - 10.0e-3^2))
        @test layer.tag == :tape_fill
        @test layer.material isa Material{Float64, :insulator}
        @test [r_in(layer.primitive), r_ex(layer.primitive)] ≈ [10.0e-3, 10.2e-3] rtol = 1e-12
        @test 1 / layer.material.rho ≈ f / 1.0e15 + (1 - f) / 1.0e13 rtol = 1e-12
        @test layer.material.eps_r ≈ f * 3.3 + (1 - f) * 2.2 rtol = 1e-12
        @test layer.material.eps_r * layer.material.tan_delta ≈
              f * 3.3 * 2.0e-3 + (1 - f) * 2.2 * 5.0e-4 rtol = 1e-12
        @test insulation.material === polymer
    end

    @testset "a group around a terminal is the conductor of that terminal" begin
        # A core in a fill is the core with a layer of the fill around it.
        filled = build(CableDesign, "filled",
            Group(Region(:core_insulation, Circle(16.0e-3), polymer),
                Terminal(:core, Region(:core, Circle(10.0e-3), copper))),
            Terminal(:screen, Region(:screen, Offset(1.0e-3), copper)))
        core, layer, screen = flatten(filled.tree).items
        @test core.name == :core
        @test [r_in(core.part.primitive), r_ex(core.part.primitive)] ≈ [0, 10.0e-3]
        @test core.part.material.rho ≈ copper.rho && core.part.material.mu_r ≈ copper.mu_r
        @test layer.tag == :core_insulation
        @test [r_in(layer.primitive), r_ex(layer.primitive)] ≈ [10.0e-3, 16.0e-3]
        @test screen.part.primitive === Annulus(16.0e-3, 17.0e-3)
        # A fill around a rod and the wires of a terminal gives a tube from the axis.
        courses = [Lay(Helix(LayRatio(12)),
                       Repeat(CircularPattern(n; r = 3.0e-3k),
                           Region(Symbol(:wires_, k), Circle(1.5e-3), copper)))
                   for (k, n) in enumerate((6, 12))]
        rod = build(CableDesign, "rod",
            Group(Region(:core_fill, Circle(7.5e-3), air),
                Stack(Region(:rod, Circle(1.5e-3), polymer),
                    Terminal(:core, Stack(courses...)))),
            Region(:insulation, Offset(3.0e-3), polymer))
        tube, insulation = flatten(rod.tree).items
        region = tube.part.part
        @test tube.name == :core
        @test [r_in(region.primitive), r_ex(region.primitive)] ≈ [0, 7.5e-3] rtol = 1e-12
        resistance = copper.rho * sqrt(1 + (π / 12)^2) / (18π * 1.5e-3^2)
        @test region.material.rho / (π * 7.5e-3^2) ≈ resistance rtol = 1e-12
    end

    @testset "foreach visits the regions of a resolved part from the inside out" begin
        design = build(CableDesign, "armoured",
            Terminal(:core, Region(:core, Circle(24.5e-3), copper)),
            Terminal(:armour,
                Group(Region(:armour_fill, Offset(4.0e-3), air),
                    Repeat(CircularPattern(40; r = 26.5e-3),
                        Region(:armour_wires, Circle(2.0e-3), steel)))))
        tags, areas = Symbol[], Float64[]
        foreach(design.tree) do placed
            push!(tags, placed.source.tag)
            push!(areas, area(placed.shape))
        end
        @test tags == [placed.source.tag for placed in design.regions]
        @test sum(areas) ≈ area(design.outer)
    end

    @testset "a layer that is no annulus is refused" begin
        design = build(CableDesign, "tube_split_dielectric",
            Region(:bore, Circle(1e-3), polymer),
            Terminal(:core, Region(:tube, Annulus(1e-3, 2e-3), copper)),
            Repeat(CircularPattern(2; r = Offset(0)),
                Region(:insulation, Rectangle(2.5π * 1e-3, 1e-3), polymer)))
        @test_throws ArgumentError flatten(design.tree)
    end

    @testset "a group about a center lies about that center when flattened" begin
        at = Placement(30e-3, 0)
        design = build(CableDesign, "aside",
            Group(Region(:core_insulation, Circle(16.0e-3, at), polymer),
                Terminal(:core, Region(:core, Circle(10.0e-3, at), copper))),
            Region(:jacket, Offset(2.0e-3), polymer))
        @test flatten(design.tree) == Stack(
            Terminal(:core, Region(:core, Circle(10.0e-3, at), copper)),
            Region(:core_insulation, Annulus(10.0e-3, 16.0e-3, at), polymer),
            Region(:jacket, Annulus(16.0e-3, 16.0e-3 + 2.0e-3, at), polymer))
    end

    @testset "a cable of single regions that are not round is its own flattened form" begin
        aluminum = Material(:conductor, 2.8264e-8, 1, 1, 20, 0.00429)
        sector = LineCableModels.Sector(
            0, 10.24e-3, deg2rad(119.0), 1.1e-3 / sind(60), 1.02e-3, Placement(0, 0))
        parts = (Terminal(:core, Region(:core, sector, aluminum)),
            Region(:insulation, Offset(1.1e-3), polymer),
            Terminal(:screen, Region(:screen, Offset(0.2e-3), copper)),
            Region(:jacket, Offset(2.0e-3), polymer))
        design = build(CableDesign, "screened_sector", parts...)
        # The screen is a layer that an offset declares: it stays the resolved part.
        @test flatten(design.tree) ==
              Stack(parts[1], parts[2], design.tree.items[3], parts[4])
        @test flatten(Stack(parts...)) == Stack(parts...)
    end

    @testset "an enclosure and an assembly keep their form around flattened members" begin
        offset = 2 * 28.0e-3 / sqrt(3)
        at = Placement(offset, 0.0)
        members = Tuple(Placement(0.0,
                            0.0,
                            θ) => Stack(
                            Terminal(name, Region(name, Circle(19.5e-3, at), copper)),
                            Region(:insulation, Offset(8.5e-3), polymer))
        for (name, θ) in zip((:a, :b, :c), (0.0, 2π / 3, 4π / 3)))
        design = build(CableDesign, "pipe_type",
            Enclosure(Region(:filler, Circle(62.0e-3), polymer), Assembly(members)),
            Terminal(:pipe, Region(:pipe, Offset(3.0e-3), steel)),
            Region(:jacket, Offset(5.0e-3), polymer))
        flat = flatten(design.tree)
        @test flat == Stack(
            Enclosure(Region(:filler, Circle(62.0e-3), polymer),
                Assembly(Tuple(Placement(0.0,
                                   0.0,
                                   θ) => Stack(
                                   Terminal(name, Region(name, Circle(19.5e-3, at), copper)),
                                   Region(:insulation,
                                       Annulus(19.5e-3, 19.5e-3 + 8.5e-3, at), polymer))
                for (name, θ) in zip((:a, :b, :c), (0.0, 2π / 3, 4π / 3))))),
            Terminal(:pipe, Region(:pipe, Annulus(62.0e-3, 62.0e-3 + 3.0e-3), steel)),
            Region(:jacket, Annulus(62.0e-3 + 3.0e-3, 62.0e-3 + 3.0e-3 + 5.0e-3), polymer))
        @test flatten(flat) == flat
        # A rod, an enclosure declared by an offset and a binder around it.
        wire = Stack(Terminal(:core, Region(:core, Circle(0.5e-3), copper)),
            Region(:insulation, Offset(0.25e-3), polymer))
        placements = [Placement(2.5e-3 * cos(θ), 2.5e-3 * sin(θ))
                      for θ in (0.0, 2π / 3, 4π / 3)]
        rod = build(CableDesign, "rod", Region(:rod, Circle(1.0e-3), polymer),
            Enclosure(Region(:fill, Offset(3.0e-3), air),
                Assembly(Tuple(at => wire for at in placements))),
            Region(:binder, Offset(0.2e-3), polymer))
        core = Stack(Terminal(:core, Region(:core, Circle(0.5e-3), copper)),
            Region(:insulation, Annulus(0.5e-3, 0.5e-3 + 0.25e-3), polymer))
        @test flatten(rod.tree) == Stack(Region(:rod, Circle(1.0e-3), polymer),
            Enclosure(Region(:fill, Annulus(1.0e-3, 1.0e-3 + 3.0e-3), air),
                Assembly(Tuple(at => core for at in placements))),
            Region(:binder, Annulus(1.0e-3 + 3.0e-3, 1.0e-3 + 3.0e-3 + 0.2e-3), polymer))
    end

    @testset "a terminal around an assembly begins at the circle that contains its members" begin
        core = Stack(Terminal(:core, Region(:core, Circle(10.0e-3), copper)),
            Region(:insulation, Offset(2.0e-3), polymer))
        design = build(CableDesign, "tight",
            Assembly((Placement(0.0, 0.0) => core,)),
            Terminal(:screen,
                Group(Region(:screen_fill, Offset(2.0e-3), air),
                    Repeat(CircularPattern(30; r = Offset(0)),
                        Region(:screen_wires, Circle(1.0e-3), copper)))),
            Region(:jacket, Offset(2.0e-3), polymer))
        assembly, screen, jacket = flatten(design.tree).items
        @test only(assembly.members).second == Stack(
            Terminal(:core, Region(:core, Circle(10.0e-3), copper)),
            Region(:insulation, Annulus(10.0e-3, 10.0e-3 + 2.0e-3), polymer))
        @test [r_in(screen.part.primitive), r_ex(screen.part.primitive)] ≈
              [12.0e-3, 14.0e-3] rtol = 1e-12
        @test r_in(jacket.primitive) ≈ r_ex(screen.part.primitive) rtol = 1e-12
    end

    @testset "a repeat keeps its placements around its flattened part" begin
        core = Stack(Terminal(:core, Region(:core, Circle(1.0e-3), copper)),
            Region(:core_insulation, Offset(0.5e-3), polymer))
        design = build(CableDesign, "three_cores",
            Enclosure(Region(:pipe_fill, Circle(6.5e-3), polymer),
                Repeat(CircularPattern(3; r = 3.5e-3), core)))
        enclosure = only(flatten(design.tree).items)
        @test enclosure.fill == Region(:pipe_fill, Circle(6.5e-3), polymer)
        @test enclosure.contents.pattern == only(design.tree.items).contents.pattern
        @test enclosure.contents.part == Stack(
            Terminal(:core, Region(:core, Circle(1.0e-3), copper)),
            Region(:core_insulation, Annulus(1.0e-3, 1.0e-3 + 0.5e-3), polymer))
        @test flatten(flatten(design.tree)) == flatten(design.tree)
    end

    # Each flattened reference is its own flattened form.
    @testset "a flattened reference flattens to itself" begin
        cases = [REFERENCES;
                 [case => include(joinpath(
                      @__DIR__, "..", "references", case, "problem.jl"))
                  for case in ("cable_30kv_630mm2_trefoil",
                     "cable_30kv_na2xs2y_630mm2_trefoil", "cable_132kv_630mm2_flathor",
                     "cable_18kv_1000mm2_trefoil", "cable_380kv_2000mm2_flatver")]]
        for (case, problem) in cases, design in unique(problem.system.designs)

            flat = flatten(design)
            @test flatten(flat) == flat
        end
    end

    # A flattened single core is a cable, and it flattens to itself.
    @testset "a flattened single core builds" begin
        cases = [REFERENCES;
                 [case => include(joinpath(
                      @__DIR__, "..", "references", case, "problem.jl"))
                  for case in ("cable_30kv_630mm2_trefoil",
                     "cable_30kv_na2xs2y_630mm2_trefoil", "cable_132kv_630mm2_flathor",
                     "cable_18kv_1000mm2_trefoil", "cable_380kv_2000mm2_flatver")]]
        for (case, problem) in cases, design in unique(problem.system.designs)

            flat = flatten(design)
            any(item -> item isa Enclosure, flat.items) && continue
            cable = build(CableDesign, design.id, flat.items...)
            @test flatten(cable.tree) == flat
            @test [row.name for row in cable.terminals] ==
                  [row.name for row in design.terminals]
        end
    end
end
