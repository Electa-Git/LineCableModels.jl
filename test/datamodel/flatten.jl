using LineCableModels.DataModel: gmr
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
        core = flatten(design.tree.items[2]).first
        @test core.name == :core
        @test core.part.tag == :tube
        @test core.part.primitive === Annulus(0.03 * (1 - 1e-8), 0.03)
        @test core.part.material === steel
        wire = build(CableDesign, "wire",
            Terminal(:core, Region(:wire, Circle(2e-3, Placement(5e-3, 0)), copper)))
        @test flatten(only(wire.tree.items)).first.part.primitive ===
              Circle(2e-3, Placement(5e-3, 0))
    end

    @testset "a round outline gives the tube up to the farthest metal" begin
        wires = [Circle(1e-3, Placement(3e-3cos(θ), 3e-3sin(θ))) for θ in (0:5) .* (π / 3)]
        tube = flatten(Circle(5e-3), 2e-3, wires)
        @test [r_in(tube), r_ex(tube)] ≈ [2e-3, 4e-3]
        solid = flatten(Circle(5e-3), 0.0, wires)
        @test r_in(solid) == 0
        @test r_ex(solid) ≈ 4e-3
    end

    @testset "any other outline gives Ametani's cylinder" begin
        bar = Rectangle(10e-3, 5e-3)
        tube = flatten(bar, 0.0, [bar])
        @test r_ex(tube) ≈ 30e-3 / 2π
        @test area(tube) ≈ area(bar)
        # The wires in the bar leave its perimeter and take their own area.
        wires = [Circle(1e-3, Placement(x, 0)) for x in (-3e-3, -1e-3, 1e-3, 3e-3)]
        stranded = flatten(bar, 0.0, wires)
        @test r_ex(stranded) == r_ex(tube)
        @test area(stranded) ≈ 4π * 1e-3^2
        # An ellipse of equal axes fills its cylinder.
        disk = LineCableModels.Ellipse(3e-3, 3e-3)
        solid = flatten(disk, 0.0, [disk])
        @test r_in(solid) == 0
        @test r_ex(solid) ≈ 3e-3
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
        core = flatten(first(design.tree.items)).first
        @test core.name == :core
        @test core.part.tag == :core
        @test r_in(core.part.primitive) == 0
        @test r_ex(core.part.primitive) ≈ 7.5e-3
        # The centre wire lies straight, and each of the 18 laid wires is longer by its
        # overlength.
        conductance = (1 + 18 / sqrt(1 + (π / 12)^2)) * π * 1.5e-3^2 / copper.rho
        @test core.part.material.rho / (π * 7.5e-3^2) ≈ inv(conductance) rtol = 1e-12
        @test core.part.material.alpha ≈ copper.alpha
        # Equal wires of unit permeability have the geometric mean radius of their
        # distances, with the weights of their conductances.
        wires = [placed.shape for placed in design.regions[1:19]]
        w = [1; fill(inv(sqrt(1 + (π / 12)^2)), 18)]
        w ./= sum(w)
        log_g = sum(w[i] * w[j] *
                    log(i == j ? gmr(wires[i], 1.0) :
                        hypot((centroid(wires[i]) .- centroid(wires[j]))...))
        for i in eachindex(wires), j in eachindex(wires))
        @test log(7.5e-3) - core.part.material.mu_r / 4 ≈ log_g
    end

    @testset "a terminal holds its conductor alone" begin
        bar = Region(:bar, Rectangle(10e-3, 5e-3), copper)
        design = build(CableDesign, "wrapped",
            Terminal(:core, Stack(bar, Region(:insulation, Offset(1e-3), polymer))))
        @test_throws ArgumentError flatten(only(design.tree.items))
        @test_throws ArgumentError flatten(design.tree)
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
        @test [r_in(core.part.primitive), r_ex(core.part.primitive)] ≈ [0, 7.5e-3]
        @test [r_in(fill_layer.primitive), r_ex(fill_layer.primitive)] ≈ [7.5e-3, 8.0e-3]
        @test [r_in(insulation.primitive), r_ex(insulation.primitive)] ≈ [8.0e-3, 11.0e-3]
        # The turns of the core go into the permeability of the layers of its interval.
        radii = [3.0e-3, 6.0e-3]
        G = [1; fill(inv(overlength(helix, radii[1])), 6);
             fill(inv(overlength(helix, radii[2])), 12)]
        turns = [0; fill(inv(helix.dir * pitch(helix, radii[1])), 6);
                 fill(inv(helix.dir * pitch(helix, radii[2])), 12)]
        N = sum(G .* turns) / sum(G)
        factor = 1 + 2π^2 * N^2 * (11.0e-3^2 - 7.5e-3^2) / log(11.0e-3 / 7.5e-3)
        @test fill_layer.material.mu_r ≈ factor rtol = 1e-12
        @test insulation.material.mu_r ≈ polymer.mu_r * factor rtol = 1e-12
        @test insulation.material.eps_r == polymer.eps_r
        # A tube of one annulus has no turns, and its interval keeps its material.
        @test sheath.part.material === copper
        @test jacket.material === polymer

        reached = build(CableDesign, "reached",
            Terminal(:core, Group(Region(:core_fill, Circle(7.5e-3), air), strands)),
            Region(:insulation, Offset(3.0e-3), polymer))
        @test length(flatten(reached.tree).items) == 2

        nested = build(CableDesign, "nested",
            Stack(Terminal(:core, Region(:wire, Circle(2e-3), copper)),
                Region(:insulation, Offset(1e-3), polymer)),
            Region(:jacket, Offset(1e-3), polymer))
        @test length(flatten(nested.tree).items) == 3
    end

    @testset "an offset around a bar is the annulus of its thickness" begin
        design = build(CableDesign, "bar",
            Terminal(:core, Region(:bar, Rectangle(10e-3, 5e-3), copper)),
            Region(:insulation, Offset(1e-3), polymer))
        core, insulation = flatten(design.tree).items
        @test r_in(insulation.primitive) == r_ex(core.part.primitive)
        @test r_ex(insulation.primitive) ≈ 30e-3 / 2π + 1e-3
        @test insulation.material === polymer
    end

    @testset "a layer that is no annulus is refused" begin
        design = build(CableDesign, "tube_split_dielectric",
            Region(:bore, Circle(1e-3), polymer),
            Terminal(:core, Region(:tube, Annulus(1e-3, 2e-3), copper)),
            Repeat(CircularPattern(2; r = Offset(0)),
                Region(:insulation, Rectangle(2.5π * 1e-3, 1e-3), polymer)))
        @test_throws ArgumentError flatten(design.tree)
    end

    # Until the coaxial input reads the tree, this test compares each terminal of the
    # references with the tube of the blueprint. On Julia 1.12.7, x86_64 Linux, the largest difference was
    # 1.4e-13, in the permeability of a jacket of two metal layers.
    @testset "the terminals of the references give the tubes of the blueprint" begin
        cases = [REFERENCES;
                 [case => include(joinpath(
                      @__DIR__, "..", "references", case, "problem.jl"))
                  for case in ("cable_30kv_630mm2_trefoil",
                     "cable_30kv_na2xs2y_630mm2_trefoil", "cable_132kv_630mm2_flathor",
                     "cable_18kv_1000mm2_trefoil", "cable_380kv_2000mm2_flatver")]]
        for (case, problem) in cases, design in unique(problem.system.designs)

            # The cores of an enclosure come first, each the first part of its member.
            members = [first(member.second.items)
                       for item in design.tree.items if item isa Enclosure
                       for member in item.contents.members]
            own = [item for item in design.tree.items if item isa Terminal]
            tubes = flatten(design).conductors
            @test length(members) + length(own) == length(tubes)
            for (terminal, tube) in zip([members; own], tubes)
                region = flatten(terminal).first.part
                @test r_in(region.primitive)≈tube.r_in rtol=1e-12
                @test r_ex(region.primitive)≈tube.r_ex rtol=1e-12
                @test region.material.rho≈tube.material.rho rtol=1e-12
                @test region.material.mu_r≈tube.material.mu_r rtol=1e-12
                @test region.material.alpha≈tube.material.alpha rtol=1e-12
                @test region.material.T0 == tube.material.T0
            end
        end
    end

    # The layers of each design without sub-cables, against the rows of the blueprint that
    # have a width. The series of the layers of an interval gives its permeability.
    @testset "the stacks of the references give the layers of the blueprint" begin
        cases = [REFERENCES;
                 [case => include(joinpath(
                      @__DIR__, "..", "references", case, "problem.jl"))
                  for case in ("cable_30kv_630mm2_trefoil",
                     "cable_30kv_na2xs2y_630mm2_trefoil", "cable_132kv_630mm2_flathor",
                     "cable_18kv_1000mm2_trefoil", "cable_380kv_2000mm2_flatver")]]
        for (case, problem) in cases, design in unique(problem.system.designs)

            length(design.assemblies) == 1 || continue
            blueprint = flatten(design)
            parts = flatten(design.tree).items
            places = findall(part -> part isa Terminal, parts)
            @test length(places) == length(blueprint.conductors)
            for (k, place) in enumerate(places)
                stop = k < length(places) ? places[k + 1] - 1 : length(parts)
                layers = parts[(place + 1):stop]
                rows = [row
                        for row in blueprint.dielectrics[blueprint.dielectric_ranges[k]]
                        if row.r_ex - row.r_in > tolerance(row.r_ex)]
                @test length(layers) == length(rows)
                for (layer, row) in zip(layers, rows)
                    @test r_in(layer.primitive)≈row.r_in rtol=1e-12
                    @test r_ex(layer.primitive)≈row.r_ex rtol=1e-12
                    @test layer.material.rho == row.material.rho
                    @test layer.material.eps_r == row.material.eps_r
                    @test layer.material.tan_delta == row.material.tan_delta
                end
                isempty(layers) && continue
                r_e = r_ex(parts[place].part.primitive)
                series = sum(
                    layer.material.mu_r *
                    log(r_ex(layer.primitive) / r_in(layer.primitive))
                for layer in layers) /
                         log(r_ex(last(layers).primitive) / r_e)
                @test series≈blueprint.mu_r_ins[k] rtol=1e-12
            end
        end
    end
end
