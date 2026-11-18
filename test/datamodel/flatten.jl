using LineCableModels.DataModel: gmr

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
        core = flatten(first(design.tree.items))
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
                region = flatten(terminal).part
                @test r_in(region.primitive)≈tube.r_in rtol=1e-12
                @test r_ex(region.primitive)≈tube.r_ex rtol=1e-12
                @test region.material.rho≈tube.material.rho rtol=1e-12
                @test region.material.mu_r≈tube.material.mu_r rtol=1e-12
                @test region.material.alpha≈tube.material.alpha rtol=1e-12
                @test region.material.T0 == tube.material.T0
            end
        end
    end
end
