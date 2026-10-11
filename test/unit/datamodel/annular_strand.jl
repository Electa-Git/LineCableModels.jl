@testitem "DataModel / a complete bent strip retains annular geometry and lay" tags=[:unit, :importexport] begin
    const DM = LineCableModels.DataModel
    const EN = LineCableModels.Engine
    const IE = LineCableModels.ImportExport
    copper = Material(:conductor, 1.72e-8, 1.0, 1.0)
    strand = Rectangle(0.3e-3, 0.1e-3)
    center = Disk(0.2e-3)
    body = stranded(copper; shape=strand, center=center,
        boundary=Disk(0.6e-3), lay=LayRatio(12))
    design = build(CableDesign, "annular-last-strip", terminal(:core, body))
    metal = filter(region -> region.source.material.kind === :conductor,
        design.geometry.regions)
    @test length(metal) == 34
    final_strip = last(metal)
    @test final_strip.source.primitive == strand
    @test final_strip.primitive isa Annulus
    @test area(final_strip.primitive) ≈ area(strand)
    # A partial course member keeps its exact curved boundary. Directional
    # support must agree with that geometric boundary after a rigid placement too.
    bent = first(member.primitive for member in metal if member.primitive isa DM.BentStrip)
    for shape in (bent, DM.resolve(Pose2(0.002, -0.003, 0.47), bent))
        @test area(shape) ≈ area(strand)
        @test DM.perimeter(shape) ≈
            shape.span * (shape.ro + shape.ri) + 2(shape.ro - shape.ri)
        @test DM.thickness(shape) ≈ shape.ro - shape.ri
        @test DM.boundary(shape) === shape
        points = DM.tessellate(shape; points_per_arc=2049)
        for angle in (shape.at.φ, shape.at.φ + pi / 2, shape.at.φ + pi, -0.63)
            projection = maximum(point -> point[1] * cos(angle) + point[2] * sin(angle), points)
            @test DM.support(shape, angle) ≈ projection atol=2e-10 rtol=0
        end
    end
    @test perimeter(final_strip.primitive) ≈
        2pi * (final_strip.primitive.ri + final_strip.primitive.ro)
    @test only(final_strip.placement.patterns).pattern.course == 5
    @test only(final_strip.placement.patterns).member == 34
    for member in metal[2:end]
        @test only(member.paths).radius ≈ (r_in(member.primitive) + r_ex(member.primitive)) / 2
        @test overlength(only(member.paths).path, only(member.paths).radius) ≈ sqrt(1 + (pi / 12)^2)
    end
    @test IE.deserialize_value(IE.serialize_value(design)) == design
    @test length(design.geometry.regions) == length(metal)
    @test outer_radius(design) == final_strip.primitive.ro
    @test outer_radius(design) < 0.6e-3
    @test all(region -> only(region.placement.patterns).pattern.boundary ==
        design.geometry.outer, metal)
    blueprint = @inferred LineCableModels.DataModel.flatten(design, Float64)
    row = only(blueprint.conductors)
    @test row.r_ex == outer_radius(design)
    @test count(region -> region.terminal === :core, design.geometry.regions) ==
          length(metal)
    @test sum(region -> area(region.primitive), metal) ≈ area(center) + 33area(strand)
    expected_resistance = copper.rho / (area(center) + 33area(strand) / sqrt(1 + (pi / 12)^2))
    resistance = DM.tubular_resistance(row.r_in, row.r_ex, row.material.rho)
    gmr = DM.tubular_gmr(row.r_ex, row.r_in, row.material.mu_r)
    @test resistance ≈ expected_resistance
    @test isfinite(gmr) && 0 < gmr < row.r_ex

    shifted = build(CableDesign, "shifted-annular-last-strip",
        assembly(at(terminal(:core, body), Pose2(0.002, -0.003, 0.47))))
    shifted_blueprint = DM.flatten(shifted, Float64)
    shifted_row = only(shifted_blueprint.conductors)
    @test DM.tubular_resistance(shifted_row.r_in, shifted_row.r_ex,
        shifted_row.material.rho) ≈ resistance rtol = 1e-12
    @test DM.tubular_gmr(shifted_row.r_ex, shifted_row.r_in, shifted_row.material.mu_r) ≈
          gmr rtol = 1e-12
    # The turns affect the results only through the permeability of the dielectric interval.
    # With insulation, the lay raises it above the insulation's own, and the shifted copy
    # keeps the same value.
    dielectric = Material(:insulator, 1e14, 2.3, 1.0)
    insulated = build(CableDesign, "insulated-annular-last-strip",
        terminal(:core, body, insulation(dielectric; t=0.2e-3)))
    insulated_shifted = build(CableDesign, "shifted-insulated-annular-last-strip",
        assembly(at(terminal(:core, body, insulation(dielectric; t=0.2e-3)),
            Pose2(0.002, -0.003, 0.47))))
    insulated_mu = only(DM.flatten(insulated, Float64).mu_r_ins)
    @test insulated_mu > dielectric.mu_r
    @test only(DM.flatten(insulated_shifted, Float64).mu_r_ins) ≈ insulated_mu rtol = 1e-12

    designs = build(CableDesign, "annular-grid", terminal(:core,
        stranded(copper; shape=strand, center=center,
            boundary=Disk(Grid((0.59e-3, 0.6e-3))), lay=LayRatio(12))))
    @test designs isa Gridspace{CableDesign}
    @test length(designs) == 2
    rows = [only(LineCableModels.DataModel.flatten(selected, Float64).conductors)
        for selected in designs]
    @test [count(region -> region.terminal === :core, selected.geometry.regions)
           for selected in designs] == [33, 34]
    @test DM.tubular_resistance(last(rows).r_in, last(rows).r_ex, last(rows).material.rho) ≈
          resistance rtol = 1e-12
    @test DM.tubular_gmr(last(rows).r_ex, last(rows).r_in, last(rows).material.mu_r) ≈
          gmr rtol = 1e-12
end

@testitem "DataModel / annular mutual GMD uses area integration, not centroid distance" tags=[:unit, :datamodel] begin
    using QuadGK
    const DM = LineCableModels.DataModel
    outer = Annulus(1e-3, 2e-3)
    integral = first(quadgk(radius -> 2radius * log(radius), outer.ri, outer.ro;
        rtol=1e-12))
    expected = exp(integral / (outer.ro^2 - outer.ri^2))
    for inner in (Disk(0.4e-3), Disk(0.4e-3, Pose2(0.2e-3, -0.3e-3)),
            Annulus(0.1e-3, 0.8e-3),
            Rectangle(0.4e-3, 0.6e-3, Pose2(0.1e-3, 0.2e-3, 0.4)))
        @test DM.geometric_mean_distance(outer, inner) ≈ expected rtol=2e-12
        @test DM.geometric_mean_distance(inner, outer) ≈ expected rtol=2e-12
    end
    distant = Annulus(0.3e-3, 0.6e-3, Pose2(0.01, 0.0))
    @test DM.geometric_mean_distance(outer, distant) ≈ 0.01
    @test DM.geometric_mean_distance(distant, outer) ≈ 0.01
    @test_throws ArgumentError DM.geometric_mean_distance(Disk(0.4e-3), Disk(0.4e-3))
    for T in (Float32, Float64)
        annulus = Annulus(T(1e-3), T(1.0001e-3))
        center = Disk(T(0.5e-3))
        distance = @inferred DM.geometric_mean_distance(annulus, center)
        @test distance isa T
        @test annulus.ri <= distance <= annulus.ro
        reference = setprecision(BigFloat, 256) do
            inner, outer = BigFloat(annulus.ri), BigFloat(annulus.ro)
            exp((outer^2 * log(outer) - inner^2 * log(inner)) /
                (outer^2 - inner^2) - big"0.5")
        end
        @test distance ≈ T(reference) rtol=8eps(T)
    end
end

@testitem "DataModel / disk plus a uniform annular strand equals a solid conductor" tags=[:unit, :parametric] begin
    const EN = LineCableModels.Engine
    copper = Material(:conductor, 1.72e-8, 1.0, 1.0)
    center = Disk(0.2e-3)
    strand = Rectangle(1.4e-3, 0.1e-3)
    radius = sqrt(center.r^2 + area(strand) / pi)
    design = build(CableDesign, "one-annular-strand", terminal(:core,
        stranded(copper; shape=strand, center=center, boundary=Disk(radius))))
    row = only(LineCableModels.DataModel.flatten(design, Float64).conductors)
    # Independent uniform-current limiting case: there is no material or
    # current-density discontinuity between the center and its one full ring.
    @test count(region -> region.terminal === :core, design.geometry.regions) == 2
    @test sum(region -> area(region.primitive),
        filter(region -> region.terminal === :core, design.geometry.regions)) ≈
          pi * radius^2
    @test LineCableModels.DataModel.tubular_resistance(row.r_in, row.r_ex,
        row.material.rho) ≈ copper.rho / (pi * radius^2) rtol=1e-12
    @test LineCableModels.DataModel.tubular_gmr(row.r_ex, row.r_in, row.material.mu_r) ≈
          radius * exp(-1 / 4) rtol=2e-12
    solid_design = build(CableDesign, "equivalent-solid", terminal(:core,
        solid(copper, Disk(radius))))
    solid_row = only(LineCableModels.DataModel.flatten(solid_design, Float64).conductors)
    @test row.material.rho ≈ solid_row.material.rho
    @test row.material.mu_r ≈ solid_row.material.mu_r rtol=2e-11

    dielectric = Material(:insulator, Inf, 2.3, 1.0)
    insulated = map((design.origin, solid_design.origin)) do origin
        build(CableDesign, "uniform-current-limit", Stack(
            origin,
            insulation(dielectric; t=0.2e-3),
            terminal(:sheath, Region(:metal_screen, Shell(0.05e-3), copper)),
        ))
    end
    constants = map(source -> compute(CableConstantsProblem(source)), insulated)
    for field in (:R, :L, :C, :G)
        @test getproperty(first(constants), field) ≈ getproperty(last(constants), field) rtol=2e-10
    end
    parameters = map(insulated) do source
        system = build(LineCableSystem, source, Pose2(0.0, -1.0);
            connections=Dict(:core=>1, :sheath=>0))
        problem = LineParametersProblem(system; earth_props=homogeneous(rho=100.0),
            frequencies=[0.1, 50.0, 1e6])
        compute(problem)
    end
    @test first(parameters).Z ≈ last(parameters).Z rtol=2e-10
    @test first(parameters).Y ≈ last(parameters).Y rtol=2e-10
end

@testitem "DataModel / flatten / a helix on the axis is straight" tags=[:unit, :datamodel] begin
    library = MaterialsLibrary(add_defaults = true)
    copper, xlpe = Material(library, :copper), Material(library, :xlpe)
    insulated(core) = build(CableDesign, "axis-path",
        Stack(core, Region(:insulation, Shell(0.002), xlpe)))
    straight = LineCableModels.DataModel.flatten(
        insulated(Group(:core, Region(:metal, Disk(0.005), copper))))
    # On the axis, every lay law gives a straight path, with zero turns and overlength 1.
    for lay in (LayRatio(10.0), LayAngle(π / 4), Pitch(0.5))
        helical = LineCableModels.DataModel.flatten(insulated(
            Group(:core, Region(:metal, Disk(0.005), copper); path = Helix(lay))))
        @test helical.mu_r_ins == straight.mu_r_ins
        @test only(helical.conductors).material.rho == only(straight.conductors).material.rho
        @test overlength(Helix(lay), 0.0) == 1
        @test LineCableModels.DataModel.turns_per_length(
            ((path = Helix(lay), radius = 0.0),), Float64) == 0
    end
end
