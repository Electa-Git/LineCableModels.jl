@testitem "DataModel / rectangular packing hands off its declared boundary" tags=[:unit, :importexport] begin
    const DM = LineCableModels.DataModel
    const EN = LineCableModels.Engine
    const IE = LineCableModels.ImportExport
    copper = Material(kind=:conductor, rho=1.72e-8)
    dielectric = Material(kind=:insulator, rho=Inf, eps_r=2.3)
    center, strand = Disk(0.2e-3), Rectangle(0.3e-3,0.1e-3)
    exact = sqrt(center.r^2 + 33area(strand)/pi)
    for limit in (0.59e-3, exact, exact/sqrt(1-1e-6), exact/sqrt(1-3e-6),
                  exact/sqrt(1-6e-6), 0.6e-3)
        body = stranded(copper; center=center, shape=strand, boundary=Disk(limit),
            lay=LayRatio(12), fill=dielectric)
        @test body isa Group
        @test body.boundary == Disk(limit)
        standalone = build(CableDesign,"rectangular-boundary",terminal(:core,body))
        regions = standalone.geometry.regions
        occupied = last(regions).primitive.ro
        @test outer_radius(standalone) == occupied == limit
        @test all(r -> r.source.material == copper && r.source.tag === :wire, regions)
        # The strips scale alike to fill the disk boundary, short of it by under one strand.
        @test area(first(regions).primitive) ≈ area(center)
        scale = area(last(regions).primitive) / area(strand)
        @test all(r -> area(r.primitive) ≈ scale * area(strand), regions[2:end])
        declared = sum(r -> area(r.source.primitive), regions)
        @test 0 <= area(Disk(limit)) - declared < area(strand)
        @test sum(area, regions) ≈ area(standalone.geometry.outer)
        @test all(r -> only(r.placement.patterns).pattern.boundary ==
            standalone.geometry.outer, regions)
        @test IE.deserialize_value(IE.serialize_value(standalone)) == standalone
        insulated = build(CableDesign,"rectangular-insulated",
            terminal(:core,body,insulation(dielectric;t=0.2e-3)))
        @test last(insulated.geometry.regions).primitive.ri == occupied
        @test outer_radius(insulated) ≈ occupied + 0.2e-3
        row = only(LineCableModels.DataModel.flatten(insulated,Float64).conductors)
        @test row.r_ex == occupied
        @test sum(area,filter(r -> r.terminal !== nothing,insulated.geometry.regions)) ≈
            sum(area,regions)
        @test all(r -> r.source.tag !== :stranded_fill, insulated.geometry.regions)
    end
    # Explicit enclosures remain physical even for sub-micrometer layers.
    body = stranded(copper;center=center,shape=strand,boundary=Disk(0.6e-3))
    for thickness in (1e-9, 4e-6)
        wrapped = build(CableDesign,"explicit-wrap",Enclosure(:paper,
            terminal(:core,body);primitive=Disk(0.6e-3+thickness),fill=dielectric))
        layer = last(wrapped.geometry.regions)
        @test layer.source.tag === :paper_fill
        @test layer.primitive isa Annulus
        @test layer.primitive.ri ≈ 0.6e-3
        @test DM.thickness(layer.primitive) ≈ thickness
    end
    rectangles = stranded(copper;center=center,
        shape=Rectangle(Grid((0.3e-3,0.31e-3)),0.1e-3),boundary=Disk(0.6e-3))
    @test rectangles isa Gridspace{AbstractCablePart}
    @test all(part -> part isa Group, rectangles)
    mixed = stranded(copper;center=center,
        shape=Grid((strand,Disk(0.1e-3))),boundary=Disk(0.6e-3))
    @test first(mixed) isa Group
    @test last(collect(mixed)) isa Enclosure
end

@testitem "DataModel / filled circular courses continue rectangular cores" tags=[:unit, :parametric] begin
    const DM = LineCableModels.DataModel
    const EN = LineCableModels.Engine
    copper = Material(kind=:conductor,rho=1.72e-8)
    matrix = Material(kind=:insulator,rho=1e10,eps_r=3.0)
    dielectric = Material(kind=:insulator,rho=Inf,eps_r=2.3)
    limit, wire_radius, count = 0.6e-3, 0.05e-3, 12
    body = stranded(copper;center=Disk(0.2e-3),shape=Rectangle(0.3e-3,0.1e-3),
        boundary=Disk(limit))
    core = terminal(:core,body)
    occupied = outer_radius(build(CableDesign,"core",core))
    for radius in (nothing, limit+wire_radius)
        ring = Group(:core,Region(:round_wire,Disk(wire_radius),copper);
            pattern=Ring(count;r=radius))
        filled = Enclosure(:ring_matrix,ring;
            primitive=Annulus(limit,limit+2wire_radius),fill=matrix)
        design = build(CableDesign,"rectangular-ring",core,filled,
            insulation(dielectric;t=0.1e-3))
        metal = filter(r -> r.source.tag === :round_wire,design.geometry.regions)
        @test length(metal) == count
        expected_radius = something(radius,occupied+wire_radius)
        @test all(r -> hypot(centroid(r.primitive)...) ≈ expected_radius,metal)
        fill_region = only(filter(r -> r.source.tag === :ring_matrix_fill,
            design.geometry.regions))
        @test fill_region.source.material == matrix
        @test fill_region.primitive isa DM.DifferenceShape
        @test fill_region.primitive.outer.ri == occupied
        @test fill_region.primitive.outer.ro == limit+2wire_radius
        @test length(fill_region.primitive.holes) == count
        @test area(fill_region.primitive) ≈
            pi*((limit+2wire_radius)^2-occupied^2)-count*pi*wire_radius^2
        @test sum(area,design.geometry.regions) ≈ area(design.geometry.outer)
        @test all(r -> r.source.tag !== :stranded_fill,design.geometry.regions)
        @test sum(area,filter(r -> r.terminal !== nothing,design.geometry.regions)) ≈
            pi*occupied^2+count*pi*wire_radius^2
        # The resistance is that of the declared rectangular strands.
        row = only(LineCableModels.DataModel.flatten(design,Float64).conductors)
        declared = area(Disk(0.2e-3))+33area(Rectangle(0.3e-3,0.1e-3))
        @test LineCableModels.DataModel.tubular_resistance(row.r_in,row.r_ex,
            row.material.rho) ≈ copper.rho/(declared+count*pi*wire_radius^2) rtol=1e-12
        # The transfer of geometry also works after a rigid placement of the complete stack.
        shifted = build(CableDesign,"rectangular-ring-placed",
            assembly(at(Stack(core,filled,insulation(dielectric;t=0.1e-3)),
                Pose2(0.01,-0.02,0.37))))
        @test sum(area,shifted.geometry.regions) ≈ area(shifted.geometry.outer)
    end
    ring = Group(:core,Region(:round_wire,Disk(wire_radius),copper);
        pattern=Ring(count;r=limit+wire_radius))
    # Only an otherwise feasible, concentric material fill can own the gap.
    @test_throws DomainError build(CableDesign,"overlapping-ring",core,
        Enclosure(:ring_matrix,ring;
            primitive=Annulus(occupied-wire_radius,limit+2wire_radius),fill=matrix))
    @test_throws DomainError build(CableDesign,"offset-ring",core,
        Enclosure(:ring_matrix,ring;at=Pose2(wire_radius,0.),
            primitive=Annulus(limit,limit+2wire_radius),fill=matrix))
    @test_throws DomainError build(CableDesign,"explicit-fill-region",core,
        Enclosure(:ring_matrix,ring;primitive=Annulus(limit,limit+2wire_radius),
            fill=Region(:specified_fill,Annulus(limit,limit+2wire_radius),matrix)))
end

@testitem "DataModel / rectangular occupied boundaries reconstruct under uncertainty" tags=[:extension, :measurements] begin
    using Measurements, Random
    const DM = LineCableModels.DataModel
    copper = Material(kind=:conductor,rho=1.72e-8)
    dielectric = Material(kind=:insulator,rho=Inf,eps_r=2.3)
    space = build(CableDesign,"rectangular-uncertain",terminal(:core,
        stranded(copper;center=Disk(0.2e-3),
            shape=Rectangle(Grid(0.3e-3,AbsoluteError(1e-6)),0.1e-3),
            boundary=Disk(Grid(0.6e-3,AbsoluteError(3e-6)))),
        insulation(dielectric;t=0.2e-3)))
    propagated = only(space)
    @test uncertainty(outer_radius(propagated)) > 0
    for design in (propagated, (rand(MersenneTwister(seed),space) for seed in 1:12)...)
        metal = filter(r -> r.source.material.kind === :conductor,design.geometry.regions)
        occupied = last(metal).primitive.ro
        @test last(design.geometry.regions).primitive.ri == occupied
        @test sum(area,metal) ≈ pi*occupied^2
        @test all(r -> r.source.tag !== :stranded_fill,design.geometry.regions)
        @test all(r -> only(r.placement.patterns).pattern.boundary.r == occupied,metal)
    end
    body = stranded(copper;center=Disk(0.2e-3),
        shape=Rectangle(Grid(0.3e-3,AbsoluteError(1e-6)),0.1e-3),
        boundary=Disk(0.6e-3))
    for radius in (nothing,0.65e-3)
        ring = Group(:core,Region(:round_wire,Disk(0.05e-3),copper);
            pattern=Ring(12;r=radius))
        filled = Enclosure(:ring_matrix,ring;
            primitive=Annulus(0.6e-3,0.7e-3),fill=dielectric)
        mixed = build(CableDesign,"uncertain-mixed-core",terminal(:core,body),filled,
            insulation(dielectric;t=0.1e-3))
        for design in (only(mixed),(rand(MersenneTwister(seed),mixed) for seed in 1:12)...)
            strips = filter(r -> r.source.primitive isa Rectangle,design.geometry.regions)
            occupied = last(strips).primitive.ro
            fill_region = only(filter(r -> r.source.tag === :ring_matrix_fill,
                design.geometry.regions))
            wires = filter(r -> r.source.tag === :round_wire,design.geometry.regions)
            @test fill_region.primitive.outer.ri == occupied
            @test all(r -> hypot(centroid(r.primitive)...) ≈
                something(radius,occupied+0.05e-3),wires)
            @test sum(area,design.geometry.regions) ≈ area(design.geometry.outer)
        end
    end
end
