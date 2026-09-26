@testitem "ParametricBuilder / Region finite inputs preserve targets and product/zip semantics" tags=[:unit] begin
    const EN = LineCableModels.Engine
    tags = (:wire_a, :wire_b)
    radii = (1e-3, 2e-3)
    resistivities = (1.72e-8, 2.82e-8)
    primitives = Disk.(radii)
    materials = Material.(Ref(:conductor), resistivities, 1.0, 1.0)
    tag_sources = (first(tags), Grid(tags), Gridspace{Symbol}(identity, (Grid(tags),)))
    primitive_sources = (first(primitives), Grid(primitives), Disk(Grid(radii)))
    material_sources = (first(materials), Grid(materials),
        Material(kind=:conductor, rho=Grid(resistivities)))

    for (tag_index, tag) in enumerate(tag_sources),
            (primitive_index, primitive) in enumerate(primitive_sources),
            (material_index, material) in enumerate(material_sources)
        regions = Region(tag, primitive, material)
        if (tag_index, primitive_index, material_index) == (1, 1, 1)
            @test regions isa Region
            @test regions.tag === first(tags)
            @test regions.primitive == first(primitives)
            @test regions.material == first(materials)
            continue
        end
        @test regions isa Gridspace{Region}
        @test Base.IteratorEltype(typeof(regions)) isa Base.HasEltype
        @test isconcretetype(eltype(regions))
        @test @inferred(first(regions)) isa eltype(regions)
        selected_tags = tag_index == 1 ? tags[1:1] : tags
        selected_radii = primitive_index == 1 ? radii[1:1] : radii
        selected_rho = material_index == 1 ? resistivities[1:1] : resistivities
        expected = vec(collect(Iterators.product(selected_tags, selected_radii, selected_rho)))
        actual = [(region.tag, region.primitive.r, region.material.rho) for region in regions]
        @test actual == expected

        zipped = Region(tag, primitive, material; combine=:zip)
        @test zipped isa Gridspace{Region}
        @test length(zipped) == 2
        @test isconcretetype(eltype(zipped))
        @test @inferred(first(zipped)) isa eltype(zipped)
        @test [(region.tag, region.primitive.r, region.material.rho) for region in zipped] ==
            [(tags[tag_index == 1 ? 1 : index], radii[primitive_index == 1 ? 1 : index],
                resistivities[material_index == 1 ? 1 : index]) for index in 1:2]
    end

    @test_throws DimensionMismatch Region(Grid(tags), Disk(Grid((1e-3, 2e-3, 3e-3))),
        first(materials); combine=:zip)
    deferred = Region(Grid((:valid, Symbol(""))), first(primitives), first(materials))
    @test deferred isa Gridspace{Region}
    @test first(deferred).tag === :valid
    @test_throws ArgumentError collect(deferred)
    empty_regions = Region(Grid(()), first(primitives), first(materials))
    @test isempty(empty_regions)
    @test eltype(empty_regions) === Any
    @test Base.IteratorEltype(typeof(empty_regions)) isa Base.EltypeUnknown

    selected = Region(Grid(tags), Disk(Grid(radii)), Grid(materials); combine=:zip)
    xlpe = Material(:insulator, Inf, 2.3, 1.0)
    designs = @cable "region-grid" begin
        @terminal :core begin
            selected
            insulation(xlpe; t=0.3e-3)
        end
        @terminal :sheath begin
            sheath(first(materials); t=0.1e-3)
        end
    end
    @test designs isa Gridspace{CableDesign}
    @test length(designs) == 2
    for (index, design) in enumerate(designs)
        @test design isa CableDesign
        metal = only(filter(region -> region.source.tag === tags[index],
            design.geometry.regions))
        @test metal.source.material.rho == resistivities[index]
        @test area(metal.primitive) ≈ pi * radii[index]^2
        blueprint = @inferred EN.flatten(LineCableModelsCoaxial(), design, Float64)
        @test length(blueprint.conductors) == 2
        @test first(blueprint.conductors).cross_section ≈ pi * radii[index]^2
        @test first(blueprint.conductors).resistance ≈
            resistivities[index] / (pi * radii[index]^2)
    end
end
