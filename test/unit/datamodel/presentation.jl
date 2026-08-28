@testitem "DataModel / v1 presentation / physical tables and library" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport,
    TestFixtures
] begin
    design=TestFixtures.mv_cable_design()
    system=TestFixtures.three_phase_system()

    for object in (design.root, first(design.geometry.regions), design, system)
        @test !isempty(sprint(show, MIME("text/plain"), object))
    end

    regions=DataFrame(design)
    @test nrow(regions) == length(design.geometry.regions)
    @test names(regions) == [
        "tag", "terminal", "primitive", "material_kind", "area",
        "centroid_x", "centroid_y", "overlength"
    ]
    base=DataFrame(CableConstants(design))
    @test nrow(base) == 1
    @test names(base) == ["R", "L", "C"]

    system_table=DataFrame(system)
    @test nrow(system_table) == ncables(system)
    @test system_table.cable_id == getproperty.(system.designs, :cable_id)

    library=CablesLibrary()
    catalogue_record=(
        designation_code = "sample",
        U0 = Float32(12),
        U = 20.0,
        resistance = 1
    )
    @test add!(library, design; catalogue = catalogue_record) === library
    @test library[design.cable_id] === design
    @test catalogue(library, design.cable_id) == catalogue_record
    @test nrow(DataFrame(library)) == 1
    @test names(DataFrame(library)) == ["cable_id", "catalogue", "terminals"]
    @test_throws ArgumentError add!(library, design)
    @test delete!(library, design.cable_id) === library
    @test !haskey(library.catalogues, design.cable_id)
end
