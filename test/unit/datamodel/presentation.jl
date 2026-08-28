@testitem "DataModel / v1 presentation / model display and result table" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport,
    TestFixtures
] begin
    design=TestFixtures.mv_cable_design()
    system=TestFixtures.three_phase_system()

    for object in (design.root, first(design.geometry.regions), design, system)
        @test !isempty(sprint(show, MIME("text/plain"), object))
    end

    base=DataFrame(CableConstants(design))
    @test nrow(base) == 1
    @test names(base) == ["R", "L", "C"]

    design_display=sprint(show, MIME("text/plain"), design)
    @test contains(design_display, design.cable_id)
    @test contains(design_display, "regions=$(length(design.geometry.regions))")
    system_display=sprint(show, MIME("text/plain"), system)
    @test contains(system_display, system.system_id)
    @test contains(system_display, "cables=$(ncables(system))")

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
    library_display=sprint(show, MIME("text/plain"), library)
    @test contains(library_display, design.cable_id)
    @test contains(library_display, "1 design")
    @test_throws ArgumentError add!(library, design)
    @test delete!(library, design.cable_id) === library
    @test !haskey(library.catalogues, design.cable_id)
end
