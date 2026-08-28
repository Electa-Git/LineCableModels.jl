@testitem "DataModel / v1 presentation / eager tables and library" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport,
    TestFixtures
] begin
    design=TestFixtures.mv_cable_design()
    system=TestFixtures.three_phase_system()

    for object in (design.root, first(design.geometry.regions), design, system)
        @test !isempty(sprint(show, MIME("text/plain"), object))
    end

    terminals=DataFrame(design)
    @test nrow(terminals) == length(design.terminal_order)
    @test terminals.terminal == design.terminal_order
    regions=DataFrame(design, :regions)
    @test nrow(regions) == length(design.geometry.regions)
    @test names(regions) == [
        "tag", "terminal", "primitive", "material_kind", "area",
        "centroid_x", "centroid_y", "overlength"
    ]
    base=DataFrame(design, :baseparams)
    @test nrow(base) == 1
    @test names(base) == ["R", "L", "C"]
    @test_throws ArgumentError DataFrame(design, :unknown)

    system_table=DataFrame(system)
    @test nrow(system_table) == ncables(system)
    @test system_table.cable_id == getproperty.(system.designs, :cable_id)

    library=CablesLibrary()
    @test add!(library, design) === library
    @test library[design.cable_id] === design
    @test nrow(DataFrame(library)) == 1
    @test names(DataFrame(library)) == ["cable_id", "nominal_data", "terminals"]
    @test_throws ArgumentError add!(library, design)
    @test delete!(library, design.cable_id) === library

    nominal=NominalData(
        designation_code = "sample",
        U0 = Float32(12),
        U = 20.0,
        resistance = 1
    )
    @test nominal isa NominalData{Float64}
    @test convert(NominalData{Float32}, nominal) isa NominalData{Float32}
end
