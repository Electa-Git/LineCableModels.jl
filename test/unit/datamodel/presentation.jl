@testitem "DataModel / presentation / libraries and tables" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport,
    TestFixtures
] begin
    design=TestFixtures.mv_cable_design()
    system=TestFixtures.three_phase_system()

    for object in (
        first(design.components).conductor_group.layers[1],
        first(design.components).conductor_group,
        first(design.components),
        design,
        system
    )
        text=sprint(show, MIME("text/plain"), object)
        @test !isempty(text)
        @test !occursin("temperature=", text)
    end

    @test nrow(DataFrame(design)) > 0
    @test nrow(DataFrame(design, :detailed)) > 0
    @test nrow(DataFrame(design, :baseparams)) == 3
    @test nrow(DataFrame(system)) == ncables(system)

    library=CablesLibrary()
    @test isempty(library)
    @test add!(library, design) === library
    @test library[design.cable_id] === design
    @test get(library, design.cable_id, nothing) === design
    sentinel=Ref(:missing)
    @test get(library, "missing", sentinel) === sentinel
    @test_throws ArgumentError add!(library, design)

    replacement=deepcopy(design)
    library[design.cable_id]=replacement
    @test library[design.cable_id] === replacement
    @test collect(keys(library)) == [design.cable_id]
    @test collect(values(library)) == [replacement]
    @test nrow(DataFrame(library)) == 1
    @test delete!(library, design.cable_id) === library
    @test isempty(library)
    @test_throws KeyError delete!(library, "missing")

    nominal=NominalData(
        designation_code = "sample",
        U0 = Float32(12),
        U = 20.0,
        resistance = 1
    )
    @test nominal isa NominalData{Float64}
    @test nominal.resistance === 1.0
    @test convert(NominalData{Float32}, nominal) isa NominalData{Float32}
end
