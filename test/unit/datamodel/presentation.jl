@testitem "DataModel / v1 presentation / model display and result table" tags=[:unit] setup=[
    UseDataModelSupport,
    TestFixtures
] begin
    design=TestFixtures.coaxial_design()
    system=TestFixtures.three_phase_system()
    @test hasfield(typeof(design), :origin)
    @test hasproperty(design, :origin)
    @test !hasfield(typeof(design), :root)
    @test !hasproperty(design, :root)
    @test_throws FieldError getproperty(design, :root)
    @test build(CableDesign, design.cable_id, design.origin;
        nominal_data=design.nominal_data) == design

    for object in (design.origin, first(design.geometry.regions), design, system)
        @test !isempty(sprint(show, MIME("text/plain"), object))
    end

    constants=CableConstants(design)
    base=DataFrame(observables(constants, (R, L, C, G)))
    @test nrow(base) == 1
    @test names(base) == ["core", "R", "L", "C", "G"]

    design_display=sprint(show, MIME("text/plain"), design)
    @test contains(design_display, design.cable_id)
    @test contains(design_display, "regions    $(length(design.geometry.regions))")
    @test contains(design_display, "origin     ")
    @test !contains(design_display, "root       ")
    system_display=sprint(show, MIME("text/plain"), system)
    @test contains(system_display, system.system_id)
    @test contains(system_display, "cables")

    library=CablesLibrary()
    @test library isa AbstractDict{String, CableDesign}
    @test keytype(typeof(library)) === String
    @test valtype(typeof(library)) === CableDesign
    @test eltype(typeof(library)) === Pair{String, CableDesign}
    @test pairs(library) === library
    catalogue_record=DatasheetInfo(
        designation_code = "sample",
        U0 = Float32(12),
        U = 20.0,
        resistance = 1
    )
    @test add!(library, design; catalogue = catalogue_record) === library
    @test library[design.cable_id] === design
    @test collect(library) isa Vector{Pair{String, CableDesign}}
    @test catalogue(library, design.cable_id) == catalogue_record
    @test get(library, "missing", nothing) === nothing
    @test get(() -> design, library, "missing") === design
    @test get(() -> design, library, design.cable_id) === design
    @test_throws MethodError get(library, "missing")
    @test !haskey(library, Symbol(design.cable_id))
    @test get(library, Symbol(design.cable_id), nothing) === nothing
    @test_throws KeyError library[Symbol(design.cable_id)]
    @test propertynames(catalogue_record) == (:designation_code, :U0, :U, :resistance)
    catalogue_display=sprint(show, MIME("text/plain"), catalogue_record)
    @test contains(catalogue_display, "U0")
    @test contains(catalogue_display, "12 kV")
    @test contains(catalogue_display, "1 Ω/km")
    library_display=sprint(show, MIME("text/plain"), library)
    @test contains(library_display, design.cable_id)
    @test contains(library_display, "1 design")

    copied=copy(library)
    @test copied isa CablesLibrary
    @test copied !== library
    @test copied.data !== library.data
    @test copied.catalogues !== library.catalogues
    @test copied[design.cable_id] === design
    @test catalogue(copied, design.cable_id) == catalogue_record
    @test delete!(copied, design.cable_id) === copied
    @test haskey(library, design.cable_id)

    blank=empty(library)
    @test blank isa CablesLibrary
    @test isempty(blank)
    @test isempty(blank.catalogues)

    @test_throws ArgumentError add!(library, design)
    @test_throws ArgumentError setindex!(library, design, "wrong-id")
    @test setindex!(library, design, design.cable_id) === library
    @test (library[design.cable_id] = design) === design
    @test catalogue(library, design.cable_id) == DatasheetInfo(design.nominal_data)
    @test delete!(library, "missing") === library
    @test delete!(library, :missing) === library
    @test delete!(library, design.cable_id) === library
    @test !haskey(library.catalogues, design.cable_id)
    @test add!(library, design; catalogue = catalogue_record) === library
    @test empty!(library) === library
    @test isempty(library)
    @test isempty(library.catalogues)
end
