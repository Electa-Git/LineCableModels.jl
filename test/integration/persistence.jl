@testitem "ImportExport / CablesLibrary / versioned JSON and trusted JLS" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures
] begin
    design=TestFixtures.mv_cable_design()
    library=CablesLibrary()
    add!(library, design)
    reference=compute(CableConstantsProblem(design), Formulation())

    mktempdir() do directory
        for extension in ("json", "jls")
            destination=joinpath(directory, "cables.$extension")
            @test save(library; file_name = destination) == destination
            @test isfile(destination)
            @test filesize(destination) > 0

            restored=CablesLibrary()
            @test load!(restored; file_name = destination) === restored
            @test collect(keys(restored)) == collect(keys(library))
            restored_design=restored[design.cable_id]
            @test restored_design !== design
            @test compute(CableConstantsProblem(restored_design), Formulation()) ==
                  reference
        end
    end
end

@testitem "ImportExport / CablesLibrary / failures are atomic" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures
] begin
    using JSON3
    using Serialization
    import LineCableModels.ImportExport as IE

    library=CablesLibrary()
    design=TestFixtures.mv_cable_design()
    add!(library, design)

    mktempdir() do directory
        function unchanged_after(path)
            before=library.data
            exception=try
                load!(library; file_name = path)
                nothing
            catch caught
                caught
            end
            @test exception !== nothing
            @test library.data === before
            @test library[design.cable_id] === design
            return exception
        end

        @test unchanged_after(joinpath(directory, "missing.json")) isa ArgumentError

        malformed=joinpath(directory, "malformed.json")
        write(malformed, "{not valid json")
        @test unchanged_after(malformed) isa Exception

        unsupported=joinpath(directory, "library.unsupported")
        write(unsupported, "fixture")
        @test unchanged_after(unsupported) isa ArgumentError

        invalid_jls=joinpath(directory, "invalid.jls")
        serialize(invalid_jls, 42)
        @test unchanged_after(invalid_jls) isa ArgumentError

        invalid_document=IE._json_document(library)
        invalid_document["cables"][design.cable_id]["type"]="Material"
        invalid_json=joinpath(directory, "invalid-entry.json")
        open(invalid_json, "w") do io
            JSON3.pretty(io, invalid_document)
        end
        @test unchanged_after(invalid_json) isa ArgumentError

        wrong_schema=IE._json_document(library)
        wrong_schema["schema"]="LineCableModels.MaterialsLibrary"
        schema_path=joinpath(directory, "wrong-schema.json")
        open(schema_path, "w") do io
            JSON3.pretty(io, wrong_schema)
        end
        @test unchanged_after(schema_path) isa ArgumentError

        legacy=joinpath(directory, "legacy.json")
        write(legacy, "{\"data\":{}}")
        before=library.data
        @test_logs (:warn, r"Unversioned LineCableModels JSON is retired") begin
            @test_throws ArgumentError load!(library; file_name = legacy)
        end
        @test library.data === before

        @test_throws ArgumentError save(
            library; file_name = joinpath(directory, "cables.archive")
        )
        @test_throws SystemError save(
            library; file_name = joinpath(directory, "absent", "cables.json")
        )
    end
end
