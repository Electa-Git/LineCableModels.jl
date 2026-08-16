@testitem "ImportExport / CablesLibrary / JSON and Julia round trips" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures,
    TestAssertions
] begin
    design=TestFixtures.mv_cable_design()
    library=CablesLibrary()
    add!(library, design)
    reference_constants=compute!(CableConstantsProblem(design), Formulation())

    mktempdir() do directory
        for extension in ("json", "jls")
            destination=joinpath(directory, "cables.$extension")
            saved_path=@test_logs (:info, r"Cables library saved") save(
                library;
                file_name = destination
            )
            @test saved_path == destination
            @test isfile(saved_path)
            @test filesize(saved_path) > 0

            restored=CablesLibrary()
            @test_logs (:info, r"loaded|loading") match_mode = :any load!(
                restored;
                file_name = saved_path
            )
            @test collect(keys(restored)) == collect(keys(library))

            restored_design=restored[design.cable_id]
            @test restored_design !== design
            @test restored_design.cable_id == design.cable_id
            @test length(restored_design.components) == length(design.components)
            @test compute!(CableConstantsProblem(restored_design), Formulation()) ==
                  reference_constants
        end
    end
end

@testitem "ImportExport / CablesLibrary / malformed and missing inputs" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport
] begin
    mktempdir() do directory
        library=CablesLibrary()
        @test_throws ErrorException load!(
            library;
            file_name = joinpath(directory, "missing.json")
        )

        malformed=joinpath(directory, "malformed.json")
        write(malformed, "{not valid json")
        returned=redirect_stderr(devnull) do
            @test_logs (:error, r"Error loading CablesLibrary") match_mode = :any load!(
                library;
                file_name = malformed
            )
        end
        @test returned === library
        @test isempty(library.data)

        unsupported=joinpath(directory, "library.unsupported")
        write(unsupported, "fixture")
        returned=redirect_stderr(devnull) do
            @test_logs (
                :warn,
                r"Unrecognized file extension"
            ) (:error, r"Error loading CablesLibrary") match_mode = :any load!(
                library;
                file_name = unsupported
            )
        end
        @test returned === library
        @test isempty(library.data)
    end
end

@testitem "ImportExport / CablesLibrary / format dispatch and partial recovery" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures
] begin
    using JSON3
    using Logging
    using Serialization

    design=TestFixtures.mv_cable_design()
    source=CablesLibrary()
    add!(source, design)

    mktempdir() do directory
        fallback_request=joinpath(directory, "library.archive")
        fallback_path=@test_logs (
            :warn,
            r"Unrecognized file extension"
        ) (:info, r"saved") save(source; file_name = fallback_request)
        @test fallback_path == fallback_request * ".json"
        @test isfile(fallback_path)

        failed_save=redirect_stderr(devnull) do
            @test_logs (:error, r"Error saving CablesLibrary") save(
                source;
                file_name = joinpath(directory, "absent", "library.json")
            )
        end
        @test failed_save === nothing

        canonical=joinpath(directory, "canonical.json")
        @test_logs (:info, r"saved") save(source; file_name = canonical)
        parsed=JSON3.read(read(canonical, String), Dict{String, Any})
        direct=joinpath(directory, "direct.json")
        open(direct, "w") do io
            JSON3.pretty(io, parsed["data"])
        end
        restored=CablesLibrary()
        @test_logs match_mode = :any (
            :info,
            r"Assuming top-level JSON"
        ) (:info, r"Successfully loaded") load!(restored; file_name = direct)
        @test haskey(restored, design.cable_id)

        invalid_binary=joinpath(directory, "invalid.jls")
        serialize(invalid_binary, 42)
        @test_logs (:error, r"Invalid data format") load!(
            restored;
            file_name = invalid_binary
        )
        @test haskey(restored, design.cable_id)

        unreadable_binary=joinpath(directory, "unreadable.jls")
        write(unreadable_binary, UInt8[])
        returned=redirect_stderr(devnull) do
            @test_logs (:error, r"Error loading CablesLibrary") load!(
                restored;
                file_name = unreadable_binary
            )
        end
        @test returned === restored

        unrecognizable=joinpath(directory, "unrecognizable.json")
        write(unrecognizable, "{\"metadata\": 3}")
        @test_logs (:error, r"does not contain a recognizable") load!(
            restored;
            file_name = unrecognizable
        )
        @test isempty(restored.data)

        partial=joinpath(directory, "partial.json")
        open(partial, "w") do io
            JSON3.pretty(io, Dict("data"=>Dict("bad-entry"=>3)))
        end
        @test_logs match_mode = :any (
            :warn,
            r"Skipping entry 'bad-entry'"
        ) (:info, r"failed to load 1") load!(restored; file_name = partial)
        @test isempty(restored.data)

        broken_design=joinpath(directory, "broken-design.json")
        open(broken_design, "w") do io
            JSON3.pretty(io, Dict(
                "data"=>Dict(
                "broken"=>Dict(
                "nominal_data"=>3,
                "components"=>Any[]
            ),
            ),
            ))
        end
        redirect_stderr(devnull) do
            @test_logs match_mode = :any (
                :error,
                r"Failed to reconstruct cable design"
            ) (:info, r"failed to load 1") load!(restored; file_name = broken_design)
        end
        @test isempty(restored.data)
    end
end
