@testset "same-origin artifact gateway" begin
    mktempdir() do directory
        digest = repeat("a", 64)
        bytes = collect(codeunits("artifact payload"))
        write(joinpath(directory, digest), bytes)
        write(
            joinpath(directory, "$digest.metadata.json"),
            "{\"media_type\":\"application/json\",\"size\":$(length(bytes)),\"sha256\":\"$digest\"}"
        )
        gateway = ArtifactGateway(directory)
        context(method, target) = (
            request=LineCableModelsPlayground.Bonito.HTTP.Request(method, target),
        )
        response = LineCableModelsPlayground.Bonito.HTTPServer.apply_handler(
            gateway,
            context("GET", "/artifacts/sha256/$digest")
        )
        @test response.status == 200
        @test response.body == bytes
        @test LineCableModelsPlayground.Bonito.HTTP.header(
            response,
            "Content-Type"
        ) == "application/json"

        range_context = (
            request=LineCableModelsPlayground.Bonito.HTTP.Request(
                "GET",
                "/artifacts/sha256/$digest",
                ["Range" => "bytes=1-7"]
            ),
        )
        partial = LineCableModelsPlayground.Bonito.HTTPServer.apply_handler(
            gateway,
            range_context
        )
        @test partial.status == 206
        @test collect(partial.body) == bytes[2:8]
        @test LineCableModelsPlayground.Bonito.HTTP.header(
            partial,
            "Content-Range"
        ) == "bytes 1-7/$(length(bytes))"

        invalid_range = (
            request=LineCableModelsPlayground.Bonito.HTTP.Request(
                "GET",
                "/artifacts/sha256/$digest",
                ["Range" => "bytes=999-1000"]
            ),
        )
        @test LineCableModelsPlayground.Bonito.HTTPServer.apply_handler(
            gateway,
            invalid_range
        ).status == 416

        head = LineCableModelsPlayground.Bonito.HTTPServer.apply_handler(
            gateway,
            context("HEAD", "/artifacts/sha256/$digest")
        )
        @test head.status == 200
        @test isempty(head.body)

        traversal = LineCableModelsPlayground.Bonito.HTTPServer.apply_handler(
            gateway,
            context("GET", "/artifacts/sha256/../../Project.toml")
        )
        @test traversal.status == 404

        rejected = LineCableModelsPlayground.Bonito.HTTPServer.apply_handler(
            gateway,
            context("POST", "/artifacts/sha256/$digest")
        )
        @test rejected.status == 405

        @test isnothing(LineCableModelsPlayground.validate_artifact_metadata(
            (
                sha256=digest,
                size=length(bytes),
                media_type="text/plain\r\nX-Injected: true",
            ),
            digest
        ))
        @test isnothing(LineCableModelsPlayground.validate_artifact_metadata(
            (sha256=digest, size=Int128(1) << 41, media_type="application/json"),
            digest
        ))
    end

    @test_throws ArgumentError LineCableModelsPlayground.S3ArtifactGateway(
        "http://127.0.0.1:9000",
        "artifacts",
        "reader",
        "secret"
    )
    @test_throws ArgumentError LineCableModelsPlayground.S3ArtifactGateway(
        "https://objects.example.invalid",
        "artifacts",
        "reader",
        "secret";
        prefix="../escape"
    )

    attempts = Ref(0)
    @test LineCableModelsPlayground.with_s3_transport_retry(attempts=3) do
        attempts[] += 1
        attempts[] < 3 && error("transient transport failure")
        :recovered
    end == :recovered
    @test attempts[] == 3
    @test_throws ArgumentError LineCableModelsPlayground.with_s3_transport_retry(
        () -> nothing;
        attempts=0
    )
end
