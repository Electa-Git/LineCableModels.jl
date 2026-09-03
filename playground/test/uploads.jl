using Bonito
using SHA

const LCMPlayground = LineCableModelsPlayground

function upload_staging_files(store)
    isdir(store.staging) || return String[]
    return filter(isfile, joinpath.(store.staging, readdir(store.staging)))
end

@testset "transactional upload store" begin
    mktempdir() do directory
        reject = Ref(false)
        validator = path -> begin
            reject[] && error("invalid specimen")
            startswith(read(path, String), "valid") || error("missing marker")
            return nothing
        end
        store = DirectoryUploadStore(directory)
        target = UploadTarget(
            store,
            "network-input";
            extensions=["CSV"],
            max_bytes=64,
            validate=validator
        )

        first_file = store_upload!(
            target,
            collect(codeunits("valid,first"));
            original_name="first.CSV",
            media_type="text/csv; charset=utf-8"
        )
        @test first_file.key == "network-input"
        @test first_file.original_name == "first.CSV"
        @test first_file.media_type == "text/csv"
        @test first_file.sha256 == bytes2hex(sha256("valid,first"))
        @test read(LCMPlayground.upload_path(target), String) == "valid,first"
        @test isempty(upload_staging_files(store))

        second_file = store_upload!(
            target,
            collect(codeunits("valid,second"));
            original_name="folder/second.csv",
            media_type="text/csv"
        )
        @test second_file.original_name == "second.csv"
        @test read(LCMPlayground.upload_path(target), String) == "valid,second"
        @test filter(!=(basename(store.staging)), readdir(directory)) == ["network-input"]

        reject[] = true
        @test_throws LCMPlayground.UploadRejected store_upload!(
            target,
            collect(codeunits("valid,rejected"));
            original_name="rejected.csv",
            media_type="text/csv"
        )
        @test read(LCMPlayground.upload_path(target), String) == "valid,second"
        @test isempty(upload_staging_files(store))

        @test_throws LCMPlayground.UploadRejected store_upload!(
            target,
            collect(codeunits("valid,bad-extension"));
            original_name="bad.json",
            media_type="application/json"
        )
        @test read(LCMPlayground.upload_path(target), String) == "valid,second"

        delete_upload!(target)
        @test !ispath(LCMPlayground.upload_path(target))
    end
end

@testset "upload target validation and stale staging cleanup" begin
    store = DirectoryUploadStore(mktempdir())
    @test_throws ArgumentError UploadTarget(store, "../escape")
    @test_throws ArgumentError UploadTarget(store, "nested/input")
    @test_throws ArgumentError UploadTarget(store, "input"; max_bytes=0)
    @test_throws ArgumentError UploadTarget(store, "input"; retention=:forever)
    @test_throws ArgumentError UploadTarget(store, "input"; media_types=["invalid"])

    mkpath(store.staging)
    abandoned = joinpath(store.staging, "abandoned")
    write(abandoned, "partial")
    @test sweep_upload_staging!(store; older_than_seconds=0) == 1
    @test !ispath(abandoned)
end

function upload_request(method, target, generation; name=nothing, body=UInt8[])
    headers = ["X-LCM-Upload-Generation" => string(generation)]
    if method == "POST"
        push!(headers, "X-LCM-Upload" => "1")
        push!(headers, "X-LCM-Original-Name" => something(name, "input.txt"))
        push!(headers, "Content-Type" => "text/plain")
    end
    return Bonito.HTTP.Request(method, target; headers, body)
end

@testset "upload gateway lifecycle" begin
    mktempdir() do directory
        registry = UploadRegistry()
        target = UploadTarget(
            DirectoryUploadStore(directory),
            "current-input";
            extensions=[".txt"],
            max_bytes=128
        )
        field = UploadField(target; registry)
        token = LCMPlayground.register_upload!(registry, target, field)
        gateway = LCMPlayground.UploadGateway(registry)
        endpoint = "/uploads/$token"

        request = upload_request(
            "POST",
            endpoint,
            1;
            name="input%20one.txt",
            body=collect(codeunits("first payload"))
        )
        response = Bonito.HTTPServer.apply_handler(gateway, (; request))
        @test response.status == 201
        @test field.status[] == :ready
        @test field.progress[] == 1.0
        @test field.value[].original_name == "input one.txt"
        @test read(LCMPlayground.upload_path(target), String) == "first payload"

        download = upload_request("GET", endpoint, 1)
        response = Bonito.HTTPServer.apply_handler(gateway, (; request=download))
        @test response.status == 200
        @test String(response.body) == "first payload"
        @test Bonito.HTTP.header(response, "Content-Type", "") == "text/plain"
        @test Bonito.HTTP.header(response, "Cache-Control", "") == "no-store"

        stale = Bonito.HTTPServer.apply_handler(gateway, (; request))
        @test stale.status == 409
        @test read(LCMPlayground.upload_path(target), String) == "first payload"

        rejected = upload_request(
            "POST",
            endpoint,
            2;
            name="wrong.csv",
            body=collect(codeunits("replacement"))
        )
        response = Bonito.HTTPServer.apply_handler(gateway, (; request=rejected))
        @test response.status == 415
        @test field.status[] == :failed
        @test read(LCMPlayground.upload_path(target), String) == "first payload"
        @test isempty(upload_staging_files(target.store))

        removal = upload_request("DELETE", endpoint, 3)
        response = Bonito.HTTPServer.apply_handler(gateway, (; request=removal))
        @test response.status == 200
        @test field.status[] == :empty
        @test isnothing(field.value[])
        @test !ispath(LCMPlayground.upload_path(target))

        response = Bonito.HTTPServer.apply_handler(gateway, (; request=download))
        @test response.status == 404

        LCMPlayground.unregister_upload!(registry, token)
        missing = upload_request("DELETE", endpoint, 4)
        response = Bonito.HTTPServer.apply_handler(gateway, (; request=missing))
        @test response.status == 404
        @test isempty(registry.entries)
    end
end

@testset "upload field render and session cleanup" begin
    mktempdir() do directory
        registry = UploadRegistry()
        target = UploadTarget(
            DirectoryUploadStore(directory),
            "temporary-input";
            extensions=[".txt"],
            retention=:session
        )
        initial = store_upload!(
            target,
            collect(codeunits("temporary"));
            original_name="temporary.txt",
            media_type="text/plain"
        )
        field = UploadField(target; initial, registry)
        session = Bonito.Session(
            Bonito.NoConnection();
            asset_server=Bonito.NoServer()
        )
        @test !isnothing(Bonito.jsrender(session, field))
        @test length(registry.entries) == 1
        @test isfile(LCMPlayground.upload_path(target))

        close(session)
        @test isempty(registry.entries)
        @test !ispath(LCMPlayground.upload_path(target))
        @test isnothing(field.value[])
    end
end
