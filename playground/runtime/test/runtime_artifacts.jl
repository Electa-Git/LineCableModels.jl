import HTTP, JSON3

function artifact_fixture_job(;worker="worker-a")
    P=RT.Protocol
    fence=P.AssignmentFence(string(uuid4()),string(uuid4()),"alice","main",worker,string(uuid4()),
        string(uuid4()),"fixture","1.0.0",repeat("a",64),1)
    P.validate(P.AssignedJob("2.0",fence,P.new_job_request("fixture.echo",Dict("value"=>3);session_id=fence.run_id),
        P.PreparedExecution(string(uuid4()),1,repeat("b",64))))
end

@testset "private artifact configuration is passive and excludes credentials" begin
    mktempdir() do dir
        root=joinpath(dir,"private-artifacts")
        location=LocalRuntimeArtifacts(root)
        @test !ispath(root)
        @test !occursin(dir,repr(location))
        @test_throws ArgumentError LocalRuntimeArtifacts(tempdir())
        @test RT.configured_runtime_artifacts(Dict(),dir)===nothing
        @test RT.configured_runtime_artifacts(Dict("artifacts"=>Dict("backend"=>"filesystem","root"=>"data")),dir).root==joinpath(dir,"data")
        @test_throws ArgumentError RT.configured_runtime_artifacts(Dict("artifacts"=>Dict("backend"=>"filesystem","root"=>".")),dir)
        @test_throws ArgumentError RT.configured_runtime_artifacts(Dict("artifacts"=>Dict("backend"=>"filesystem","root"=>"data","url"=>"http://bad")),dir)
        credentials=joinpath(dir,"credentials.toml")
        write(credentials,"access_key_id = \"fixture-access\"\nsecret_access_key = \"fixture-secret\"\n");chmod(credentials,0o600)
        remote=S3RuntimeArtifacts("https://objects.invalid","lcm-private","runtime-v1",credentials)
        @test !occursin(dir,repr(remote)) && !occursin("objects.invalid",repr(remote))
        @test !occursin("fixture-secret",repr(RT.S3EndpointConfig(remote.endpoint,"fixture-access","fixture-secret")))
        @test_throws ArgumentError S3RuntimeArtifacts("http://objects.invalid","lcm-private","runtime-v1",credentials;allow_loopback_plaintext=true)
        @test_throws ArgumentError S3RuntimeArtifacts("http://127.0.0.1:9999","lcm-private","runtime-v1",credentials)
        @test_throws ArgumentError S3RuntimeArtifacts("https://user:secret@objects.invalid","lcm-private","runtime-v1",credentials)
        @test_throws ArgumentError S3RuntimeArtifacts("https://objects.invalid/path","lcm-private","runtime-v1",credentials)
        @test_throws ArgumentError S3RuntimeArtifacts("https://objects.invalid","lcm-private","../public",credentials)
        @test_throws ArgumentError S3RuntimeArtifacts("https://objects.invalid","lcm-private","",credentials)
        data=Dict("artifacts"=>Dict("backend"=>"s3","endpoint"=>remote.endpoint,"bucket"=>remote.bucket,
            "prefix"=>remote.prefix,"credentials_file"=>"credentials.toml"))
        @test RT.configured_runtime_artifacts(data,dir).credentials_file==credentials
        chmod(credentials,0o644)
        @test_throws ArgumentError RT.configured_runtime_artifacts(data,dir)
        @test !ispath(root)
    end
end

@testset "private filesystem artifacts preserve existing layout without public hash access" begin
    mktempdir() do dir
        location=LocalRuntimeArtifacts(joinpath(dir,"private"))
        job=artifact_fixture_job()
        bytes=collect(codeunits(JSON3.write(Dict("values"=>collect(1:30000)))))
        reference=RT.store_job_artifact!(location,job,bytes)
        @test reference.size>RT.ASSIGNED_INLINE_BYTES && reference.sha256==bytes2hex(RT.SHA.sha256(bytes))
        @test RT.read_job_artifact(location,job,reference)==bytes
        @test RT.store_job_artifact!(location,job,bytes)==reference
        @test RT.read_job_artifact(location,artifact_fixture_job(),reference)===nothing
        @test RT.read_job_artifact(location,artifact_fixture_job(;worker="worker-b"),reference)===nothing
        @test !ispath(joinpath(location.root,reference.sha256))
        scope=RT.artifact_scope(job.fence,job.request.job_id)
        parent=RT.private_artifact_directory(location,scope)
        @test sort(readdir(parent))==sort([reference.sha256,reference.sha256*".metadata.json"])
        @test stat(parent).mode & 0o777==0o700
        @test stat(joinpath(parent,reference.sha256)).mode & 0o777==0o600
        write(joinpath(parent,reference.sha256),"corrupt")
        @test_throws ArtifactUnavailable RT.read_job_artifact(location,job,reference)
        @test_throws ArtifactUnavailable RT.store_job_artifact!(location,job,bytes) # no silent replacement of committed data
        @test_throws ArtifactUnavailable RT.store_job_artifact!(location,job,zeros(UInt8,RT.ASSIGNED_ARTIFACT_BYTES+1))
        linked=LocalRuntimeArtifacts(joinpath(dir,"linked"))
        symlink(location.root,linked.root)
        @test_throws ArtifactUnavailable RT.store_job_artifact!(linked,job,bytes)
        @test read(joinpath(parent,reference.sha256),String)=="corrupt"
    end
end

@testset "bounded signed object transport uses job scopes and metadata-last commit" begin
    mktempdir() do dir
        objects=Dict{String,Vector{UInt8}}(); requests=Tuple{String,String}[]
        fail_metadata=Ref(false); redirect=Ref(false); oversized=Ref(false)
        server=HTTP.serve!("127.0.0.1",0;max_body_bytes=RT.ASSIGNED_ARTIFACT_BYTES) do request
            path=String(request.target)
            payload=request.body isa HTTP.EmptyBody ? UInt8[] : copy(request.body)
            push!(requests,(request.method,path))
            @test startswith(HTTP.header(request,"Authorization",""),"AWS4-HMAC-SHA256 Credential=fixture-access/")
            @test HTTP.header(request,"x-amz-content-sha256","")==bytes2hex(RT.SHA.sha256(payload))
            if redirect[]
                return HTTP.Response(307,["Location"=>"http://127.0.0.1:1/trap"])
            elseif oversized[]
                return HTTP.Response(200,fill(UInt8('x'),RT.ASSIGNED_ARTIFACT_BYTES+1))
            elseif request.method=="PUT"
                fail_metadata[] && occursin("/metadata/",path) && return HTTP.Response(503)
                objects[path]=copy(request.body);return HTTP.Response(200)
            elseif request.method=="DELETE"
                delete!(objects,path);return HTTP.Response(204)
            end
            HTTP.Response(haskey(objects,path) ? 200 : 404,get(objects,path,UInt8[]))
        end
        credentials=joinpath(dir,"credentials.toml")
        write(credentials,"access_key_id = \"fixture-access\"\nsecret_access_key = \"fixture-secret\"\n");chmod(credentials,0o600)
        location=S3RuntimeArtifacts("http://127.0.0.1:$(HTTP.port(server))","lcm-private","runtime-v1",credentials;allow_loopback_plaintext=true)
        job=artifact_fixture_job()
        bytes=collect(codeunits(JSON3.write(Dict("values"=>collect(1:30000)))))
        try
            reference=RT.store_job_artifact!(location,job,bytes)
            @test RT.read_job_artifact(location,job,reference)==bytes
            @test length(objects)==2
            @test all(path->startswith(path,"/lcm-private/runtime-v1/workers/worker-a/runs/$(job.fence.run_id)/"),keys(objects))
            @test endswith(last(requests)[2],"/sha256/"*reference.sha256)
            @test RT.read_job_artifact(location,artifact_fixture_job(),reference)===nothing
            @test length(objects)==2
            failed=artifact_fixture_job()
            fail_metadata[]=true
            @test_throws ArtifactUnavailable RT.store_job_artifact!(location,failed,bytes)
            @test length(objects)==2 # only the failed job's exact objects were removed
            @test count(p->p[1]=="DELETE",requests)==2
            fail_metadata[]=false
            redirect[]=true
            count_before=length(requests)
            @test_throws ArtifactUnavailable RT.read_job_artifact(location,job,reference)
            @test length(requests)==count_before+1 # no redirect or hidden transport retry
            redirect[]=false;oversized[]=true
            @test_throws ArtifactUnavailable RT.read_job_artifact(location,job,reference)
            oversized[]=false
            write(credentials,"secret_access_key = BROKEN_SECRET_CONTENT")
            error=try RT.read_job_artifact(location,job,reference);nothing catch e;e end
            @test error isa ArtifactUnavailable && !occursin("SECRET",sprint(showerror,error))
            for operation in (() -> RT.read_job_artifact(location,job,reference),
                    () -> RT.store_job_artifact!(location,job,bytes))
                failed=@async operation()
                @test_throws TaskFailedException fetch(failed)
                diagnostic=sprint(showerror,TaskFailedException(failed))
                @test !occursin("BROKEN_SECRET_CONTENT",diagnostic)
                @test !occursin(credentials,diagnostic) && !occursin(location.endpoint,diagnostic)
            end
        finally
            close(server)
        end
    end
end
