@testset "actual TLS object storage enforces private worker and reader identities" begin
    fence=P.AssignmentFence(string(uuid4()),string(uuid4()),"alice","main","worker-a",string(uuid4()),
        string(uuid4()),"fixture","1.0.0",repeat("a",64),1)
    job=P.AssignedJob("2.0",fence,P.new_job_request("fixture.echo",Dict("value"=>3);session_id=fence.run_id),
        P.PreparedExecution(string(uuid4()),1,repeat("b",64)))
    bytes=collect(codeunits(JSON3.write(Dict("values"=>collect(1:60000)))))
    writer=artifact_location("worker-a");reader=artifact_location("coordinator")
    reference=RT.store_job_artifact!(writer,job,bytes)
    @test RT.read_job_artifact(reader,job,reference)==bytes
    @test RT.store_job_artifact!(writer,job,bytes)==reference
    @test_throws ArtifactUnavailable RT.read_job_artifact(writer,job,reference) # write-only
    @test_throws ArtifactUnavailable RT.store_job_artifact!(reader,job,bytes) # read-only
    @test_throws ArtifactUnavailable RT.store_job_artifact!(artifact_location("worker-b"),job,bytes)
    @test RT.read_job_artifact(reader,job,reference)==bytes # failed foreign writes/cleanup cannot remove it
    other=P.AssignedJob("2.0",fence,P.new_job_request("fixture.echo",Dict("value"=>3);session_id=fence.run_id),job.execution)
    @test RT.read_job_artifact(reader,other,reference)===nothing # knowing a digest is not sufficient
    untrusted=S3RuntimeArtifacts(reader.endpoint,reader.bucket,reader.prefix,reader.credentials_file)
    @test_throws ArtifactUnavailable RT.read_job_artifact(untrusted,job,reference)
    scope=RT.artifact_scope(job.fence,job.request.job_id)
    public_url=reader.endpoint*"/"*reader.bucket*"/"*reader.prefix*"/"*scope*"/sha256/"*reference.sha256
    client=RT.HTTP.Client(transport=RT.HTTP.Transport(tls_config=RT.HTTP.TLS.Config(ca_file=reader.ca_file),proxy=nothing))
    try
        response=RT.HTTP.request("GET",public_url;client,retry=false,redirect=false,
            status_exception=false,request_timeout=5)
        @test response.status==403
    finally
        close(client)
    end
end
