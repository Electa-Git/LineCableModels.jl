using UUIDs
const Protocol = LineCableModelsPlaygroundProtocol
replace_record(value::T; changes...) where T =
    T((get(changes, field, getfield(value, field)) for field in fieldnames(T))...)

@testset "scientific commands are bounded and readiness is finite evidence" begin
    fence = AssignmentFence(string(uuid4()), string(uuid4()), "alice", "parameters",
        "worker-a", string(uuid4()), string(uuid4()), "line-parameters", "1.0.0", repeat("a",64), 1)
    command = ScientificCommand("2.0",string(uuid4()),fence,1,"prepare",Dict{String,Any}("length"=>10),nothing)
    status = replace_record(command; action="status",parameters=Dict{String,Any}())
    cancel = replace_record(status; action="cancel",target_id=string(uuid4()))
    report = ScientificReport("2.0",status.request_id,fence,1,true,"accepted","idle","ready",
        string(uuid4()),1,nothing,repeat("b",64),1000,1250,3,nothing,3000)
    cancel_job=replace_record(cancel;action="cancel_job")
    for record in (command,status,cancel,cancel_job,report)
        @test decode_runtime_message(typeof(record),encode_message(record)) == record
    end
    for changes in ((action="eval",), (revision=0,), (protocol_version="1.0",),
            (target_id=string(uuid4()),), (parameters=Dict{String,Any}("large"=>repeat("x",65536)),))
        @test_throws ArgumentError validate(replace_record(command;changes...))
    end
    @test_throws ArgumentError validate(replace_record(status;parameters=Dict{String,Any}("code"=>"1+1")))
    @test_throws ArgumentError validate(replace_record(cancel;target_id=nothing))
    @test_throws ArgumentError validate(replace_record(cancel;parameters=command.parameters))
    @test_throws ArgumentError validate(replace_record(cancel_job;target_id=nothing))
    @test_throws ArgumentError validate(replace_record(cancel_job;parameters=command.parameters))
    for changes in ((accepted=false,), (phase="executing",), (executor_id=nothing,), (executor_generation=0,),
            (preparation_key=nothing,), (valid_for_ms=0,), (valid_for_ms=5001,),
            (progress_milli=1001,), (output_lines=1_000_001,), (elapsed_ms=-1,),
            (reason="/private/path",), (failure="secret exception",), (preparation="installed",))
        @test_throws ArgumentError validate(replace_record(report;changes...))
    end
    cold = replace_record(report;preparation="cold",preparation_key=nothing,valid_for_ms=0)
    @test validate(cold) == cold
    @test_throws ArgumentError validate(replace_record(cold;valid_for_ms=1))
    @test_throws ArgumentError validate(replace_record(cold;preparation_key=report.preparation_key))
    for (record,field,value) in ((command,"revision",true),(report,"accepted",1),
            (report,"valid_for_ms",1.5),(report,"executor_generation",true),(command,"source","eval"))
        data=JSON3.read(encode_message(record),Dict{String,Any}); data[field]=value
        @test_throws ArgumentError decode_runtime_message(typeof(record),JSON3.write(data))
    end
    duplicate=replace(encode_message(status),"\"action\":\"status\""=>"\"action\":\"prepare\",\"action\":\"status\"")
    @test_throws ArgumentError decode_runtime_message(ScientificCommand,duplicate)
end

@testset "assigned runtime protocol is separate and strictly passive" begin
    fence = AssignmentFence(string(uuid4()), string(uuid4()), "alice", "parameters",
        "worker-a", string(uuid4()), string(uuid4()), "line-parameters", "1.0.0",
        repeat("a", 64), 1)
    profile = ProfileAdvertisement("line-parameters", "1.0.0", repeat("a", 64))
    report = WorkerAnnouncement("2.0", "worker-a", fence.worker_boot,
        fence.coordinator_id, string(uuid4()), 1, 2, [profile])
    probe = WorkerProbe("2.0", "worker-a", fence.coordinator_id, report.challenge)
    grant = LeaseControl("2.0", string(uuid4()), "grant", fence, 1, 30_000)
    ack = LeaseAcknowledgement("2.0", grant.request_id, fence, 1, true, "accepted")
    execution=PreparedExecution(string(uuid4()),1,repeat("b",64))
    job = AssignedJob("2.0", fence, new_job_request("system.echo",
        Dict("value"=>2); session_id=fence.run_id, timeout=Second(10)),execution)
    for value in (profile, report, probe, fence, grant, ack, job,execution)
        @test decode_runtime_message(typeof(value), encode_message(value)) == value
    end
    @test job.request.protocol_version == "1.0"
    @test_throws ArgumentError validate(replace_record(execution;executor_generation=0))
    @test_throws ArgumentError validate(replace_record(execution;executor_id="not-a-uuid"))
    @test_throws ArgumentError validate(replace_record(execution;preparation_key="not-a-key"))
    for mutate in (d->delete!(d,"execution"), d->(d["execution"]["executor_generation"]=true),
            d->(d["execution"]["privileged"]=true))
        data=JSON3.read(encode_message(job),Dict{String,Any}); mutate(data)
        @test_throws ArgumentError decode_runtime_message(AssignedJob,JSON3.write(data))
    end
    @test startswith(assigned_job_subject(fence), "lcm.jobs.v2.worker-a.")
    @test assigned_job_subject(fence) != assigned_job_subject(replace_record(fence; worker_id="worker-b"))
    @test assigned_job_subject(fence) != assigned_job_subject(replace_record(fence; generation=2))
    @test_throws ArgumentError validate(replace_record(job; protocol_version="1.0"))
    @test_throws ArgumentError validate(replace_record(job;
        request=replace_record(job.request; session_id="unrelated-run")))
    for value in ("worker.*", "worker.a", "WorkerA", "../worker")
        @test_throws ArgumentError validate(replace_record(report; worker_id=value))
    end
    for value in ("1", string(uuid4()) * "x", uppercase(fence.worker_boot))
        @test_throws ArgumentError validate(replace_record(fence; worker_boot=value))
    end
    @test_throws ArgumentError validate(replace_record(fence; generation=0))
    @test_throws ArgumentError validate(replace_record(fence; owner="alice\nadmin"))
    @test_throws ArgumentError validate(replace_record(profile; fingerprint="not-a-digest"))
    @test_throws ArgumentError validate(replace_record(profile; version="1"))
    @test_throws ArgumentError validate(replace_record(report; profiles=[profile, profile]))
    @test_throws ArgumentError validate(replace_record(report; sequence=0))
    @test_throws ArgumentError validate(replace_record(report; capacity=257))
    @test isempty(validate(replace_record(report; profiles=ProfileAdvertisement[])).profiles)
    @test_throws ArgumentError validate(replace_record(grant; action="eval"))
    @test_throws ArgumentError validate(replace_record(grant; duration_ms=60_001))
    @test_throws ArgumentError validate(replace_record(grant; duration_ms=0))
    @test_throws ArgumentError validate(replace_record(grant; action="release"))
    @test validate(replace_record(grant; action="release", duration_ms=0)).action == "release"
    @test_throws ArgumentError validate(replace_record(ack; reason="/private/path"))
    @test_throws ArgumentError decode_runtime_message(LeaseControl, fill(UInt8(' '), 256 * 1024 + 1))
    @test_throws ArgumentError decode_runtime_message(LeaseControl, "{")
    for change in (
            data -> (data["authority"]="administrator"),
            data -> delete!(data, "request_id"),
            data -> (data["revision"]=true),
            data -> (data["duration_ms"]=1.5),
            data -> (data["fence"]["generation"]=true),
            data -> (data["fence"]["image"]="arbitrary-image"),
        )
        data = JSON3.read(encode_message(grant), Dict{String,Any})
        change(data)
        @test_throws ArgumentError decode_runtime_message(LeaseControl, JSON3.write(data))
    end
    duplicate = replace(encode_message(grant),
        "\"action\":\"grant\"" => "\"action\":\"release\",\"action\":\"grant\"")
    @test_throws ArgumentError decode_runtime_message(LeaseControl, duplicate)
    false_acceptance = replace(encode_message(ack), "\"accepted\":true" => "\"accepted\":1")
    @test_throws ArgumentError decode_runtime_message(LeaseAcknowledgement, false_acceptance)
end

@testset "assigned results preserve provenance and reject foreign worker environments" begin
    fence = AssignmentFence(string(uuid4()), string(uuid4()), "alice", "parameters",
        "worker-a", string(uuid4()), string(uuid4()), "line-parameters", "1.0.0",
        repeat("a", 64), 1)
    request = new_job_request("system.echo", Dict("message"=>"hello"); session_id=fence.run_id)
    result = JobResult("1.0", request.job_id, request.operation, "1.0", request.input_hash,
        "fixture", fence.fingerprint, fence.worker_id, "miss", utc_timestamp(), utc_timestamp(),
        Dict{String,Any}("value"=>2.0), nothing, nothing, String[])
    artifact = ArtifactReference("sha256:" * repeat("b", 64), "application/json", 128,
        repeat("b", 64), "filesystem", "/artifacts/sha256/" * repeat("b", 64))
    failure = FailureInfo("fixture", "Expected fixture failure", string(uuid4()), false)
    variants = (result, replace_record(result; inline_result=nothing, artifact=artifact),
        replace_record(result; inline_result=nothing, failure=failure))
    execution=PreparedExecution(string(uuid4()),1,repeat("b",64))
    for variant in variants
        outcome = AssignedResult("2.0", fence, variant,execution)
        @test decode_runtime_message(AssignedResult, encode_message(outcome)) == outcome
    end
    outcome = AssignedResult("2.0", fence, result,execution)
    @test startswith(assigned_result_subject(fence, request.job_id), "lcm.results.v2.worker-a.")
    @test endswith(assigned_result_subject(fence, request.job_id), "." * request.job_id)
    @test_throws ArgumentError assigned_result_subject(fence, "job.*")
    @test_throws ArgumentError validate(replace_record(outcome; protocol_version="1.0"))
    @test_throws ArgumentError validate(replace_record(outcome;
        result=replace_record(result; worker_id="worker-b")))
    @test_throws ArgumentError validate(replace_record(outcome;
        result=replace_record(result; environment_fingerprint=repeat("b", 64))))
    malformed = JSON3.read(encode_message(outcome), Dict{String,Any})
    malformed["result"]["owner"] = "bob"
    @test_throws ArgumentError decode_runtime_message(AssignedResult, JSON3.write(malformed))
    bad_size = JSON3.read(encode_message(AssignedResult("2.0", fence, variants[2],execution)), Dict{String,Any})
    bad_size["result"]["artifact"]["size"] = true
    @test_throws ArgumentError decode_runtime_message(AssignedResult, JSON3.write(bad_size))
    bad_retry = JSON3.read(encode_message(AssignedResult("2.0", fence, variants[3],execution)), Dict{String,Any})
    bad_retry["result"]["failure"]["retryable"] = 1
    @test_throws ArgumentError decode_runtime_message(AssignedResult, JSON3.write(bad_retry))
end
