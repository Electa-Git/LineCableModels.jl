import HTTP
include("scientific_driver_fixture.jl")

@testset "TLS scientific subjects cannot self-authorize or cross worker identities" begin
    coordinator=BrokerControl(endpoint("coordinator"),CoordinatorIdentity();worker_ids=("worker-a","worker-b"),traffic=:science)
    a=b=nothing
    try
        a=BrokerControl(endpoint("worker-a"),WorkerIdentity("worker-a");traffic=:science)
        b=BrokerControl(endpoint("worker-b"),WorkerIdentity("worker-b");traffic=:science)
        fence=P.AssignmentFence(string(uuid4()),string(uuid4()),"alice","main","worker-a",string(uuid4()),
            string(uuid4()),"fixture","1.0.0",repeat("a",64),1)
        command=P.ScientificCommand("2.0",string(uuid4()),fence,1,"status",Dict{String,Any}(),nothing)
        send_control!(coordinator,command)
        @test only(records(a)).record==command
        @test isempty(poll_control!(b))
        report=RT.unavailable_scientific_report(command,"cold")
        send_control!(a,report)
        @test only(records(coordinator)).record==report
        @test_throws AccessDenied send_control!(b,report)
        NATS.publish(a.connection,"lcm.science.v2.worker-a.command",P.encode_message(command))
        NATS.ping(a.connection;measure=false)
        @test isempty(records(a;timeout=0.2))
        NATS.publish(a.connection,"lcm.science.v2.worker-b.report",P.encode_message(report))
        NATS.ping(a.connection;measure=false)
        @test isempty(records(coordinator;timeout=0.2))
        NATS.publish(a.connection,"lcm.science.v2.worker-a.report",P.encode_message(
            P.ScientificReport(report.protocol_version,report.request_id,
                P.AssignmentFence(fence.lease_id,fence.run_id,fence.owner,fence.role,"worker-b",fence.worker_boot,
                    fence.coordinator_id,fence.profile_id,fence.profile_version,fence.fingerprint,fence.generation),
                report.revision,report.accepted,report.reason,report.phase,report.preparation,report.executor_id,
                report.executor_generation,report.current_request_id,report.preparation_key,report.progress_milli,
                report.elapsed_ms,report.output_lines,report.failure,report.valid_for_ms)))
        @test isempty(records(coordinator;timeout=0.2))
        @test coordinator.rejected==1
    finally
        a===nothing || close(a); b===nothing || close(b); close(coordinator)
    end
end

@testset "protected API prepares and cancels a remote finite child through production schedulers" begin
    # The child is explicitly a test fixture, not a claim of OS-quota admission.
    # Control/API/transport/resource ownership below are the production code.
    mktempdir() do dir
        store=RuntimeStore(joinpath(dir,"runtime.sqlite"))
        profiles,applications=ProfileRegistry(),ApplicationRegistry()
        project=normpath(joinpath(@__DIR__,"..","..","worker","core"))
        fingerprint=native_environment_fingerprint(project)
        register!(profiles,ProfileDefinition("fixture",project,fingerprint.digest;
            operations=("fixture.echo","fixture.delay","fixture.evict","fixture.fail","fixture.large"),
            budget=ResourceBudget(prepare_seconds=30,job_seconds=30)))
        driver=ScientificDriverFixture(profiles,Dict("fixture"=>fixture_command(project,
            joinpath(@__DIR__,"scientific_child.jl"),"300")),
            Dict{String,Tuple{P.AssignmentFence,RT.ExecutionCore.ExecutorSupervisor}}(),0,0,false,true,false)
        resources=ScientificResources(driver)
        definition=ApplicationDefinition("science","Science",:workbench,"/science";
            requirements=(RuntimeRequirement("main",("fixture",)),))
        register!(applications,definition)
        alice=Principal("alice"); operator=Principal("operator";administrator=true)
        trust=WorkerTrust("worker-a","credential-worker-a",("fixture",))
        enroll_worker!(store,operator,trust)
        set_registration_state!(store,operator,"worker-a",:approved;expected_revision=1)
        control=ControlService(ControlConfig(endpoint("coordinator"),profiles,[trust];artifacts=artifact_location("coordinator")),store,applications)
        agent=AgentService(AgentConfig("worker-a",endpoint("worker-a"),profiles,joinpath(dir,"agent");artifacts=artifact_location("worker-a")),resources)
        supervisor=UIHostSupervisor(store,applications,joinpath(dir,"hosts"))
        key=repeat("fixture-proxy-key-",3)
        policy=ProxyIdentity("https://lcm.test",["127.0.0.1"],key)
        server=start_gateway(supervisor,policy;control)
        jobs=BrokerJobs(endpoint("coordinator"),CoordinatorIdentity())
        base="http://127.0.0.1:$(HTTP.port(server))"
        headers(owner)=["X-LCM-Proxy-Key"=>key,"X-LCM-Principal"=>owner]
        mutation(owner)=[headers(owner);"Origin"=>"https://lcm.test";"X-LCM-Request"=>"1"]
        request(method,path,hs=headers("alice"),body=nothing)=HTTP.request(method,base*path,hs,body;
            proxy=nothing,retry=false,status_exception=false,request_timeout=10)
        snapshot(path)=JSON3.read(request("GET",path).body)
        try
            start_agent!(agent); start_control!(control)
            @test request("GET","/health",Pair{String,String}[]).status==200
            @test timedwait(()->haskey(control.inventory.presence,"worker-a") &&
                control.inventory.presence["worker-a"].report.sequence>=2 &&
                control.science.state==:online && agent.science.state==:online,25)==:ok
            @test control.science.connection.connection!==control.link.control.connection
            @test agent.science.connection.connection!==agent.link.control.connection
            run=reserve_run!(store,alice,definition)
            lease=reserve_assignment!(control.coordinator.assignments,alice,run.id,"main","fixture";
                placement=PinnedPlacement("worker-a"))
            id=UUID(lease.fence.lease_id)
            grant_assignment!(control.coordinator,alice,id)
            @test timedwait(()->assignment_usable(control.coordinator,alice,id),5)==:ok
            path="/runtime/api/assignments/$id/science"
            @test request("GET",path,headers("bob")).status==404
            @test timedwait(()->snapshot(path).preparation=="cold",8;pollint=0.2)==:ok
            @test isempty(driver.processes) && isempty(supervisor.handles)
            jobs_path="/runtime/api/assignments/$id/jobs"
            cold_job=JSON3.write((operation="fixture.echo",parameters=Dict("value"=>3),request_id=string(uuid4())))
            @test request("POST",jobs_path,mutation("alice"),cold_job).status==409
            @test request("POST",jobs_path,mutation("bob"),"{").status==404
            @test isempty(list_jobs(store,alice,run.id))
            prepare_id=string(uuid4())
            body=JSON3.write((action="prepare",parameters=Dict("seconds"=>0.1),request_id=prepare_id))
            @test request("POST",path,headers("alice"),body).status==403
            accepted=request("POST",path,mutation("alice"),body)
            @test accepted.status==202 && JSON3.read(accepted.body).preparation!="ready"
            @test timedwait(()->snapshot(path).preparation=="ready",35;pollint=0.2)==:ok
            ready=snapshot(path)
            @test ready.preparation=="ready" && 0<ready.valid_for_ms<=5000
            @test ready.executor_generation==1 && ready.executor_id!==nothing
            @test agent_lease_usable(agent.ledger,lease.fence)
            before=control.science.lanes[string(id)].revision
            @test request("POST",path,mutation("alice"),body).status==202
            @test control.science.lanes[string(id)].revision==before
            @test length(control.science.lanes[string(id)].mutations)==1

            # The production root agent owns a separate durable consumer. It
            # must execute the exact process prepared through the protected API.
            @test timedwait(()->agent.jobs.state==:online,5)==:ok
            @test agent.jobs.connection.connection!==agent.science.connection.connection
            @test agent.jobs.connection.connection!==agent.link.control.connection
            @test timedwait(()->control.jobs.state==:online,5)==:ok
            @test control.jobs.connection.connection!==control.science.connection.connection
            submission_id=string(uuid4())
            inputs=JSON3.write((operation="fixture.echo",parameters=Dict("value"=>3),request_id=submission_id))
            @test request("POST",jobs_path,headers("alice"),inputs).status==403
            @test request("POST",jobs_path,mutation("alice"),JSON3.write((operation="fixture.echo",
                parameters=Dict("value"=>3),request_id=submission_id,executor_id=string(uuid4())))).status==400
            submitted=request("POST",jobs_path,mutation("alice"),inputs)
            @test submitted.status==202
            public_job=JSON3.read(submitted.body)
            job=get_job(store,alice,UUID(public_job.id)).job
            job_path="/runtime/api/jobs/$(public_job.id)"
            @test request("GET",job_path,headers("bob")).status==404
            @test request("GET",job_path*"/result",headers("bob")).status==404
            @test request("GET","/runtime/api/runs/$(run.id)/jobs",headers("bob")).status==404
            @test request("POST",job_path*"/cancel",mutation("bob"),"{").status==404
            @test timedwait(()->snapshot(job_path).state=="succeeded",12;pollint=0.1)==:ok
            reply=snapshot(job_path*"/result")
            outcome=P.decode_runtime_message(P.AssignedResult,JSON3.write(reply.result))
            @test outcome.result.failure===nothing
            @test outcome.result.inline_result!==nothing && outcome.result.inline_result["value"]==3
            @test outcome.result.schema_version=="1.2" && outcome.execution==job.execution
            @test outcome.result.engine_version=="environment-sha256:"*fingerprint.digest
            @test reply.job.current_assignment
            replay=request("POST",jobs_path,mutation("alice"),inputs)
            @test replay.status==202 && JSON3.read(replay.body).id==public_job.id
            @test get_job(store,alice,UUID(public_job.id)).job==job # deadline and target are not regenerated
            changed=JSON3.write((operation="fixture.echo",parameters=Dict("value"=>99),request_id=submission_id))
            @test request("POST",jobs_path,mutation("alice"),changed).status==409
            @test length(snapshot("/runtime/api/runs/$(run.id)/jobs"))==1
            @test request("GET","/health",Pair{String,String}[]).status==200

            # An operator's wrong storage identity fails the job explicitly,
            # without misreporting successful upload or losing the prepared child.
            @test timedwait(()->isempty(agent.jobs.flights),5)==:ok
            @test timedwait(()->snapshot(path).preparation=="ready",8;pollint=0.1)==:ok
            writer=agent.jobs.artifacts
            agent.jobs.artifacts=artifact_location("coordinator") # deliberately read-only
            try
                failed_reply=request("POST",jobs_path,mutation("alice"),JSON3.write((operation="fixture.large",
                    parameters=Dict("count"=>60000),request_id=string(uuid4()))))
                @test failed_reply.status==202
                failed_path="/runtime/api/jobs/$(JSON3.read(failed_reply.body).id)"
                @test timedwait(()->snapshot(failed_path).state=="failed",20;pollint=0.1)==:ok
                failed_result=snapshot(failed_path*"/result").result.result
                @test failed_result.failure.category=="artifact_unavailable"
                @test failed_result.artifact===nothing && failed_result.inline_result===nothing
                @test request("GET",failed_path*"/artifact").status==404
            finally
                agent.jobs.artifacts=writer
            end
            @test snapshot(path).executor_id==ready.executor_id

            # Payload exceeds the broker's 256 KiB message limit. The real agent
            # stores it over TLS with a write-only identity; the coordinator reads
            # through another identity, without a shared artifact filesystem.
            @test timedwait(()->snapshot(path).preparation=="ready",8;pollint=0.1)==:ok
            large_reply=request("POST",jobs_path,mutation("alice"),JSON3.write((operation="fixture.large",
                parameters=Dict("count"=>60000),request_id=string(uuid4()))))
            @test large_reply.status==202
            large_path="/runtime/api/jobs/$(JSON3.read(large_reply.body).id)"
            @test timedwait(()->snapshot(large_path).state=="succeeded",20;pollint=0.1)==:ok
            large_result=snapshot(large_path*"/result").result.result
            @test large_result.inline_result===nothing && large_result.artifact.size>262144
            @test request("GET",large_path*"/artifact",headers("bob")).status==404
            artifact_response=request("GET",large_path*"/artifact")
            @test artifact_response.status==200
            @test HTTP.header(artifact_response,"Cache-Control")=="no-store"
            @test HTTP.header(artifact_response,"X-Content-Type-Options")=="nosniff"
            @test HTTP.header(artifact_response,"Content-Disposition")=="attachment; filename=\"result.json\""
            values=JSON3.read(artifact_response.body).values
            @test length(values)==60000 && first(values)==1 && last(values)==60000
            @test bytes2hex(RT.SHA.sha256(artifact_response.body))==large_result.artifact.sha256
            @test request("GET",large_result.artifact.retrieval_reference).status==404
            @test request("GET",job_path*"/artifact").status==404 # inline result has no artifact
            @test request("HEAD",large_path*"/artifact").status==200
            @test !isdir(joinpath(dir,"agent")) || isempty(readdir(joinpath(dir,"agent"))) # no local artifact staging

            @test timedwait(()->snapshot(path).preparation=="ready",8;pollint=0.1)==:ok
            delayed_reply=request("POST",jobs_path,mutation("alice"),JSON3.write((operation="fixture.delay",
                parameters=Dict("seconds"=>10),request_id=string(uuid4()))))
            @test delayed_reply.status==202
            delayed=get_job(store,alice,UUID(JSON3.read(delayed_reply.body).id)).job
            @test timedwait(()->scientific_status(resources,lease.fence).phase==:executing,5)==:ok
            sequence=control.inventory.presence["worker-a"].report.sequence
            @test timedwait(()->control.inventory.presence["worker-a"].report.sequence>sequence,5)==:ok
            @test request("GET","/health",Pair{String,String}[]).status==200
            delayed_path="/runtime/api/jobs/$(delayed.request.job_id)"
            cancel_job=JSON3.write((request_id=string(uuid4()),))
            @test request("POST",delayed_path*"/cancel",headers("alice"),cancel_job).status==403
            @test request("POST",delayed_path*"/cancel",mutation("alice"),cancel_job).status==202
            @test timedwait(()->snapshot(delayed_path).state=="canceled",12;pollint=0.1)==:ok
            @test assigned_result(jobs,lease.fence,delayed.request.job_id).result.failure.category=="canceled"
            @test scientific_status(resources,lease.fence).preparation!=:ready

            slow_id=string(uuid4())
            slow=JSON3.write((action="prepare",parameters=Dict("seconds"=>10),request_id=slow_id))
            @test request("POST",path,mutation("alice"),slow).status==202
            @test timedwait(()->snapshot(path).preparation=="preparing",8;pollint=0.1)==:ok
            @test snapshot(path).current_request_id==slow_id
            sequence=control.inventory.presence["worker-a"].report.sequence
            @test timedwait(()->control.inventory.presence["worker-a"].report.sequence>sequence,5)==:ok
            @test request("GET","/health",Pair{String,String}[]).status==200
            cancel=JSON3.write((action="cancel",target_id=slow_id,request_id=string(uuid4())))
            @test request("POST",path,mutation("alice"),cancel).status==202
            @test timedwait(()->snapshot(path).failure=="canceled",10;pollint=0.1)==:ok
            @test snapshot(path).preparation!="ready"
            @test assignment_usable(control.coordinator,alice,id)

            # Reserve and cancel under the coordinator's short admission lock to
            # deterministically precede publication, without a special pause API.
            queue_prepare=JSON3.write((action="prepare",parameters=Dict{String,Any}(),request_id=string(uuid4())))
            @test request("POST",path,mutation("alice"),queue_prepare).status==202
            @test timedwait(()->snapshot(path).preparation=="ready",35;pollint=0.1)==:ok
            queued=lock(control.jobs.lock) do
                receipt=submit_job!(control.jobs,alice,id,"fixture.echo",Dict("value"=>42))
                cancel_job!(control.jobs,alice,UUID(receipt.job.request.job_id))
                receipt
            end
            queued_path="/runtime/api/jobs/$(queued.job.request.job_id)"
            @test timedwait(()->snapshot(queued_path).state=="canceled",12;pollint=0.1)==:ok
            @test snapshot(queued_path).cancel_requested && snapshot(queued_path).cancel_acknowledged
            canceled_outcome=snapshot(queued_path*"/result").result
            @test canceled_outcome.result.failure.category=="canceled"
            @test resources.handles[lease.fence.lease_id].request_id!=queued.job.request.job_id
            @test scientific_status(resources,lease.fence).preparation==:ready
            events=snapshot("/runtime/api/control/events").records
            @test any(e->e.job_id==public_job.id && e.stage=="succeeded" && e.executor_id==job.execution.executor_id,events)
            @test all(e->!haskey(e,:parameters) && !haskey(e,:owner),events)

            # Losing only durable transport does not stop the control/UI path.
            close(control.jobs.connection)
            @test request("GET",job_path*"/result").status==503
            @test request("GET",job_path*"/result",headers("bob")).status==404
            @test request("GET","/health",Pair{String,String}[]).status==200

            # Losing only science transport cannot stop lease/heartbeat control.
            close(agent.science.connection)
            @test timedwait(()->agent.science.state==:offline,5)==:ok
            sequence=control.inventory.presence["worker-a"].report.sequence
            @test timedwait(()->control.inventory.presence["worker-a"].report.sequence>sequence,5)==:ok
            @test request("GET","/health",Pair{String,String}[]).status==200
            @test isempty(supervisor.handles)
            release_assignment!(control.coordinator,alice,id)
            @test timedwait(()->get_assignment(store,alice,id).state==:released,10)==:ok
            @test request("GET",path).status==409
            @test isempty(driver.processes)
        finally
            close(server); close(control); close(agent); close(jobs); close(supervisor); close(store)
        end
        @test agent.cleanup_complete && driver.closed && isempty(resources.handles)
        @test istaskdone(control.science.task) && istaskdone(agent.science.task)
        @test istaskdone(agent.jobs.task) && agent.jobs.state==:stopped && isempty(agent.jobs.flights)
        @test istaskdone(control.jobs.task) && control.jobs.state==:stopped && isempty(control.jobs.flights)
    end
end
