@testset "durable agent retries recover results without repeating uncertain computation" begin
    mktempdir() do dir
        # Finite test-only child driver. Real TLS, stream and job owner; not an
        # assertion that the host has admitted mandatory numerical OS quotas.
        project=normpath(joinpath(@__DIR__,"..","..","worker","core"))
        fingerprint=native_environment_fingerprint(project)
        profiles=ProfileRegistry()
        register!(profiles,ProfileDefinition("fixture",project,fingerprint.digest;
            operations=("fixture.echo","fixture.delay","fixture.evict","fixture.fail"),
            budget=ResourceBudget(prepare_seconds=30,job_seconds=60)))
        driver=ScientificDriverFixture(profiles,Dict("fixture"=>fixture_command(project,
            joinpath(@__DIR__,"scientific_child.jl"),"300")),
            Dict{String,Tuple{P.AssignmentFence,RT.ExecutionCore.ExecutorSupervisor}}(),0,0,false,true,false)
        resources=ScientificResources(driver)
        agent=AgentService(AgentConfig("worker-a",endpoint("worker-a"),profiles,joinpath(dir,"owned")),resources)
        coordinator=BrokerJobs(endpoint("coordinator"),CoordinatorIdentity())
        jobs=agent.jobs
        jobs.connection=BrokerJobs(endpoint("worker-a"),WorkerIdentity("worker-a"))
        control_id=string(uuid4())
        receive_probe!(agent.ledger,P.WorkerProbe("2.0","worker-a",control_id,string(uuid4())))
        fence=P.AssignmentFence(string(uuid4()),string(uuid4()),"alice","main","worker-a",
            agent.ledger.boot_id,control_id,"fixture","1.0.0",fingerprint.digest,1)
        @test handle_lease_control!(agent.ledger,P.LeaseControl("2.0",string(uuid4()),"grant",fence,1,60000)).accepted
        # This test manually controls pulls to reproduce lost acknowledgements.
        # Supply the independent two-second coordinator presence that the root
        # scheduler exercises over TLS in science_tls.jl; never freeze its clock.
        probing=Ref(true)
        probes=@async while probing[]
            receive_probe!(agent.ledger,P.WorkerProbe("2.0","worker-a",control_id,string(uuid4()))) ||
                error("fixture coordinator presence was rejected")
            for _ in 1:40
                probing[] || break
                sleep(0.05)
            end
        end
        function pull(job)
            NATS.JetStream.stream_publish(coordinator.connection,P.assigned_job_subject(fence),P.encode_message(job))
            delivery=poll_assigned_job!(jobs.connection,agent.ledger)
            @test delivery isa AssignedDelivery && delivery.job==job
            return delivery
        end
        function finish(delivery)
            flight=accept_assigned_delivery!(jobs,delivery)
            @test timedwait(()->istaskdone(flight.task),10;pollint=0.025)==:ok
            fetch(flight.task)
            @test flight.phase==:completed && flight.outcome!==nothing
            tick_jobs!(jobs)
            @test isempty(jobs.flights)
            return assigned_result(coordinator,fence,delivery.job.request.job_id)
        end
        try
            recover_owned!(resources)
            fetch(prepare_assigned!(resources,fence,Dict{String,Any}()))
            handle=resources.handles[fence.lease_id]
            process=handle.supervisor.process
            job=assigned_fixture_job(resources,fence)
            delivery=pull(job)
            @test RT.delivery_count(jobs.connection,delivery)==1
            @test progress_assigned_delivery!(jobs.connection,agent.ledger,delivery)===nothing
            @test assigned_result(coordinator,fence,job.request.job_id)===nothing # progress is not completion
            changed=P.AssignedJob(job.protocol_version,job.fence,
                P.new_job_request("fixture.echo",Dict("value"=>99);session_id=fence.run_id),job.execution)
            @test_throws AccessDenied accept_assigned_delivery!(jobs,AssignedDelivery(changed,delivery.message))
            @test isempty(jobs.flights) && handle.request_id!=job.request.job_id
            malformed=NATS.Msg(delivery.message.subject,delivery.message.sid,delivery.message.reply_to,0,Vector{UInt8}("{"))
            @test_throws AccessDenied accept_assigned_delivery!(jobs,AssignedDelivery(job,malformed))
            badreply=NATS.Msg(delivery.message.subject,delivery.message.sid,
                "\$JS.ACK.$(RT.job_stream("worker-a")).agent.bad.1.1.1.0",0,Vector{UInt8}(P.encode_message(job)))
            @test_throws AccessDenied accept_assigned_delivery!(jobs,AssignedDelivery(job,badreply))
            outcome=finish(delivery)
            @test outcome.result.inline_result["value"]==3 && outcome.result.schema_version=="1.2"
            @test handle.supervisor.process===process

            # Hold a real read-only inspection in progress while accepting a
            # job. The bounded admission wait cannot reject valid prepared work.
            probe_job=assigned_fixture_job(resources,fence)
            probe_delivery=pull(probe_job)
            pending_flight=lock(handle.supervisor.lock) do
                inspection=refresh_preparation!(resources,fence)
                yield()
                flight=accept_assigned_delivery!(jobs,probe_delivery)
                @test timedwait(()->flight.phase==:admitting,3)==:ok
                @test !istaskdone(inspection) && !istaskdone(flight.task)
                flight
            end
            @test timedwait(()->istaskdone(pending_flight.task),8)==:ok
            @test pending_flight.phase==:completed
            @test assigned_result(coordinator,fence,probe_job.request.job_id).result.failure===nothing
            tick_jobs!(jobs)

            # Computation and durable result exist, but the input ACK is lost.
            # Redelivery must use that result, not the supplied calculation.
            saved_job=assigned_fixture_job(resources,fence,"fixture.echo",Dict("value"=>22))
            original=pull(saved_job)
            saved=RT.job_outcome(saved_job,P.utc_timestamp(),ScientificOutput(Dict{String,Any}("value"=>22),"1.2",String[]))
            subject=P.assigned_result_subject(fence,saved_job.request.job_id)
            RT.job_request(jobs.connection,NATS.JetStream.PubAck,subject,
                (P.encode_message(saved),["Nats-Msg-Id"=>subject]))
            RT.job_request(jobs.connection,NATS.Msg,original.message.reply_to,"-NAK")
            again=poll_assigned_job!(jobs.connection,agent.ledger)
            @test RT.delivery_count(jobs.connection,again)==2
            previous=handle.request_id
            @test finish(again)==saved
            @test handle.request_id==previous # no child request, even an inspection

            # No saved result means previous execution is uncertain, never an
            # instruction to repeat it automatically.
            uncertain=assigned_fixture_job(resources,fence)
            original=pull(uncertain)
            RT.job_request(jobs.connection,NATS.Msg,original.message.reply_to,"-NAK")
            again=poll_assigned_job!(jobs.connection,agent.ledger)
            @test RT.delivery_count(jobs.connection,again)==2
            failure=finish(again)
            @test failure.result.failure.category=="execution_uncertain"
            @test !failure.result.failure.retryable && handle.request_id==previous

            target=P.PreparedExecution(job.execution.executor_id,job.execution.executor_generation+1,job.execution.preparation_key)
            stale=P.AssignedJob("2.0",fence,P.new_job_request("fixture.echo",Dict("value"=>1);session_id=fence.run_id),target)
            failure=finish(pull(stale))
            @test failure.result.failure.category=="not_admitted" && handle.request_id==previous
            rejected=assigned_fixture_job(resources,fence,"fixture.fail",Dict{String,Any}())
            failure=finish(pull(rejected))
            @test failure.result.failure.category=="operation_rejected"
            @test !occursin("/private/fixture",P.encode_message(failure))
            @test scientific_status(resources,fence).preparation==:ready

            # Cross the provisioned 30-second AckWait using the real heartbeat.
            # A second pull still must not receive a concurrent calculation.
            @test handle_lease_control!(agent.ledger,P.LeaseControl("2.0",string(uuid4()),"renew",fence,2,60000)).accepted
            long=assigned_fixture_job(resources,fence,"fixture.delay",Dict("seconds"=>32);timeout=Second(60))
            delivery=pull(long)
            flight=accept_assigned_delivery!(jobs,delivery)
            @test accept_assigned_delivery!(jobs,delivery)===flight
            @test timedwait(()->scientific_status(resources,fence).phase==:executing,5)==:ok
            started=time_ns()
            while (time_ns()-started)/1e9<31 && !istaskdone(flight.task)
                sleep(0.2)
            end
            @test !istaskdone(flight.task)
            @test poll_assigned_job!(jobs.connection,agent.ledger)===nothing
            @test timedwait(()->istaskdone(flight.task),8)==:ok
            @test flight.phase==:completed
            finished=assigned_result(coordinator,fence,long.request.job_id)
            @test finished!==nothing && finished.result.inline_result!==nothing && finished.result.inline_result["done"]
            @test process_running(process)
            tick_jobs!(jobs)
            stopping=assigned_fixture_job(resources,fence,"fixture.delay",Dict("seconds"=>20))
            flight=accept_assigned_delivery!(jobs,pull(stopping))
            @test timedwait(()->scientific_status(resources,fence).phase==:executing,5)==:ok
            close(jobs) # independent stop cancels only this owner's exact job
            @test istaskdone(flight.task) && !process_running(process)
            @test agent_lease_usable(agent.ledger,fence) # closing job transport is not lease release
        finally
            probing[]=false
            wait(probes)
            close(agent); close(coordinator)
        end
        @test jobs.closed && jobs.state==:stopped && isempty(jobs.flights)
        @test agent.cleanup_complete && isempty(driver.processes) && isempty(resources.handles)
    end
end
