import HTTP, JSON3

@testset "shared worker-control API enforces configuration, identity and lease ownership" begin
    mktempdir() do dir
        store = RuntimeStore(joinpath(dir, "runtime.sqlite"))
        applications, profiles = ApplicationRegistry(), ProfileRegistry()
        definition = ApplicationDefinition("study", "Study", :workbench, "/study";
            requirements=(RuntimeRequirement("main", ("line-parameters",)),))
        register!(applications, definition)
        register!(profiles, ProfileDefinition("line-parameters", "/approved/environment", repeat("a",64);
            operations=("system.echo",)))
        trust = WorkerTrust("worker-a", "credential-a", ("line-parameters",))
        config = ControlConfig(BrokerEndpoint("tls://broker.invalid", "/private/secret"), profiles, [trust])
        clock = Ref(0.0)
        service = ControlService(config, store, applications; clock=()->clock[])
        supervisor = UIHostSupervisor(store, applications, joinpath(dir, "hosts"))
        key = repeat("proxy-private-", 4)
        policy = ProxyIdentity("https://lcm.test", ["127.0.0.1"], key; administrators=["operator"])
        server = start_gateway(supervisor, policy; control=service)
        @test isempty(supervisor.handles) && service.task === service.connector === nothing
        @test isempty(list_runs(store, Principal("operator"; administrator=true)))
        @test isempty(list_assignments(store, Principal("operator"; administrator=true)))
        base = "http://127.0.0.1:$(HTTP.port(server))"
        headers(user) = ["X-LCM-Proxy-Key"=>key, "X-LCM-Principal"=>user]
        mutation(user) = [headers(user); "Origin"=>"https://lcm.test"; "X-LCM-Request"=>"1"]
        request(method, path, hs=Pair{String,String}[], body=nothing) =
            HTTP.request(method, base * path, hs, body; proxy=nothing,
                retry=false, status_exception=false, request_timeout=10)
        payload(x) = JSON3.write(x)
        alice, bob = Principal("alice"), Principal("bob")
        try
            @test JSON3.read(request("GET", "/runtime/api/capabilities").body).worker_control
            @test JSON3.read(request("GET", "/runtime/api/capabilities").body).preparation_control
            capabilities = JSON3.read(request("GET", "/runtime/api/capabilities").body)
            @test capabilities.assigned_execution && capabilities.private_terminal
            @test !RT.control_capabilities(nothing).private_terminal
            @test !RT.control_snapshot(nothing, alice).assigned_execution
            @test request("GET", "/runtime/api/control").status == 401
            @test request("GET", "/runtime/api/control", ["X-LCM-Principal"=>"operator"]).status == 401
            @test request("GET", "/health").status == 200
            client_asset = request("GET", "/runtime/assets/runtime-client.js")
            @test client_asset.status == 200
            @test occursin("LineCableModelsRuntimeClient", String(copy(client_asset.body)))
            @test request("GET", "/runtime/control").status == 401
            @test request("GET", "/runtime/control", headers("alice")).status == 200
            @test request("GET", "/runtime/assets/runtime-controls.js").status == 200
            @test request("GET", "/runtime/assets/runtime-controls.css").status == 200
            snapshot = JSON3.read(request("GET", "/runtime/api/control", headers("alice")).body)
            @test snapshot.enabled && snapshot.broker == "unavailable"
            @test snapshot.assigned_execution == capabilities.assigned_execution
            @test snapshot.preparation_control == capabilities.preparation_control
            @test snapshot.private_terminal == capabilities.private_terminal
            @test isempty(snapshot.provisioned) && isempty(snapshot.workers)
            @test !snapshot.administrator && length(snapshot.profiles) == 1
            @test !occursin("/private", payload(snapshot))
            @test !occursin("/approved", payload(snapshot))
            @test !occursin("broker.invalid", payload(snapshot))
            @test isempty(supervisor.handles) && service.task === nothing
            @test length(JSON3.read(request("GET", "/runtime/api/control", headers("operator")).body).provisioned) == 1

            enrollment = payload((worker_id="worker-a", request_id=string(uuid4())))
            @test request("POST", "/runtime/api/workers", mutation("alice"), enrollment).status == 403
            @test request("POST", "/runtime/api/workers", headers("operator"), enrollment).status == 403
            @test request("POST", "/runtime/api/workers", mutation("operator"),
                payload((worker_id="not-provisioned", request_id=string(uuid4())))).status == 400
            @test request("POST", "/runtime/api/workers", mutation("operator"),
                payload((worker_id="worker-a", request_id=string(uuid4()), password="forbidden"))).status == 400
            @test request("POST", "/runtime/api/workers", mutation("operator"),
                "{\"worker_id\":\"worker-a\",\"worker_id\":\"worker-a\",\"request_id\":\"$(uuid4())\"}").status == 400
            enrolled = request("POST", "/runtime/api/workers", mutation("operator"), enrollment)
            @test enrolled.status == 200
            @test JSON3.read(enrolled.body).state == "pending"
            @test request("POST", "/runtime/api/workers", mutation("operator"), enrollment).body == enrolled.body
            worker = JSON3.read(only(JSON3.read(request("GET", "/runtime/api/workers", headers("alice")).body)).registration |> payload)
            @test worker.worker_id == "worker-a"
            approved_request = payload((state="approved", expected_revision=1, request_id=string(uuid4())))
            @test request("PATCH", "/runtime/api/workers/worker-a", mutation("alice"), approved_request).status == 403
            @test request("PATCH", "/runtime/api/workers/worker-a", mutation("operator"),
                payload((state="approved", expected_revision=true, request_id=string(uuid4())))).status == 400
            approved = request("PATCH", "/runtime/api/workers/worker-a", mutation("operator"), approved_request)
            @test approved.status == 200 && JSON3.read(approved.body).revision == 2
            @test request("PATCH", "/runtime/api/workers/worker-a", mutation("operator"), approved_request).body == approved.body

            run = reserve_run!(store, alice, definition)
            foreign = reserve_run!(store, bob, definition)
            control_page = request("GET", "/runtime/control?run=$(run.id)", headers("alice"))
            @test control_page.status == 200
            @test occursin("data-lcm-runtime-controls", String(copy(control_page.body)))
            @test occursin("line-parameters", String(copy(control_page.body)))
            @test !occursin(string(foreign.id), String(copy(control_page.body)))
            @test request("GET", "/runtime/control?run=$(run.id)", headers("bob")).status == 404
            @test request("GET", "/runtime/control?run=$(run.id)&run=$(run.id)", headers("alice")).status == 400
            @test isempty(supervisor.handles) && service.task === nothing
            route = "/runtime/api/runs/$(run.id)/assignments"
            grant_body = payload((role="main", profile="line-parameters",
                placement=(mode="pinned", worker_id="worker-a"), request_id=string(uuid4())))
            @test request("GET", route, headers("bob")).status == 404
            @test request("POST", route, mutation("bob"), "{}").status == 404
            @test request("POST", route, mutation("alice"), grant_body).status == 409 # no report, not false ready
            agent = AgentLeaseLedger("worker-a", profiles; clock=()->clock[])
            probe = probe_worker!(service.inventory, "worker-a")
            @test receive_probe!(agent, probe)
            report = RT.Protocol.WorkerAnnouncement("2.0", "worker-a", agent.boot_id,
                probe.coordinator_id, probe.challenge, 1, 1,
                [RT.Protocol.ProfileAdvertisement("line-parameters", "1.0.0", repeat("a",64))])
            RT.accept_control_record!(service, RT.ControlEnvelope("worker-a", report))
            response = request("POST", route, mutation("alice"), grant_body)
            @test response.status == 202
            assignment = JSON3.read(response.body)
            id = UUID(assignment.id)
            @test assignment.state == "reserving" && !assignment.usable
            @test assignment.preparation == "unknown"
            @test !occursin(probe.coordinator_id, String(copy(response.body)))
            @test request("POST", route, mutation("alice"), grant_body).body == response.body
            @test length(list_assignments(store, alice)) == 1
            ack = handle_lease_control!(agent, service.coordinator.flights[id].command)
            RT.accept_control_record!(service, RT.ControlEnvelope("worker-a", ack))
            @test JSON3.read(request("GET", "/runtime/api/assignments/$id", headers("alice")).body).usable
            @test request("GET", "/runtime/api/assignments/$id", headers("bob")).status == 404
            @test request("DELETE", "/runtime/api/assignments/$id", mutation("bob"), "{}").status == 404
            @test only(JSON3.read(request("GET", route, headers("alice")).body)).state == "active"
            science_route="/runtime/api/assignments/$id/science"
            prepare_body=payload((action="prepare",parameters=Dict(),request_id=string(uuid4())))
            @test request("GET",science_route).status==401
            @test request("GET",science_route,headers("bob")).status==404
            @test request("POST",science_route,mutation("bob"),"not even JSON").status==404
            @test request("POST",science_route,headers("alice"),prepare_body).status==403
            @test request("POST",science_route,mutation("alice"),
                payload((action="prepare",parameters=Dict(),request_id=string(uuid4()),image="forbidden"))).status==400
            @test request("POST",science_route,mutation("alice"),
                payload((action="eval",parameters=Dict(),request_id=string(uuid4())))).status==400
            unavailable=JSON3.read(request("GET",science_route,headers("alice")).body)
            @test unavailable.preparation=="unknown" && unavailable.channel=="offline"
            @test request("POST",science_route,mutation("alice"),prepare_body).status==503
            @test isempty(supervisor.handles) && service.science.connection===nothing
            @test only(JSON3.read(request("GET", "/runtime/api/workers", headers("alice")).body)).occupied == 1
            foreign_events = request("GET", "/runtime/api/control/events", headers("bob"))
            @test !occursin(string(run.id), String(copy(foreign_events.body)))
            @test !occursin(string(id), String(copy(foreign_events.body)))
            @test request("GET", "/runtime/api/control/events?after=bad", headers("alice")).status == 400
            @test request("GET", "/runtime/api/control/events?after=0&after=1", headers("alice")).status == 400

            drained = request("PATCH", "/runtime/api/workers/worker-a", mutation("operator"),
                payload((state="draining", expected_revision=2, request_id=string(uuid4()))))
            @test drained.status == 200
            # Replaying a prior approval returns that reply, without undoing drain.
            @test request("PATCH", "/runtime/api/workers/worker-a", mutation("operator"), approved_request).body == approved.body
            @test only(list_registrations(store, alice)).state == :draining
            @test !assignment_usable(service.coordinator, alice, id) # rejects new starts
            released = request("DELETE", "/runtime/api/assignments/$id", mutation("alice"),
                payload((request_id=string(uuid4()),)))
            @test released.status == 200 && JSON3.read(released.body).state == "releasing"
            @test handle_lease_control!(agent, service.coordinator.flights[id].command) === nothing
            ack = complete_agent_cleanup!(agent, service.coordinator.flights[id].command.fence)
            RT.accept_control_record!(service, RT.ControlEnvelope("worker-a", ack))
            @test get_assignment(store, alice, id).state == :released
            @test request("GET", "/health").status == 200
            @test isempty(supervisor.handles)
            changed = ControlConfig(config.endpoint, profiles,
                [WorkerTrust("worker-a", "changed-credential", ("line-parameters",))])
            @test_throws ArgumentError ControlService(changed, store, applications)
            close(service.jobs)
            @test !JSON3.read(request("GET", "/runtime/api/capabilities").body).assigned_execution
            @test !JSON3.read(request("GET", "/runtime/api/control", headers("alice")).body).assigned_execution
            close(service)
            @test !RT.control_snapshot(service, alice).enabled
            @test !RT.control_capabilities(service).worker_control
        finally
            close(server)
            close(service); close(service)
            close(supervisor)
            close(store)
        end
        @test service.state == :stopped
    end
end
