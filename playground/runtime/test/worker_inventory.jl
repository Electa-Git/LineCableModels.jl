@testset "challenged presence is not registration or remembered warmth" begin
    mktempdir() do directory
        store = RuntimeStore(joinpath(directory, "runtime.sqlite"))
        profiles = ProfileRegistry()
        register!(profiles, ProfileDefinition("line-parameters", "/approved/project", repeat("a", 64);
            operations=("line.evaluate",)))
        admin = Principal("operator"; administrator=true)
        alice = Principal("alice")
        clock = Ref(0.0)
        inventory = WorkerInventory(store, profiles; clock=()->clock[])
        try
            trust = WorkerTrust("worker-a", "credential-a", ("line-parameters",); capacity=2)
            enroll_worker!(store, admin, trust)
            @test only(worker_inventory(inventory, alice)).liveness == "unknown"
            @test_throws AccessDenied probe_worker!(inventory, "worker-a")
            set_registration_state!(store, admin, "worker-a", :approved; expected_revision=1)
            boot = string(uuid4())
            installed = [RT.Protocol.ProfileAdvertisement("line-parameters", "1.0.0", repeat("a", 64))]
            report(probe; id="worker-a", boot_id=boot, sequence=1, capacity=2,
                advertised=installed, coordinator=inventory.coordinator_id) =
                RT.Protocol.WorkerAnnouncement("2.0", id, boot_id, coordinator,
                    probe.challenge, sequence, capacity, advertised)
            probe = probe_worker!(inventory, "worker-a")
            @test_throws AccessDenied accept_report!(inventory, "worker-b", report(probe))
            @test_throws AccessDenied accept_report!(inventory, "worker-a", report(probe; capacity=3))
            @test_throws AccessDenied accept_report!(inventory, "worker-a",
                report(probe; advertised=[RT.Protocol.ProfileAdvertisement("line-parameters", "1.0.0", repeat("b", 64))]))
            @test_throws AccessDenied accept_report!(inventory, "worker-a", report(probe; coordinator=string(uuid4())))
            accepted = accept_report!(inventory, "worker-a", report(probe))
            @test accepted.report.boot_id == boot
            public_report = only(worker_inventory(inventory, alice)).report
            @test public_report.boot_id == boot
            @test !haskey(public_report, :challenge) && !haskey(public_report, :coordinator_id)
            @test only(worker_inventory(inventory, alice)).liveness == "online"
            @test_throws AccessDenied accept_report!(inventory, "worker-a", report(probe)) # challenge consumed
            clock[] = 5
            @test only(worker_inventory(inventory, alice)).liveness == "stale"
            clock[] = 11
            @test only(worker_inventory(inventory, alice)).liveness == "offline"
            delayed = probe_worker!(inventory, "worker-a")
            clock[] = 22
            @test_throws AccessDenied accept_report!(inventory, "worker-a", report(delayed; sequence=2))
            @test only(worker_inventory(inventory, alice)).liveness == "offline"
            fresh = probe_worker!(inventory, "worker-a")
            @test_throws AccessDenied accept_report!(inventory, "worker-a", report(fresh))
            @test accept_report!(inventory, "worker-a", report(fresh; sequence=2)).report.sequence == 2
            replacement = probe_worker!(inventory, "worker-a")
            new_boot = string(uuid4())
            @test accept_report!(inventory, "worker-a", report(replacement; boot_id=new_boot)).report.boot_id == new_boot
            retired = probe_worker!(inventory, "worker-a")
            @test_throws AccessDenied accept_report!(inventory, "worker-a", report(retired; sequence=3))
            restarted = WorkerInventory(store, profiles; clock=()->clock[])
            @test only(worker_inventory(restarted, alice)).liveness == "unknown"
            @test_throws AccessDenied accept_report!(restarted, "worker-a", report(retired))
            set_registration_state!(store, admin, "worker-a", :draining; expected_revision=2)
            row = only(worker_inventory(inventory, alice))
            @test row.registration.state == "draining" && row.liveness == "online"
            @test !haskey(row.registration, :prepared)
            limited = WorkerInventory(store, profiles; max_workers=1, clock=()->clock[])
            probe_worker!(limited, "worker-a")
            enroll_worker!(store, admin, WorkerTrust("worker-b", "credential-b", ("line-parameters",)))
            set_registration_state!(store, admin, "worker-b", :approved; expected_revision=1)
            @test_throws CapacityUnavailable probe_worker!(limited, "worker-b")
        finally
            close(store)
        end
    end
end
