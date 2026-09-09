using LineCableModelsRuntime, UUIDs
const RT = LineCableModelsRuntime
database, owner, run_id, boot_id, barrier = ARGS
store = RuntimeStore(database)
try
    profiles, applications = ProfileRegistry(), ApplicationRegistry()
    register!(profiles, ProfileDefinition("line-parameters", "/approved/project", repeat("a", 64);
        operations=("system.echo",)))
    register!(applications, ApplicationDefinition("study", "Study", :workbench, "/study";
        requirements=Tuple(RuntimeRequirement(role, ("line-parameters",)) for role in ("main", "aux", "third"))))
    inventory = WorkerInventory(store, profiles; clock=()->0.0)
    probe = probe_worker!(inventory, "worker-a")
    accept_report!(inventory, "worker-a", RT.Protocol.WorkerAnnouncement("2.0", "worker-a",
        boot_id, inventory.coordinator_id, probe.challenge, 1, 1,
        [RT.Protocol.ProfileAdvertisement("line-parameters", "1.0.0", repeat("a", 64))]))
    manager = AssignmentManager(inventory, applications)
    write(joinpath(barrier, owner * ".ready"), "ready")
    timedwait(() -> isfile(joinpath(barrier, "go")), 45; pollint=0.01) == :ok ||
        error("allocation race barrier timed out")
    result = try
        reserve_assignment!(manager, Principal(owner), UUID(run_id), "main", "line-parameters";
            placement=PinnedPlacement("worker-a"))
        "reserved"
    catch error
        error isa AccessDenied && error.status == 409 ? "unavailable" : rethrow()
    end
    write(joinpath(barrier, owner * ".result"), result)
finally
    close(store)
end

