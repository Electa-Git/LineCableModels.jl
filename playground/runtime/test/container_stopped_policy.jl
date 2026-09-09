using Test, UUIDs, LineCableModelsRuntime
const StoppedPolicyRT = LineCableModelsRuntime

function test_stopped_container_policy()
    image_ref = get(ENV,"LCM_TEST_STOPPED_CONTAINER_IMAGE","")
    occursin(r"^[A-Za-z0-9][A-Za-z0-9._:/-]*@sha256:[a-f0-9]{64}$",image_ref) ||
        error("Set LCM_TEST_STOPPED_CONTAINER_IMAGE to an already cached digest-pinned image.")
    directory = mktempdir(;prefix="lcm-stopped-policy-",cleanup=false)
    runner = CommandRunner()
    journal = nothing
    passed = false
    cleaned = false
    try
        host = check_container_host(runner;requested=get(ENV,"LCM_TEST_CONTAINER_RUNTIME","auto"))
        inspect = StoppedPolicyRT.scoped_container_command(runner,host,["image","inspect","--format","{{json .}}",image_ref])
        inspect.exitcode == 0 || error("The test image must already exist locally; no pull is permitted.")
        image = StoppedPolicyRT.JSON3.read(inspect.output)
        image_ref in image.RepoDigests || error("Local digest identity differs.")
        inventory() = begin
            result = StoppedPolicyRT.scoped_container_command(runner,host,["container","ls","--all","--no-trunc","--format","{{.ID}}"])
            result.exitcode == 0 || error("Container inventory unavailable.")
            sort(split(strip(result.output),'\n'))
        end
        baseline = inventory()
        journal = ResourceJournal(joinpath(directory,"journal"),"policy-audit";capacity=1)
        scope = container_scope(runner,host)
        for kind in (:scientific,:terminal)
            profile = ProfileDefinition("policy-audit",image_ref,last(split(image_ref,"@sha256:"));kind,isolation=:container,
                operations=kind == :scientific ? ("fixture.echo",) : ())
            fence = StoppedPolicyRT.Protocol.AssignmentFence(string(uuid4()),string(uuid4()),"audit","policy",
                "policy-audit",string(uuid4()),string(uuid4()),profile.id,string(profile.version),profile.fingerprint,1)
            receipt = reserve_resource!(journal,fence,host.engine.name,scope)
            policy = ContainerPolicy(profile,receipt)
            # Deliberately create-only: the base test image need not contain LCM
            # entry files. This does NOT pass image admission or start a process.
            result = StoppedPolicyRT.scoped_container_command(runner,host,container_create_arguments(policy))
            result.exitcode == 0 || error("Stopped policy acquisition failed.")
            id = strip(result.output)
            inspected = StoppedPolicyRT.scoped_container_command(runner,host,["container","inspect","--format","{{json .}}",id])
            inspected.exitcode == 0 || error("Stopped policy inspection failed.")
            object = StoppedPolicyRT.JSON3.read(inspected.output)
            evidence_path = joinpath(directory,"$kind-inspection.json")
            open(evidence_path,"w") do io
                chmod(evidence_path,0o600)
                write(io,inspected.output)
            end
            @testset "actual stopped $kind policy on $(host.engine.name)" begin
                @test object.State.Running === false
                @test object.State.Pid == 0
                @test StoppedPolicyRT.verify_created_container(policy,object,
                    (;id=image.Id,environment=StoppedPolicyRT.container_env_map(image.Config.Env))) === nothing
                @test remove_owned_container!(journal,runner,host,receipt)
                @test isempty(resource_receipts(journal))
                @test inventory() == baseline
            end
        end
        passed = true
    finally
        try
            if journal !== nothing
                host = check_container_host(runner;requested=get(ENV,"LCM_TEST_CONTAINER_RUNTIME","auto"))
                recover_containers!(journal,runner,host)
                cleaned = isempty(resource_receipts(journal))
            else
                cleaned = true
            end
        finally
            journal === nothing || close(journal)
            close(runner)
            if passed && cleaned
                realpath(directory) == directory && startswith(basename(directory),"lcm-stopped-policy-") ||
                    error("Refusing unexpected diagnostic cleanup.")
                rm(directory;recursive=true)
            else
                println(stderr,"Stopped policy diagnostics retained: ",directory)
            end
        end
    end
    println("Both newly owned stopped containers removed; no process started, no effective-isolation claim.")
end

test_stopped_container_policy()
