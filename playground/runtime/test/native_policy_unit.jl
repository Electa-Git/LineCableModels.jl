using Test, UUIDs, LineCableModelsRuntime
include("managed_agent.jl")

@testset "native service directives pass local inert systemd validation" begin
    directory = mktempdir(;prefix="lcm-native-unit-",cleanup=false)
    path,config,_ = managed_config(directory;worker="audit-"*string(uuid4()))
    runner = CommandRunner()
    journal = ResourceJournal(ManagedRT.agent_journal_root(config),config.worker_id)
    try
        # Pure rendering only: this synthetic scope is never sent to a manager
        # for start, and no claim about effective host quotas is made.
        project = normpath(joinpath(@__DIR__,"..","..","worker","profiles","line-parameters"))
        environment = native_environment_fingerprint(project)
        profile = ProfileDefinition("fixture",project,environment.digest;operations=["fixture.echo"])
        fence = ManagedRT.Protocol.AssignmentFence(string(uuid4()),string(uuid4()),"audit","main",config.worker_id,
            string(uuid4()),string(uuid4()),profile.id,string(profile.version),profile.fingerprint,1)
        receipt = reserve_resource!(journal,fence,:native,repeat("a",64))
        policy = NativePolicy(profile,receipt)
        manager = ManagedAgentIdentity(ManagedRT.agent_unit_name(config.worker_id),repeat("b",32),"/private/agent",
            path,journal.root,config.worker_id)
        host = NativeHostCheck(receipt.scope,Sys.which("systemd-run"),1000,1000,())
        properties = ManagedRT.native_service_properties(policy,manager)
        unit_keys = ("Description=","BindsTo=","After=")
        unit_property(value) = any(prefix->startswith(value,prefix),unit_keys)
        text = "[Unit]\n" * join(filter(unit_property,properties),'\n') * "\n[Service]\n" *
            join(filter(!unit_property,properties),'\n') * "\nExecStart=" *
            ManagedRT.systemd_command(ManagedRT.native_exec_arguments(policy,host,environment)) *
            "\nStandardInput=null\nStandardOutput=null\nStandardError=null\n"
        generated = joinpath(directory,resource_name(receipt))
        for (file,content) in ((generated,text),(joinpath(directory,manager.unit),agent_service_unit(path)))
            open(file,"w") do io; chmod(file,0o600); write(io,content); end
        end
        command = setenv(Cmd([Sys.which("systemd-analyze"),"--user","verify",generated]),ManagedRT.container_command_environment())
        result = run_owned_command!(runner,command)
        diagnostic = joinpath(directory,"verify.log")
        open(diagnostic,"w") do io; chmod(diagnostic,0o600); write(io,result.output,result.diagnostic); end
        @test result.exitcode == 0
        @test isempty(runner.active)
        @test only(resource_receipts(journal)).physical_id === nothing
        println("Inert native service validation: ",directory)
    finally
        for receipt in resource_receipts(journal); forget_resource!(journal,receipt); end
        close(journal);close(runner)
    end
end
