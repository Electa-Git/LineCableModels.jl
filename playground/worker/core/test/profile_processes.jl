using Test, Dates, TOML, LineCableModelsExecutionCore
const E = LineCableModelsExecutionCore
const worker_root = normpath(joinpath(@__DIR__, "..", ".."))

function context()
    ExecutionContext("profile-process",CancellationToken(),(_,_,_)->nothing,_->nothing,
        String[],PreparedResourceCache(),nothing)
end

function isolated_command(project, expression)
    julia = joinpath(Sys.BINDIR,Base.julia_exename())
    command = `$julia --startup-file=no --history-file=no --compiled-modules=existing --project=$project -e $expression`
    # Trusted native test, not a terminal sandbox. No inherited broker secrets.
    return setenv(`setpriv --pdeathsig KILL $command`, ["PATH"=>ENV["PATH"],
        "JULIA_DEPOT_PATH"=>join(DEPOT_PATH,':'),"JULIA_LOAD_PATH"=>"@:@stdlib",
        "JULIA_NUM_THREADS"=>"1","OPENBLAS_NUM_THREADS"=>"1"])
end

@testset "profile environments do not combine numerical engines" begin
    line = TOML.parsefile(joinpath(worker_root,"profiles","line-parameters","Manifest.toml"))["deps"]
    flow = TOML.parsefile(joinpath(worker_root,"profiles","power-flow","Manifest.toml"))["deps"]
    @test haskey(line,"LineCableModels") && !haskey(line,"PowerImpedance")
    @test haskey(flow,"PowerImpedance") && !haskey(flow,"LineCableModels")
    @test all(!haskey(environment,package) for environment in (line,flow) for package in ("Bonito","NATS","AWS","AWSS3"))
end

@testset "real line profile prepares only its own disposable process" begin
    project = joinpath(worker_root,"profiles","line-parameters")
    expression = "using LineCableModelsLineParameters; @assert !any(m -> nameof(m) in (:LineCableModels,:PowerImpedance,:NATS,:Bonito), values(Base.loaded_modules)); LineCableModelsLineParameters.main()"
    supervisor = ExecutorSupervisor(project; command=isolated_command(project,expression), startup_timeout_seconds=120)
    inputs = Dict{String,Any}("frequencies_hz"=>[50.0,100.0,150.0],"separation_m"=>0.5,
        "depth_m"=>1.0,"earth_resistivity_ohm_m"=>100.0,"line_length_m"=>1000.0)
    try
        start_executor!(supervisor,context())
        boot = supervisor.generation
        cold = @elapsed first = prepare_supervised!(supervisor,context(),inputs;timeout_seconds=300)
        warm = @elapsed second = prepare_supervised!(supervisor,context(),inputs;timeout_seconds=300)
        @test first["cache_status"] == "miss"
        @test second["cache_status"] == "hit"
        @test first["evidence"] == second["evidence"]
        @test first["evidence"]["frequencies_hz"] == [50.0,150.0]
        @test supervisor.generation == boot && process_running(supervisor.process)
        spec = OperationSpec("line.frequency_scan",identity,(_,_) -> nothing;timeout_seconds=300,execution_mode=:supervised)
        result = execute_supervised!(supervisor,spec,context(),inputs)
        @test result["frequencies_hz"] == [50.0,100.0,150.0]
        @test haskey(result,"series_impedance_ohm_per_m") && haskey(result,"shunt_admittance_s_per_m")
        foreign = OperationSpec("powerflow.prepare",identity,(_,_) -> nothing;execution_mode=:supervised)
        @test_throws PermanentOperationError execute_supervised!(supervisor,foreign,context(),Dict{String,Any}())
        stop_executor!(supervisor)
        restarted = prepare_supervised!(supervisor,context(),inputs;timeout_seconds=300)
        @test supervisor.generation == boot + 1
        @test restarted["cache_status"] == "miss"
        println("Line profile preparation seconds: cold=$cold warm=$warm (representative workload, not a compilation guarantee)")
    finally
        stop_executor!(supervisor)
    end
    @test supervisor.process === nothing
    @test !any(m -> nameof(m) in (:LineCableModels,:PowerImpedance,:NATS,:Bonito), values(Base.loaded_modules))
end

@testset "power-flow preparation does not block another profile" begin
    line_project = joinpath(worker_root,"profiles","line-parameters")
    flow_project = joinpath(worker_root,"profiles","power-flow")
    line = ExecutorSupervisor(line_project; startup_timeout_seconds=120,
        command=isolated_command(line_project,"using LineCableModelsLineParameters; LineCableModelsLineParameters.main()"))
    flow = ExecutorSupervisor(flow_project; startup_timeout_seconds=120,
        command=isolated_command(flow_project,"using LineCableModelsPowerFlow; @assert !any(m -> nameof(m) in (:LineCableModels,:PowerImpedance,:NATS,:Bonito), values(Base.loaded_modules)); LineCableModelsPowerFlow.main()"))
    inputs = Dict{String,Any}("frequencies_hz"=>[50.0],"separation_m"=>0.5,
        "depth_m"=>1.0,"earth_resistivity_ohm_m"=>100.0,"line_length_m"=>1000.0)
    specification = Dict{String,Any}("specification"=>Dict{String,Any}("earth_resistivity_ohm_m"=>100.0))
    task = nothing
    try
        prepare_supervised!(line,context(),inputs;timeout_seconds=300)
        task = @async prepare_supervised!(flow,context(),specification;timeout_seconds=600)
        yield()
        spec = OperationSpec("line.frequency_scan",identity,(_,_) -> nothing;timeout_seconds=30,execution_mode=:supervised)
        elapsed = @elapsed result = execute_supervised!(line,spec,context(),inputs)
        @test result["frequencies_hz"] == [50.0]
        @test elapsed < 30 # finite line deadline, not the power-flow preparation budget
        prepared = fetch(task)
        @test prepared["cache_status"] == "miss"
        @test prepared["evidence"]["preparation_kind"] == "solved_and_linearized_model"
        @test length(prepared["evidence"]["prepared_resource_key"]) == 64
        generation = flow.generation
        warm = @elapsed cached = prepare_supervised!(flow,context(),specification;timeout_seconds=600)
        @test cached["cache_status"] == "hit" && flow.generation == generation
        @test cached["evidence"] == prepared["evidence"]
        stop_executor!(flow)
        @test process_running(line.process)
        @test execute_supervised!(line,spec,context(),inputs)["frequencies_hz"] == [50.0]
        println("Independent line request seconds=$elapsed; repeated power-flow preparation seconds=$warm")
    finally
        stop_executor!(flow)
        stop_executor!(line)
        task === nothing || try wait(task) catch end
    end
end
