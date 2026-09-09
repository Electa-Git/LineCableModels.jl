using Bonito, UUIDs, JSON3

@testset "shared runtime component contracts" begin
    owner = RuntimeClient(uuid4())
    @test RuntimeClient().run_id === nothing
    selector = WorkerSelector(owner, :parameters; profiles=("line-parameters",))
    @test selector.role == "parameters"
    @test selector.profiles == ("line-parameters",)
    @test_throws ArgumentError WorkerSelector(owner, "../role"; profiles=("line-parameters",))
    @test_throws ArgumentError WorkerSelector(owner, :main; profiles=())
    @test_throws ArgumentError WorkerSelector(owner, :main; profiles=("p", "p"))
    @test_throws ArgumentError WorkerSelector(owner, :main; profiles=("UPPERCASE",))
    @test_throws ArgumentError PreparationStatus(owner, "bad role")
    @test PreparationStatus(owner, :main; parameters=Dict("length"=>2)).parameters["length"]==2
    @test_throws ArgumentError PreparationStatus(owner, :main; parameters=Dict("text"=>repeat("x",65536)))
    private = PreparationStatus(owner, :main; parameters=Dict("scenario"=>"private-input-fixture"))
    @test !occursin("private-input-fixture", repr(ComponentXRay.inspection(private)))
    @test_throws ArgumentError WorkerControlPanel(owner, selector, selector)
    @test_throws ArgumentError WorkerControlPanel(RuntimeClient(uuid4()), selector)
    job = ScientificJob(owner, :parameters, "system.echo"; parameters=Dict("scenario"=>"private-input-fixture"))
    @test job.result[] === nothing
    @test !occursin("private-input-fixture", repr(job))
    @test !occursin("private-input-fixture", repr(ComponentXRay.inspection(job)))
    @test_throws ArgumentError ScientificJob(owner, "../role", "system.echo")
    @test_throws ArgumentError ScientificJob(owner, :main, "Core.eval")
    @test_throws ArgumentError ScientificJob(owner, :main, "system.echo"; parameters=[1,2])
    @test_throws ArgumentError ScientificJob(owner, :main, "system.echo"; parameters=Dict("value"=>NaN))
    @test_throws ArgumentError ScientificJob(owner, :main, "system.echo"; parameters=Dict("text"=>repeat("x",65536)))
    controls = (selector, PreparationStatus(owner, :parameters), WorkerDiagnostics(owner), WorkerControlPanel(owner, selector), job)
    for control in controls
        descriptor = ComponentXRay.inspection(control)
        @test isempty(descriptor.bindings)
        @test !isempty(descriptor.parameters)
        @test ".lc-runtime-controls" in descriptor.css_scopes
        @test descriptor.source.file == "src/widgets/RuntimeControls.jl"
        session = Session()
        try
            # jsrender registers assets and browser hooks but opens no runtime
            # connection, broker or executor in this Julia process.
            @test Bonito.jsrender(session, control) !== nothing
        finally
            close(session)
        end
    end
    @test first.(LineCableModelsPlayground.WIDGET_ROUTES) |> routes -> "/widgets/runtime-controls" in routes
    css = LineCableModelsPlayground.RUNTIME_CONTROLS_CSS
    @test !occursin(r"#[0-9a-fA-F]{3,8}\b", css)
    @test !occursin(r"--lc-[a-z-]+\s*:", css)
end

@testset "session-owned scientific display projections" begin
    LCM = LineCableModelsPlayground
    owner = RuntimeClient(uuid4())
    inputs = Observable{Any}(Dict("value"=>3))
    job = ScientificJob(owner, :main, "system.echo"; parameters=inputs)
    @test job.parameters === inputs
    fence = AssignmentFence(string(uuid4()), string(owner.run_id), "alice", "main",
        "worker-a", string(uuid4()), string(uuid4()), "fixture", "1.0.0", repeat("a",64), 1)
    request = new_job_request(job.operation, inputs[]; session_id=fence.run_id)
    result = JobResult("1.0", request.job_id, request.operation, "1.0", request.input_hash,
        "fixture", fence.fingerprint, fence.worker_id, "miss", utc_timestamp(), utc_timestamp(),
        Dict{String,Any}("private-display-value"=>3), nothing, nothing, String[])
    envelope = AssignedResult("2.0", fence, result, PreparedExecution(string(uuid4()),1,repeat("b",64)))
    provenance = JSON3.read(encode_message(envelope), Dict{String,Any})
    draft = string(uuid4())
    packet = Dict("draft_id"=>draft, "current"=>true, "provenance"=>provenance,
        "value"=>Dict("private-display-value"=>3))
    apply() = LCM.apply_job_projection!(job, JSON3.write(packet), draft)
    @test !LCM.apply_job_projection!(job, "not json", draft)
    @test !LCM.apply_job_projection!(job, JSON3.write(packet), string(uuid4()))
    @test apply()
    @test job.result[].current
    @test job.result[].provenance == envelope
    @test job.result[].value["private-display-value"] == 3
    @test !occursin("private-display-value", repr(job.result[]))
    @test !occursin("private-display-value", repr(ComponentXRay.inspection(job)))
    previous = job.result[]
    packet["unexpected"] = true
    @test !apply()
    delete!(packet, "unexpected")
    packet["current"] = 1
    @test !apply()
    packet["current"] = true
    packet["value"] = [1,2]
    @test !apply()
    packet["value"] = Dict("private-display-value"=>3)
    provenance["fence"]["run_id"] = string(uuid4())
    @test !apply()
    provenance["fence"]["run_id"] = string(owner.run_id)
    provenance["fence"]["role"] = "foreign"
    @test !apply()
    provenance["fence"]["role"] = "main"
    @test job.result[] === previous
    inputs[] = Dict("value"=>4)
    @test !apply() # Browser projection does not match the current Julia inputs.
    packet["current"] = false
    @test apply() # Retain and invalidate the same previous successful result.
    @test !job.result[].current
    @test job.result[].value === previous.value
    job.result[] = nothing
    @test !apply() # Cannot introduce a new historical result through invalidation.
    inputs[] = Dict("value"=>3)
    packet["current"] = true
    @test apply()
    session = Session()
    try
        Bonito.jsrender(session, job)
        @test !job.result[].current # Initial render establishes a fresh draft.
        @test apply()
        inputs[] = Dict("value"=>3)
        @test job.result[].current # An unchanged input echo is not a new draft.
        inputs[] = Dict("value"=>5)
        @test !job.result[].current
        inputs[] = Dict("value"=>NaN)
        @test !job.result[].current
        @test ComponentXRay.inspection(job) !== nothing # Invalid inputs stay private.
    finally
        close(session)
    end
    @test isempty(inputs.listeners)
end
