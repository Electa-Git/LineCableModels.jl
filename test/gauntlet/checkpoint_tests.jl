@testitem "Gauntlet / queue and formulation checkpoints survive source edits" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using LineCableModels, JLD2, SHA, TOML
    const calls=Float64[]
    const fail=Ref(true)
    const source=Ref("")
    struct CheckpointBackend <: LineCableModels.Grammar.AbstractFormulation
        factor::Float64
    end
    Base.NamedTuple(formulation::CheckpointBackend)=(factor=formulation.factor,)
    function LineCableModels.compute(problem::LineParametersProblem, formulation::CheckpointBackend;options=(;))
        push!(calls,formulation.factor)
        # The declaration is edited while a solver is active, then disappears.
        if isfile(source[])
            write(source[],"# unrelated edit during execution\n")
            rm(source[])
        end
        fail[] && formulation.factor==3 && error("interrupted formulation")
        z=fill(complex(formulation.factor),2,2,length(problem.frequencies))
        result=LineParameters(PhaseDomain,z,z,copy(problem.frequencies))
        haskey(options,:on_result) && options.on_result(problem,1,result)
        result
    end
    model=load_case(:two_insulated_wires;variation=ExactOverrides(frequencies=[1.,37.]))
    grid=Gridspace{CheckpointBackend}(CheckpointBackend,(Grid((2.,3.,4.)),))
    events=Int[]
    candidate=BenchmarkCalculation(:candidate,model.problem,grid;
        options=(on_result=(problem,index,result)->push!(events,index),))
    reference=BenchmarkCalculation(:reference,model.problem,CheckpointBackend(1.))
    mktempdir() do parent
        source[]=joinpath(parent,"declaration.jl")
        write(source[],"# initial declaration\n")
        first_definition=benchmark_definition(:first,model.id,:fixture,source[],model,
            reference,candidate,(bands=(:all,),),(;))
        queued=benchmark_definition(:queued,model.id,:fixture,source[],model,
            reference,BenchmarkCalculation(:candidate,model.problem,CheckpointBackend(5.)),
            (bands=(:all,),),(;))
        directory=joinpath(parent,"campaign")
        @test_throws r"interrupted formulation" run_campaign(directory,[first_definition,queued];on_error=:fail)
        @test !isfile(source[])
        state=TOML.parsefile(joinpath(directory,"first","state.toml"))
        attempt=joinpath(directory,"first",state["attempt"])
        reference_path=joinpath(attempt,"reference","calculation.jld2")
        before=read(reference_path)
        initial_session=read_calculation(reference_path).metadata.session
        @test isfile(joinpath(attempt,"candidate","points","1","complete.toml"))
        @test !isfile(joinpath(attempt,"candidate","points","2","complete.toml"))
        @test TOML.parsefile(joinpath(directory,"queued","state.toml"))["state"]=="pending"
        @test calls==[1.,2.,3.]
        @test events==[1]
        fail[]=false
        outcomes=resume_campaign(directory)
        @test all(row -> row.state===:complete,outcomes)
        @test calls==[1.,2.,3.,3.,4.,1.,5.]
        @test read(reference_path)==before
        @test read_calculation(reference_path).metadata.session==initial_session
        completed=first(outcomes).result
        @test completed.timings.execution.reference.reused
        @test completed.timings.execution.candidate.compute.reused_points == 1
        session_state=TOML.parsefile(joinpath(directory,"first","state.toml"))
        snapshot=TOML.parsefile(joinpath(directory,"sessions",session_state["session"]*".progress.toml"))
        observation=only(row for row in snapshot["benchmarks"] if row["id"]=="first")
        @test observation["candidate"]["reused"] == 1
        @test observation["candidate"]["completed"] == 3
        @test observation["candidate"]["timing_reused"]
        @test completed.metadata.session.id != initial_session.id
        points=read_calculation(joinpath(attempt,"candidate","calculation.jld2")).metadata.point_sessions
        @test points[1].id==initial_session.id
        @test points[2].id==completed.metadata.session.id
        @test numerical_input_sha256(only(completed.candidate_result.axes.problems))==numerical_input_sha256(model.problem)
        @test length(completed.candidate_result)==3
        # Archived implementation metadata is not part of legacy reuse identity.
        marker=joinpath(attempt,"reference","complete.toml")
        record=TOML.parsefile(marker);record["signature"]="legacy-source-tree-signature"
        open(io->TOML.print(io,record),marker,"w")
        count=length(calls)
        @test all(row -> row.state===:complete,resume_campaign(directory))
        @test length(calls)==count
        # A changed formulation cannot reuse the old scalar checkpoint.
        changed=BenchmarkCalculation(:candidate,model.problem,CheckpointBackend(99.))
        @test_throws r"inputs changed" Gauntlet._execute(changed;
            directory=joinpath(attempt,"candidate","points","1"),model)
        # Damage to retained data is still fatal, independently of source edits.
        open(io->write(io,"damage"),reference_path,"a")
        @test_throws r"integrity|checksum" resume_campaign(directory)
    end
end

@testitem "Gauntlet / native execution checkpoints retain nested cable geometry" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using LineCableModels, JLD2
    model=load_case(:cable_18kv_1000mm2_trefoil;
        variation=ExactOverrides(frequencies=[1.,37.]))
    formulation=Formulation()
    definition=benchmark_definition(model;id=:nested,source_file=@__FILE__,
        reference=formulation,formulations=formulation)
    mktempdir() do directory
        path=joinpath(directory,"declarations.jld2")
        JLD2.jldsave(path;definitions_bytes=Gauntlet._execution_bytes([definition]))
        restored=only(Gauntlet._read_execution(path,"definitions"))
        @test calculation_record(restored.reference)==calculation_record(definition.reference)
        @test calculation_record(restored.candidate)==calculation_record(definition.candidate)
        @test restored.model.port_order==model.port_order
    end
end

@testitem "Gauntlet / report and timing failures keep completed numerical work" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using LineCableModels, TOML
    const calls=Ref(0)
    const fail_at=Ref(typemax(Int))
    struct TimingBackend <: LineCableModels.Grammar.AbstractFormulation end
    function LineCableModels.compute(problem::LineParametersProblem,::TimingBackend;options=(;))
        calls[]+=1
        calls[]==fail_at[] && error("timing failed")
        z=fill(1.0+1im,2,2,length(problem.frequencies))
        LineParameters(PhaseDomain,z,z,copy(problem.frequencies))
    end
    model=load_case(:two_insulated_wires;variation=ExactOverrides(frequencies=[1.,37.]))
    reference=BenchmarkCalculation(:reference,model.problem,TimingBackend())
    candidate=BenchmarkCalculation(:candidate,model.problem,TimingBackend())
    definition(settings,tolerances=(;))=benchmark_definition(:stages,model.id,:fixture,
        @__FILE__,model,reference,candidate,settings,tolerances)
    mktempdir() do directory
        # A moment report cannot consume these ordinary phase-domain results.
        invalid=definition((statistics=(:mean,:std),))
        @test_throws Exception run_benchmark(invalid;directory)
        @test calls[]==2
        @test isfile(joinpath(directory,"reference","complete.toml"))
        @test isfile(joinpath(directory,"candidate","complete.toml"))
        fixed=run_benchmark(definition((bands=(:all,),));directory)
        @test calls[]==2
        @test fixed.timings.execution.reference.reused
        @test fixed.timings.execution.candidate.reused
    end
    mktempdir() do directory
        calls[]=0;fail_at[]=3
        timed=definition((bands=(:all,),),
            (performance=(minimum_speedup=2.,samples=1,seconds=1.),))
        @test_throws r"timing failed" run_benchmark(timed;directory)
        @test calls[]==3
        @test !isempty(read_benchmark(directory).analyses)
        fail_at[]=typemax(Int)
        retried=run_benchmark(timed;directory)
        @test calls[]==7 # Warmup and sample for each operand; numerical checkpoints are reused.
        @test retried.timings.execution.reference.reused
        @test retried.timings.execution.candidate.reused
    end
end
