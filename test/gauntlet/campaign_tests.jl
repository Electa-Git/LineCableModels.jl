@testitem "Gauntlet / declarations remain authoritative through compute and recovery" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using LineCableModels, JLD2, SHA, TOML
    const calls=NamedTuple[]
    const fail_candidate=Ref(false)
    struct SpyBackend <: LineCableModels.Grammar.AbstractFormulation
        factor::Float64
    end
    function LineCableModels.compute(problem::LineParametersProblem, formulation::SpyBackend;options=(;))
        push!(calls,(;problem,formulation,options))
        fail_candidate[] && formulation.factor == 10 && error("deliberate interrupted operand")
        z=fill(complex(formulation.factor),2,2,length(problem.frequencies))
        y=fill(complex(0,formulation.factor),2,2,length(problem.frequencies))
        result=LineParameters(PhaseDomain,z,y,copy(problem.frequencies))
        haskey(options,:on_result) && options.on_result(problem,1,result)
        return result
    end
    model=load_case(:two_insulated_wires;variation=ExactOverrides(frequencies=[.01,3.,17.,400.],temperature=73.))
    problem=model.problem
    @test problem.frequencies == [.01,3.,17.,400.]
    events=Int[]
    callback=(problem,index,result)->push!(events,index)
    options=(tolerance=1e-9,custom_control=:unchanged,on_result=callback)
    reference=BenchmarkCalculation(:reference,problem,SpyBackend(1.);options)
    candidate=BenchmarkCalculation(:candidate,problem,SpyBackend(10.);options)
    definition=benchmark_definition(:authority,model.id,:fixture,@__FILE__,model,reference,candidate,
        (; quantities=(:Z,:Y,:G)),(;))
    direct=compute(problem,reference.formulation;options)
    empty!(calls);empty!(events)
    mktempdir() do parent
        directory=joinpath(parent,"campaign")
        fail_candidate[]=true
        failed=run_campaign(directory,[definition])
        @test only(failed).state === :failed
        state=TOML.parsefile(joinpath(directory,"authority","state.toml"))
        attempt=joinpath(directory,"authority",state["attempt"])
        @test isfile(joinpath(attempt,"reference","complete.toml"))
        @test !isfile(joinpath(attempt,"candidate","complete.toml"))
        @test only(campaign_status(directory)).state === :failed
        @test all(row -> row.problem === problem && row.options === options,calls)
        @test all(row -> row.problem.temperature == 73. && row.problem.frequencies == [.01,3.,17.,400.],calls)
        fail_candidate[]=false
        completed=resume_campaign(directory)
        @test only(completed).state === :complete
        value=only(completed).result
        @test value.reference.Z == direct.Z
        @test value.timings.execution.reference.reused
        @test !value.timings.execution.candidate.reused
        progress_state=TOML.parsefile(joinpath(directory,"authority","state.toml"))
        snapshot=TOML.parsefile(joinpath(directory,"sessions",progress_state["session"]*".progress.toml"))
        observation=only(snapshot["benchmarks"])
        @test observation["reference"]["saved_result"]
        @test !haskey(observation["reference"],"compute_seconds")
        @test observation["candidate"]["execution_seconds"] == value.timings.execution.candidate.seconds
        @test observation["candidate"]["compute_seconds"] == value.timings.execution.candidate.compute.seconds
        @test value.passes === nothing # Large cross-model differences are observations.
        @test all(==(9),only(row.error.relative for row in value.comparison if row.quantity === :Z && row.error.details.band === :all))
        @test only(campaign_status(directory)).state === :complete
        old_attempt=attempt
        fail_candidate[]=true
        replacement=only(run_campaign(directory,[definition]))
        @test replacement.state === :failed
        @test only(campaign_status(directory)).previous
        @test_throws r"latest draft" read_benchmark(joinpath(directory,"authority"))
        @test read_benchmark(joinpath(directory,"authority");previous=true).reference.result.Z == direct.Z
        @test isdir(old_attempt)
        fail_candidate[]=false
        recovered=only(resume_campaign(directory))
        @test recovered.state === :complete
        @test recovered.result.timings.execution.reference.reused
        @test !isdir(old_attempt)
        state=TOML.parsefile(joinpath(directory,"authority","state.toml"))
        attempt=joinpath(directory,"authority",state["current"])
        count=length(calls)
        before=read(joinpath(attempt,"reference","calculation.jld2"))
        @test only(resume_campaign(directory)).state === :complete
        @test length(calls) == count
        @test read(joinpath(attempt,"reference","calculation.jld2")) == before
        # New RMS bands retain the completed numerical operands.
        changed=benchmark_definition(:authority,model.id,:fixture,@__FILE__,model,reference,candidate,
            (; quantities=(:Z,:G),bands=((3.,17.),)),(;))
        run_benchmark(changed;directory=attempt)
        @test length(calls) == count
        @test read(joinpath(attempt,"reference","calculation.jld2")) == before
        reported=report(LineCableModels.ReportBuilder.BenchmarkTableDefinition(),read_benchmark(joinpath(directory,"authority")))
        @test length(unique(reported.table.maxima.snapshot))==2
        @test length(only(campaign_status(directory)).identity)==64
        bundle=lock_campaign(directory,joinpath(parent,"bundle"))
        moved=joinpath(parent,"elsewhere");mv(bundle.path,moved)
        rm(directory;recursive=true)
        retained=only(read_campaign(moved))
        @test retained.reference.result.Z == direct.Z
        @test length(retained.analyses) == 2
        @test_throws ArgumentError resume_campaign(moved)
        @test_throws ArgumentError lock_campaign(moved,moved)
        open(joinpath(moved,"authority","reference","calculation.jld2"),"a") do io
            write(io,"corruption")
        end
        @test_throws ArgumentError read_campaign(moved)
    end
end

@testitem "Gauntlet / numerical changes reject stale reuse and terminal changes reject comparison" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using LineCableModels
    model=load_case(:two_insulated_wires;variation=ExactOverrides(frequencies=[1.,3.]))
    problem=model.problem
    formulation=Formulation(options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
    first=BenchmarkCalculation(:a,problem,formulation)
    second=BenchmarkCalculation(:b,problem,formulation)
    definition=benchmark_definition(:changing,model.id,:fixture,@__FILE__,model,first,second,(;),(;))
    mktempdir() do directory
        outcome=run_benchmark(definition;directory)
        @test all(iszero,only(row.error.absolute for row in outcome.comparison if row.quantity === :Z && row.error.details.band === :all))
        changed=deepcopy(problem); changed.frequencies[2]=4.
        altered=benchmark_definition(:changing,model.id,:fixture,@__FILE__,model,
            BenchmarkCalculation(:a,changed,formulation),second,(;),(;))
        @test_throws ArgumentError run_benchmark(altered;directory)
        reordered=deepcopy(problem);reverse!(reordered.system.connection_order)
        mismatch=benchmark_definition(:different_ports,model.id,:fixture,@__FILE__,model,first,
            BenchmarkCalculation(:b,reordered,formulation),(;),(;))
        @test_throws ArgumentError run_benchmark(mismatch)
        original=read(joinpath(directory,"candidate","calculation.jld2"))
        open(joinpath(directory,"candidate","calculation.jld2"),"a") do io
            write(io,"damage")
        end
        @test_throws ArgumentError run_benchmark(definition;directory)
        write(joinpath(directory,"candidate","calculation.jld2"),original)
        @test all(iszero,only(row.error.absolute for row in run_benchmark(definition;directory).comparison if row.quantity === :Z && row.error.details.band === :all))
    end
end
