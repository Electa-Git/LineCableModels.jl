@testitem "Gauntlet / incremental work skips cases and reuses only matching operands" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using LineCableModels, TOML, JLD2
    const calls=Float64[]
    struct SelectiveBackend <: LineCableModels.Grammar.AbstractFormulation
        factor::Float64
    end
    function LineCableModels.compute(problem::LineParametersProblem, f::SelectiveBackend;options=(;))
        push!(calls,f.factor)
        z=fill(complex(f.factor),2,2,length(problem.frequencies))
        return LineParameters(PhaseDomain,z,z,copy(problem.frequencies))
    end
    model=load_case(:two_insulated_wires;variation=ExactOverrides(frequencies=[1.,37.]))
    definition(id,a,b;bands=(:all,),tolerances=(;))=benchmark_definition(id,model.id,:fixture,
        @__FILE__,model,BenchmarkCalculation(:reference,model.problem,SelectiveBackend(a)),
        BenchmarkCalculation(:candidate,model.problem,SelectiveBackend(b)),(;bands),tolerances)
    mktempdir() do parent
        directory=joinpath(parent,"campaign")
        first=definition(:first,1.,2.)
        second=definition(:second,3.,4.)
        dry=run_campaign(directory,[first,second];dry_run=true,progress=:off)
        @test !ispath(directory)
        @test isempty(calls)
        @test all(row->row.action==:new && row.reference==:compute && row.candidate==:compute,dry)
        @test_throws r"unknown benchmark" run_campaign(directory,[first];benchmark=[:missing])
        @test !ispath(directory)
        @test all(row->row.state==:complete,run_campaign(directory,[first,second];progress=:off))
        @test calls==[1.,2.,3.,4.]
        first_state=joinpath(directory,"first","state.toml")
        before=read(first_state)
        old=TOML.parsefile(first_state)
        old_attempt=joinpath(directory,"first",old["attempt"])
        checkpoint=joinpath(old_attempt,"reference","calculation.jld2")
        bytes=read(checkpoint)
        session=read_calculation(checkpoint).metadata.session
        skipped=run_campaign(directory,[first,second];progress=:off)
        @test all(row->row.skipped && row.result===nothing,skipped)
        @test calls==[1.,2.,3.,4.]
        @test read(first_state)==before
        @test all(row->row.skipped,resume_campaign(directory;progress=:off))
        @test read(first_state)==before
        @test_throws r"inputs changed" run_campaign(directory,[definition(:first,1.,99.)];
            resume=true,dry_run=true,progress=:off)
        changed=definition(:first,1.,5.)
        plan=only(run_campaign(directory,[changed];dry_run=true,progress=:off))
        @test plan.reference==:reuse && plan.candidate==:compute
        @test read(first_state)==before
        updated=only(run_campaign(directory,[changed];progress=:off))
        @test calls==[1.,2.,3.,4.,5.]
        @test updated.result.timings.execution.reference.reused
        @test updated.result.timings.execution.reference.session==session
        @test isdir(old_attempt)
        @test read(checkpoint)==bytes
        new_state=TOML.parsefile(first_state)
        new_attempt=joinpath(directory,"first",new_state["attempt"])
        @test read(joinpath(new_attempt,"reference","calculation.jld2"))==bytes
        # A report-only change does not execute either operand.
        reported=only(run_campaign(directory,[definition(:first,1.,5.;bands=(:all,:dc))];progress=:off))
        @test calls==[1.,2.,3.,4.,5.]
        @test all(e->e.reused,reported.result.timings.execution)
        @test only(run_campaign(directory,[second];force=true,progress=:off)).state==:complete
        @test calls==[1.,2.,3.,4.,5.,3.,4.]
        @test_throws ArgumentError run_campaign(directory,[second];force=true,resume=true)
        @test only(resume_campaign(directory;benchmark=[:second],progress=:off)).id==:second
        @test_throws r"unknown benchmark" resume_campaign(directory;benchmark=[:missing],progress=:off)
        # Report-only changes must also retain controlled timing observations,
        # not quietly reexecute the numerical workload through the timing path.
        tolerances=(performance=(minimum_speedup=2.,samples=1,seconds=0.01),)
        measured=only(run_campaign(directory,
            [definition(:measured,6.,7.;tolerances)];progress=:off))
        measured_count=length(calls)
        measured_state=TOML.parsefile(joinpath(directory,"measured","state.toml"))
        measured_path=joinpath(directory,"measured",measured_state["current"],"performance.jld2")
        measured_session=JLD2.load(measured_path,"session")
        report_only=only(run_campaign(directory,
            [definition(:measured,6.,7.;bands=(:all,:dc),tolerances)];progress=:off))
        @test length(calls)==measured_count
        @test report_only.result.performance==measured.result.performance
        @test all(execution->execution.reused,report_only.result.timings.execution)
        report_state=TOML.parsefile(joinpath(directory,"measured","state.toml"))
        report_path=joinpath(directory,"measured",report_state["current"],"performance.jld2")
        @test JLD2.load(report_path,"session")==measured_session
        @test isfile(measured_path)
        changed_timing=(performance=merge(tolerances.performance,(samples=2,)),)
        run_campaign(directory,[definition(:measured,6.,7.;tolerances=changed_timing)];progress=:off)
        @test length(calls)>measured_count
    end
end

@testitem "Gauntlet / complete resume never restores executable declarations" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using LineCableModels, JLD2, TOML, SHA
    model=load_case(:two_insulated_wires;variation=ExactOverrides(frequencies=[1.,37.]))
    formulation=Formulation()
    definition=benchmark_definition(model;id=:skip,source_file=@__FILE__,
        reference=formulation,formulations=formulation)
    mktempdir() do directory
        run_campaign(directory,[definition];progress=:off)
        state_path=joinpath(directory,"skip","state.toml")
        before=read(state_path)
        state=TOML.parsefile(state_path)
        attempt=joinpath(directory,"skip",state["attempt"])
        path=joinpath(attempt,"declarations.jld2")
        # Readable integrity-checked evidence with deliberately unexecutable bytes.
        # A skip cannot depend on native deserialization or current builder methods.
        document=JLD2.load(path)
        document["definitions_bytes"]=UInt8[0x00]
        JLD2.jldsave(path;(Symbol(key)=>value for (key,value) in document)...)
        write(path*".sha256",bytes2hex(open(sha256,path)))
        outcome=only(resume_campaign(directory;progress=:plain))
        @test outcome.skipped
        @test read(state_path)==before
        snapshot=last(sort(filter(p->endswith(p,".progress.toml"),
            readdir(joinpath(directory,"sessions");join=true));by=mtime))
        row=only(TOML.parsefile(snapshot)["benchmarks"])
        @test row["state"]=="skipped"
        @test all(row[role]["saved_result"] for role in ("reference","candidate"))
        @test all(row[role]["completed"]==1 for role in ("reference","candidate"))
        @test all(!haskey(row[role],"compute_seconds") for role in ("reference","candidate"))
        @test !isdir(joinpath(attempt,"sources"))
        open(io->write(io,"damage"),joinpath(attempt,"reference","calculation.jld2"),"a")
        @test_throws r"checksum" resume_campaign(directory;progress=:off)
    end
end

@testitem "Gauntlet / changed sweeps retain matching formulation points" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using LineCableModels
    const calls=Float64[]
    struct SweepReuseBackend <: LineCableModels.Grammar.AbstractFormulation
        factor::Float64
    end
    Base.NamedTuple(formulation::SweepReuseBackend)=(factor=formulation.factor,)
    function LineCableModels.compute(problem::LineParametersProblem,f::SweepReuseBackend;options=(;))
        push!(calls,f.factor)
        z=fill(complex(f.factor),2,2,length(problem.frequencies))
        LineParameters(PhaseDomain,z,z,copy(problem.frequencies))
    end
    model=load_case(:two_insulated_wires;variation=ExactOverrides(frequencies=[1.,37.]))
    definition(factors)=benchmark_definition(model;id=:sweep,source_file=@__FILE__,
        reference=SweepReuseBackend(1.),
        formulations=Gridspace{SweepReuseBackend}(SweepReuseBackend,(Grid(factors),)))
    mktempdir() do directory
        run_campaign(directory,[definition((2.,3.,4.))];progress=:off,on_error=:fail)
        updated=only(run_campaign(directory,[definition((2.,5.,4.))];progress=:off,on_error=:fail))
        @test calls==[1.,2.,3.,4.,5.]
        @test updated.result.timings.execution.candidate.compute.reused_points==2
        @test length(updated.result.candidate_result)==3
        reordered=only(run_campaign(directory,[definition((4.,2.,6.,5.))];progress=:off,on_error=:fail))
        @test calls==[1.,2.,3.,4.,5.,6.]
        @test reordered.result.timings.execution.candidate.compute.reused_points==3
        @test length(reordered.result.candidate_result)==4
    end
end
