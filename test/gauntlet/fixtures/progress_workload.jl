# Reusable warmed owned-MC worker for the manual overhead measurement.
using LineCableModels, Measurements, Logging, LinearAlgebra, TOML
include(joinpath(@__DIR__, "..", "..", "..", "gauntlet", "Gauntlet.jl"))
using .Gauntlet
BLAS.set_num_threads(1)
const model=load_case(:two_insulated_wires;variation=compose_variations(
    ExactOverrides(frequencies=collect(10.0 .^ range(0,5;length=48))),
    RelativeStandardUncertainty(2.;tags=(:geometry,:cable_layer))))
const problem=ParametricProblem(model.problem)
const inner=Gauntlet.uq_inner_formulation()
calculation(n)=BenchmarkCalculation(:reference,problem,MonteCarlo(inner;trials=n,seed=0x1234,
    return_samples=false,return_histograms=false))
function execute(mode,root,n)
    calc=calculation(n)
    tracker_ref=Ref{Any}(nothing)
    result=Ref{Any}()
    elapsed=@elapsed Logging.with_logger(Logging.NullLogger()) do
        Gauntlet._with_campaign_progress(root,[:owned_mc],"measure";progress=mode,io=devnull) do tracker
            tracker_ref[]=tracker
            receiver=LineCableModels.progress_receiver()
            if tracker!==nothing
                row=tracker.rows[1];row["resolved"]=true;row["case"]="48-frequency owned MC"
                merge!(row["reference"],Dict(string(k)=>v for (k,v) in pairs(Gauntlet._calculation_progress(calc))))
                merge!(row["candidate"],Dict("state"=>"skipped","estimate"=>0.0,"total"=>0))
                LineCableModels.report_progress(receiver,(kind=:benchmark,benchmark=:owned_mc,attempt="one",state=:running))
                LineCableModels.report_progress(receiver,(kind=:operand,benchmark=:owned_mc,attempt="one",role=:reference,
                    state=:running,backend="Owned",mode="Monte Carlo"))
            end
            seconds=@elapsed result[]=LineCableModels.with_progress_scope(benchmark=:owned_mc,attempt="one",role=:reference) do
                Gauntlet._compute_calculation(calc)
            end
            if receiver!==nothing
                LineCableModels.report_progress(receiver,(kind=:operand,benchmark=:owned_mc,attempt="one",role=:reference,state=:complete,seconds))
                LineCableModels.report_progress(receiver,(kind=:benchmark,benchmark=:owned_mc,attempt="one",state=:complete,seconds))
            end
        end
    end
    # Verify the declared number of accepted scans after the measured interval.
    @assert only(result[].trial_counts)==n
    digest=Gauntlet.semantic_sha256((Z=observe(only(result[]),Z),Y=observe(only(result[]),Y)))
    return elapsed,tracker_ref[]===nothing ? 0 : tracker_ref[].revision,digest
end
for mode in (:off,:auto)
    mktempdir() do root
        execute(mode,root,8)
    end
end
probe=mktempdir(root->first(execute(:off,root,40)))
const trials=max(40,round(Int,6.0/probe*40))
# Warm the exact allocation size and both collector paths before any reported run.
for mode in (:off,:auto)
    mktempdir(root->execute(mode,root,trials))
end
println("READY\t",trials,"\t",VERSION,"\t",Threads.nthreads(),"\t",BLAS.get_num_threads())
flush(stdout)
for line in eachline(stdin)
    mode,root=split(line,'\t';limit=2)
    GC.gc()
    elapsed,revisions,digest=execute(Symbol(mode),root,trials)
    println("RESULT\t",elapsed,"\t",revisions,"\t",digest)
    flush(stdout)
end
