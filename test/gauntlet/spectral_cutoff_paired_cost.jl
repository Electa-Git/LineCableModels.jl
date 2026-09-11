# Interleave warmed Float64 engine calls from two source trees in one process.
# Both namespaces use the active environment's dependency versions and BLAS.
# julia --project=. test/gauntlet/spectral_cutoff_paired_cost.jl <baseline-tree> <output.json>
using LineCableModels, LinearAlgebra, Statistics, JSON3, SHA
length(ARGS)==2 || error("expected baseline source tree and output JSON path")
module SpectralCutoffBaseline end
const BASELINE=abspath(ARGS[1])
Base.include(SpectralCutoffBaseline,joinpath(BASELINE,"src/LineCableModels.jl"))
const ENGINES=(SpectralCutoffBaseline.LineCableModels.Engine,LineCableModels.Engine)
const ROOTS=(BASELINE,pkgdir(LineCableModels))
const FILES=("Engine.jl","integration.jl","spectralsampling.jl","earthbounds.jl",
    "compleximages.jl","earthreturn.jl","earthkernels.jl")
source_hashes(root)=Dict(file=>(isfile(joinpath(root,"src/engine",file)) ?
    bytes2hex(sha256(read(joinpath(root,"src/engine",file)))) : nothing) for file in FILES)
const HASHES=source_hashes.(ROOTS)

function cases()
    result=NamedTuple[(name="buried2",x=[0.,1.],h=[-1.,-1.],r=fill(.0425,2),
        f=f,rho=.1,er=1.,mr=1.,Γ=0.0im,reference=:deep) for f in (1.,1e4,1e6)]
    push!(result,(name="buried3",x=[0.,1.,2.],h=fill(-1.,3),r=fill(.0425,3),
        f=1e6,rho=.1,er=1.,mr=1.,Γ=0.0im,reference=:deep))
    for (name,h) in (("overhead",[1.2,1.4]),("mixed",[1.2,-.9])), reference in (:deep,:interface,:scalar)
        push!(result,(;name,x=[0.,1.],h,r=[.01,.025],f=1e4,rho=10.,er=8.,mr=3.,Γ=1e-4+2e-4im,reference))
    end
    return result
end
function prepare(E,c,method;rtol=method===:quad ? 1e-9 : 1e-6)
    geometry=E.EarthReturnGeometry(c.x,c.h,c.r)
    s=2pi*c.f*im
    sigma=[0.,inv(c.rho)];epsilon=8.8541878128e-12.*[1.,c.er];mu=4pi*1e-7.*[1.,c.mr]
    state=(jω=s,Γ=c.Γ,sigma,epsilon,mu,gamma_medium_squared=s.*mu.*(sigma.+s.*epsilon))
    controls=E.computation_options(E.SpectralIntegral,(;method,options=(;rtol)))
    return (;E,workspace=E.EarthReturnWorkspace(geometry),state,controls,reference=c.reference)
end
function solve(job)
    job.E.unified_earth!(job.workspace,job.state,job.controls;reference=job.reference)
end
function measured(job,mode)
    w=job.workspace
    mode===:construction && empty!(w.scratch.numerical.cim.fits)
    timed=@timed solve(job)
    report=w.scratch.report
    return (seconds=timed.time,bytes=timed.bytes,evaluations=report.evaluations[],
        construction_samples=haskey(report,:samples) ? report.samples[] : nothing,
        pencils=w.scratch.numerical.cim.statistics.pencils[],
        images=sum(fit->length(fit.images),w.scratch.numerical.cim.fits;init=0))
end
function check(actual,reference)
    ratios=Float64[]
    for key in (:Ze,:Pe,:Ye), (a,b) in zip(getproperty(actual,key),getproperty(reference,key)), component in (real,imag)
        floor=key===:Ye ? 1e-10 : 1e-14
        ratio=abs(component(a-b))/(floor+1e-5*abs(component(b)))
        @assert ratio<=1 (key,a,b,ratio)
        push!(ratios,ratio)
    end
    return maximum(ratios)
end
function benchmark()
    jobs=NamedTuple[]
    # Complete compilation and workspace preparation before any timed pairs.
    for c in cases()
        reference=prepare(last(ENGINES),c,:quad;rtol=1e-10)
        solve(reference)
        for method in (:quad,:trapz,:cim)
            pair=map(E->prepare(E,c,method),ENGINES)
            foreach(solve,pair)
            foreach(job->check(job.workspace,reference.workspace),pair)
            push!(jobs,(;case=c.name,frequency_Hz=c.f,reference=c.reference,method,pair,expected=reference.workspace))
            println((stage=:warmup,case=c.name,frequency=c.f,reference=c.reference,method));flush(stdout)
        end
    end
    rows=NamedTuple[]
    for job in jobs, mode in (:construction,:reuse)
        records=map(_->NamedTuple[],ENGINES)
        for repetition in 1:5
            for index in (isodd(repetition) ? (1,2) : (2,1))
                push!(records[index],measured(job.pair[index],mode))
            end
        end
        errors=map(item->check(item.workspace,job.expected),job.pair)
        summary=map(records) do values
            fastest=values[argmin(getproperty.(values,:seconds))]
            merge(fastest,(median_seconds=median(getproperty.(values,:seconds)),))
        end
        ratios=[records[1][i].seconds/records[2][i].seconds for i in 1:5]
        row=(case=job.case,frequency_Hz=job.frequency_Hz,reference=job.reference,method=job.method,mode,
            before=summary[1],after=summary[2],median_paired_speedup=median(ratios),
            paired_speedups=ratios,maximum_component_tolerance_ratios=errors)
        push!(rows,row)
        println((case=row.case,frequency=row.frequency_Hz,reference=row.reference,method=row.method,mode,speedup=row.median_paired_speedup));flush(stdout)
    end
    @assert source_hashes.(ROOTS)==HASHES "Sources changed during measurement"
    open(ARGS[2],"w") do io
        JSON3.pretty(io,(julia=string(VERSION),blas_threads=BLAS.get_num_threads(),roots=ROOTS,sources=HASHES,
            repetitions=5,all_warmups_precede_measurements=true,rows))
    end
end
benchmark()
