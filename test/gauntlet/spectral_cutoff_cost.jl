# Compare warmed earth-matrix construction and reuse on the active package tree.
# To measure another checkpoint, run this same script with --project=<checkpoint>.
# julia --project=. test/gauntlet/spectral_cutoff_cost.jl <output.json>
using LineCableModels, LinearAlgebra, Statistics, JSON3, SHA
const E=LineCableModels.Engine
const ROOT=pkgdir(LineCableModels)
const SOURCE_FILES=("Engine.jl","integration.jl","spectralsampling.jl","earthbounds.jl",
    "compleximages.jl","earthreturn.jl","earthkernels.jl")
source_hashes()=Dict(file=>(isfile(joinpath(ROOT,"src/engine",file)) ?
    bytes2hex(sha256(read(joinpath(ROOT,"src/engine",file)))) : nothing) for file in SOURCE_FILES)
const INITIAL_HASHES=source_hashes()
function benchmark()
    cases=NamedTuple[
      (name="buried2",x=[0.,1.],h=[-1.,-1.],r=[.0425,.0425],f=f,rho=.1,er=1.,mr=1.,Γ=0.0im,reference=:deep) for f in (1.,1e4,1e6)]
    push!(cases,(name="buried3",x=[0.,1.,2.],h=[-1.,-1.,-1.],r=fill(.0425,3),f=1e6,rho=.1,er=1.,mr=1.,Γ=0.0im,reference=:deep))
    for (name,h) in (("overhead",[1.2,1.4]),("mixed",[1.2,-.9]))
      for reference in (:deep,:interface,:scalar)
        push!(cases,(;name,x=[0.,1.],h,r=[.01,.025],f=1e4,rho=10.,er=8.,mr=3.,Γ=1e-4+2e-4im,reference))
      end
    end
    rows=[]
    for c in cases
      geometry=E.EarthReturnGeometry(c.x,c.h,c.r)
      s=2pi*c.f*im
      sigma=[0.,1/c.rho]; epsilon=8.8541878128e-12.*[1.,c.er]; mu=4pi*1e-7.*[1.,c.mr]
      state=(jω=s,Γ=c.Γ,sigma,epsilon,mu,gamma_medium_squared=s.*mu.*(sigma.+s.*epsilon))
      reference=E.unified_earth!(E.EarthReturnWorkspace(geometry),state,E.computation_options(E.SpectralIntegral,(method=:quad,options=(rtol=1e-10,)));reference=c.reference)
      for method in (:quad,:trapz,:cim)
        workspace=E.EarthReturnWorkspace(geometry)
        controls=E.computation_options(E.SpectralIntegral,(;method,options=(rtol=method==:quad ? 1e-9 : 1e-6,)))
        print((;case=c.name,f=c.f,reference=c.reference,method));flush(stdout)
        try
          E.unified_earth!(workspace,state,controls;reference=c.reference)
          for mode in (:construction,:reuse)
            measurements=map(1:3) do _
              mode==:construction && empty!(workspace.scratch.numerical.cim.fits)
              timed=@timed E.unified_earth!(workspace,state,controls;reference=c.reference)
              (;seconds=timed.time,bytes=timed.bytes,evaluations=workspace.scratch.report.evaluations[],construction_samples=haskey(workspace.scratch.report,:samples) ? workspace.scratch.report.samples[] : nothing,cutoff=workspace.scratch.report.cutoff[],pencils=workspace.scratch.numerical.cim.statistics.pencils[],images=sum(fit->length(fit.images),workspace.scratch.numerical.cim.fits;init=0))
            end
            fastest=measurements[argmin(getproperty.(measurements,:seconds))]
            measurement=merge(fastest,(median_seconds=median(getproperty.(measurements,:seconds)),
                maximum_seconds=maximum(getproperty.(measurements,:seconds))))
            errors=map((:Ze,:Pe,:Ye)) do key
              a=getproperty(workspace,key);b=getproperty(reference,key)
              for (v,w) in zip(a,b), component in (real,imag)
                @assert isapprox(component(v),component(w);rtol=1e-5,atol=key==:Ye ? 1e-10 : 1e-14) (c,method,key,v,w)
              end
              norm(a-b)/max(norm(b),floatmin(Float64))
            end
            push!(rows,(;case=c.name,f=c.f,reference=c.reference,method,mode,measurement...,errors))
          end
          println(" passed");flush(stdout)
        catch err
          println(" FAILED ",sprint(showerror,err));flush(stdout)
          push!(rows,(;case=c.name,f=c.f,reference=c.reference,method,error=sprint(showerror,err)))
        end
      end
      open(ARGS[1],"w") do io
        JSON3.pretty(io,(;julia=string(VERSION),root=ROOT,blas_threads=BLAS.get_num_threads(),sources=INITIAL_HASHES,rows))
      end
    end
    @assert source_hashes()==INITIAL_HASHES "Sources changed during measurement"
    @assert all(row->!haskey(row,:error),rows) "A backend failed accuracy or convergence checks"
end
benchmark()
