# Qualification prototype: approximate the interpolation density by six native
# geometric strips. Each strip preserves its endpoints and cell count exactly.
include("modal.jl")

function coefficient_guard_grid(grid;relative_variation=.15,interpolation_cells=72)
    density=max.(grid.density.*(interpolation_cells/grid.integral[end]),
        abs.(stretch_derivative.(grid.u,grid.b)./stretch.(grid.u,grid.b))./relative_variation)
    integral=zeros(length(density))
    for i in 2:length(density)
        integral[i]=integral[i-1]+(grid.u[i]-grid.u[i-1])*(density[i]+density[i-1])/2
    end
    (;grid...,density,integral)
end

function strip_grid(grid,n;strip_count=6)
    target=equidistribute(grid,n)
    boundaries=round.(Int,range(0,n;length=strip_count+1)).+1
    strips=NamedTuple[]; nodes=[0.]
    for (i,j) in zip(boundaries[1:end-1],boundaries[2:end])
        count=j-i
        ratio=((target[j]-target[j-1])/(target[i+1]-target[i]))^(1/(count-1))
        start=target[i];stop=target[j]
        push!(strips,(;start,stop,count,ratio))
        g=count*log(ratio)
        localnodes=abs(g)<1e-10 ? collect(1:count)./count : expm1.(g.*(1:count)./count)./expm1(g)
        append!(nodes,start.+(stop-start).*localnodes)
    end
    (;strips,nodes)
end

function guarded_strips(grid;relative_variation=.15)
    prescribed=coefficient_guard_grid(grid;relative_variation)
    # Put a native patch boundary at the physical stretch transition; a single
    # geometric progression cannot represent the spacing minimum across it.
    ordinary=equidistribute(prescribed,6)
    boundaries=sort!(unique([ordinary;min(1.,grid.b^(-1/3))]))
    function cumulative(u)
        i=clamp(searchsortedlast(prescribed.u,u),1,length(prescribed.u)-1)
        t=(u-prescribed.u[i])/(prescribed.u[i+1]-prescribed.u[i])
        (1-t)*prescribed.integral[i]+t*prescribed.integral[i+1]
    end
    function position(target)
        i=clamp(searchsortedlast(prescribed.integral,target),1,length(prescribed.integral)-1)
        t=(target-prescribed.integral[i])/(prescribed.integral[i+1]-prescribed.integral[i])
        (1-t)*prescribed.u[i]+t*prescribed.u[i+1]
    end
    strips=NamedTuple[];nodes=[0.]
    for (start,stop) in zip(boundaries[1:end-1],boundaries[2:end])
        lo,hi=cumulative(start),cumulative(stop)
        count=max(3,ceil(Int,hi-lo))
        first_size=position(lo+(hi-lo)/count)-start
        last_size=stop-position(hi-(hi-lo)/count)
        ratio=(last_size/first_size)^(1/(count-1))
        push!(strips,(;start,stop,count,ratio))
        g=count*log(ratio)
        localnodes=abs(g)<1e-10 ? collect(1:count)./count : expm1.(g.*(1:count)./count)./expm1(g)
        append!(nodes,start.+(stop-start).*localnodes)
    end
    (;strips,nodes)
end

function qualify_strips(root)
    println(Dates.now()," STRIPS START: six native geometric strips; Gamma=0 production scope");flush(stdout)
    quadrature32=gauss_rule(32)
    open(joinpath(root,"strips.csv"),"w") do io
        println(io,"frequency_hz,rho,domain_depths,thickness_depths,direction,count,weighted_discretization,maximum_discretization,finite_wall_reflection,interpolation_estimate,minimum_modal_attenuation,quadrature_change")
        for (f,rho) in ((.1,.1),(21.544346900318832,.1),(1e6,.1),(1e6,100.)),
                (domain_depths,thickness_depths) in ((24.,24.),(2.,24.),(2.,2.))
            case=physical_case(f,rho;domain_depths,thickness_depths)
            designs=Dict{String,Any}()
            for direction in (:side,:top,:bottom)
                grid=density_grid(case,direction)
                fit=strip_grid(grid,72)
                metric=assess(fit.nodes,grid,case)
                change=maximum(grid.spectrum) do mode
                    fine=discrete_response(fit.nodes,mode.Q,grid.b;quadrature=quadrature32)
                    coarse=discrete_response(fit.nodes,mode.Q,grid.b)
                    abs(fine-coarse)/abs(finite_response(mode.Q,grid.b))
                end
                println(io,join((f,rho,domain_depths,thickness_depths,direction,values(metric)...,change),','))
                designs[string(direction)]=[Dict(string(k)=>v for (k,v) in pairs(strip)) for strip in fit.strips]
            end
            path=joinpath(root,"strips",string(f,"-",rho,"-",domain_depths,"-",thickness_depths,".toml"))
            mkpath(dirname(path));open(out->TOML.print(out,designs),path,"w")
            println(Dates.now()," STRIPS DONE ",basename(path));flush(stdout)
        end
    end
    println(Dates.now()," STRIPS COMPLETE; zero coupled solves");flush(stdout)
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    qualify_strips(isempty(ARGS) ? joinpath(pwd(),".linecablemodels/fem/pml-physical-mesh") : abspath(only(ARGS)))
end
