# Qualification mathematics only. No production includes, Gmsh or cable solves.
using LinearAlgebra, Printf, TOML, Dates

stretch(u,b)=1+(1-im)*b*u^3
stretch_derivative(u,b)=3(1-im)*b*u^2
coordinate(u,b)=u+(1-im)*b*u^4/4

function gauss_rule(n)
    k=collect(1:n-1)
    eigenproblem=eigen(SymTridiagonal(zeros(n),k./sqrt.(4 .* k.^2 .-1)))
    return (eigenproblem.values.+1)./2, eigenproblem.vectors[1,:].^2
end
const QUADRATURE=gauss_rule(16)

# -d/du[(1/s) dv/du] + Q^2 s v = 0; v(0)=1, v(1)=0.
# Q=q_normal*physical_thickness. The dimensional boundary response is T/L.
function discrete_response(nodes,Q,b;quadrature=QUADRATURE)
    n=length(nodes)-1
    diagonal=zeros(ComplexF64,n+1); off=zeros(ComplexF64,n)
    xi,weights=quadrature
    for i in 1:n
        h=nodes[i+1]-nodes[i]
        for (x,w) in zip(xi,weights)
            s=stretch(nodes[i]+h*x,b)
            stiffness=w/(h*s); mass=w*h*Q^2*s
            diagonal[i]+=stiffness+mass*(1-x)^2
            diagonal[i+1]+=stiffness+mass*x^2
            off[i]+=-stiffness+mass*x*(1-x)
        end
    end
    # Exact elimination of internal P1 degrees of freedom; no iterative solver.
    response=diagonal[n]
    for i in n-1:-1:1
        response=diagonal[i]-off[i]^2/response
    end
    response
end
function finite_response(Q,b)
    z=coordinate(1.,b)
    iszero(Q) && return inv(z)
    Q*(1+exp(-2Q*z))/(-expm1(-2Q*z))
end
function curvature(u,Q,b)
    s=stretch(u,b); ds=stretch_derivative(u,b)
    z=coordinate(u,b); zend=coordinate(1.,b)
    iszero(Q) && return -ds/zend
    den=-expm1(-2Q*zend)
    tail=exp(-2Q*(zend-z)); outgoing=exp(-Q*z)
    v=outgoing*(-expm1(-2Q*(zend-z)))/den
    dv=-Q*s*outgoing*(1+tail)/den
    Q^2*s^2*v+ds/s*dv
end

function physical_case(f,rho;domain_depths=24.,thickness_depths=24.,gamma_fraction=0.)
    mu=4pi*1e-7; eps0=8.8541878128e-12; omega=2pi*f
    delta=sqrt(rho/(pi*f*mu))
    D=max(5.,domain_depths*delta)
    L=max(5.,thickness_depths*delta)
    air=complex(0.,omega*sqrt(mu*eps0))
    earth=sqrt(complex(-omega^2*mu*eps0,omega*mu/rho))
    Gamma=gamma_fraction*earth
    A=-log(1e-10)/2
    return (;f,rho,delta,D,L,air,earth,Gamma,
        side_strength=4A/(imag(air)*L),bottom_strength=4A/(imag(earth)*L),
        clearance=D-1.085)
end
function modes(case,direction)
    media=direction==:side ? (case.air,case.earth) :
        direction==:top ? (case.air,) : (case.earth,)
    result=NamedTuple[]
    # Zero normal propagation isolates the transformed static compliance. It
    # has no finite outgoing relative-reflection certificate.
    push!(result,(name="static",Q=0.0im,weight=1.))
    for (m,gamma) in enumerate(media)
        tangential=sort!(unique([0.; imag(case.air).*[.5,.9,.99,.999,1.,1.001,1.01,1.1,2.];
            exp.(range(log(.01/case.clearance),log(12/case.clearance);length=30))]))
        for (j,k) in enumerate(tangential)
            q=sqrt(complex(k^2+gamma^2-case.Gamma^2))
            real(q)<0 && (q=-q)
            weight=exp(-2real(q)*case.clearance)
            push!(result,(name="medium$m-mode$j",Q=q*case.L,weight))
        end
    end
    result
end

function density_grid(case,direction)
    b=direction==:bottom ? case.bottom_strength : case.side_strength
    spectrum=modes(case,direction)
    # Dense quadrature grid for an a-priori interpolation estimate. Resolves the
    # physical stretch transition b*u^3≈1; this is not a sequence of FEM meshes.
    transition=min(1.,b^(-1/3))
    u=sort!(unique([0.;exp.(range(log(transition*1e-5),0.;length=2501));
        collect(range(0.,1.;length=2501))]))
    density=zeros(length(u))
    for mode in spectrum
        scale=abs(finite_response(mode.Q,b))
        for i in eachindex(u)
            w=mode.weight*abs2(curvature(u[i],mode.Q,b))/(abs(stretch(u[i],b))*scale)
            density[i]=max(density[i],cbrt(w))
        end
    end
    integral=zeros(length(u))
    for i in 2:length(u)
        integral[i]=integral[i-1]+(u[i]-u[i-1])*(density[i]+density[i-1])/2
    end
    return (;u,density,integral,b,spectrum)
end
function equidistribute(grid,n)
    nodes=[0.]
    for target in range(0.,grid.integral[end];length=n+1)[2:end-1]
        j=searchsortedfirst(grid.integral,target)
        fraction=(target-grid.integral[j-1])/(grid.integral[j]-grid.integral[j-1])
        push!(nodes,grid.u[j-1]+fraction*(grid.u[j]-grid.u[j-1]))
    end
    push!(nodes,1.)
end
legacy_nodes(n)=expm1.(((192/191)*log(1536)).*(0:n)./n)./expm1((192/191)*log(1536))

function assess(nodes,grid,case)
    worst=0.; truncation=0.; weighted=0.; min_growth=Inf
    for mode in grid.spectrum
        exact=finite_response(mode.Q,grid.b)
        value=discrete_response(nodes,mode.Q,grid.b)
        error=abs(value-exact)/abs(exact)
        worst=max(worst,error);weighted=max(weighted,mode.weight*error)
        if !iszero(mode.Q)
            truncation=max(truncation,mode.weight*abs((exact-mode.Q)/(exact+mode.Q)))
        end
        min_growth=min(min_growth,minimum(real(mode.Q*stretch(u,grid.b)) for u in nodes))
    end
    n=length(nodes)-1
    estimate=grid.integral[end]^3/(12n^2)
    (;n,weighted_discretization=weighted,maximum_discretization=worst,
        finite_wall_reflection=truncation,interpolation_estimate=estimate,
        minimum_modal_attenuation=min_growth)
end

function run_modal(root)
    mkpath(root)
    println(Dates.now()," MODAL START: P1 boundary response, finite-wall and outgoing references; no cable solves");flush(stdout)
    # Uniform unstretched validation: finite interval exact response and P1
    # convergence. A static field must be represented exactly on any grid.
    for Q in (0im,1+im,0.2+2im)
        exact=finite_response(Q,0.)
        errors=[abs(discrete_response(collect(range(0.,1.;length=n+1)),Q,0.)-exact) for n in (32,64,128)]
        if iszero(Q)
            @assert maximum(errors)<1e-11
        else
            @assert errors[2]<.26errors[1] && errors[3]<.26errors[2]
        end
    end
    open(joinpath(root,"modal.csv"),"w") do io
        println(io,"frequency_hz,rho,domain_depths,thickness_depths,gamma_fraction,direction,distribution,count,weighted_discretization,maximum_discretization,finite_wall_reflection,interpolation_estimate,minimum_modal_attenuation")
        for (f,rho) in ((.1,.1),(21.544346900318832,.1),(1e6,.1),(1e6,100.)),
                (domain_depths,thickness_depths) in ((24.,24.),(2.,24.),(2.,2.)), gamma_fraction in (0.,.99)
            case=physical_case(f,rho;domain_depths,thickness_depths,gamma_fraction)
            for direction in (:side,:top,:bottom)
                grid=density_grid(case,direction)
                for (label,n) in (("legacy",144),("legacy",192),("interpolation",48),("interpolation",72),("interpolation",96))
                    nodes=label=="legacy" ? legacy_nodes(n) : equidistribute(grid,n)
                    metrics=assess(nodes,grid,case)
                    println(io,join((f,rho,domain_depths,thickness_depths,gamma_fraction,direction,label,values(metrics)...),','))
                    if gamma_fraction==0 && label=="interpolation" && n==72
                        directory=joinpath(root,"nodes",string(f,"-",rho,"-",domain_depths,"-",thickness_depths))
                        mkpath(directory)
                        open(joinpath(directory,string(direction,".csv")),"w") do out
                            println(out,"u,distance_m")
                            for u in nodes; println(out,u,',',u*case.L); end
                        end
                    end
                end
            end
            flush(io)
            println(Dates.now()," MODAL DONE f=",f," rho=",rho," D/delta=",domain_depths,
                " L/delta=",thickness_depths," Gamma/gamma_earth=",gamma_fraction);flush(stdout)
        end
    end
    println(Dates.now()," MODAL COMPLETE; no production sources changed and no coupled cable solves used");flush(stdout)
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    run_modal(isempty(ARGS) ? joinpath(pwd(),".linecablemodels/fem/pml-physical-mesh") : abspath(only(ARGS)))
end
