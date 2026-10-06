# Explicit manual FEM/reference comparisons. No production refinement or acceptance rule.
# Invoke the command shown in the adjacent README with an output directory and suite selection.
using LineCableModels, Gmsh, JSON3, SHA, Printf, LinearAlgebra
include("closed_forms.jl")
const FEM=Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
const REDUCTIONS=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false)
const FREQUENCIES=[.001,.01,.1,1.,10.,100.,1e3,1e4,1e5,2e5,5e5,1e6,2e6,5e6,1e7,3e7,1e8,3e8]
const PLACEMENTS=[("overhead_pair",[1.,1.]),("mixed_pair",[1.,-1.]),("buried_pair",[-1.,-1.])]
const COPPER=Material(MaterialsLibrary(add_defaults=true),:copper)

function make_problem(heights,rho,eps_r,frequency;radius=.0425)
    designs=[build(CableDesign,"wire-$i",terminal(:core,Region(:metal,Disk(radius),COPPER))) for i in eachindex(heights)]
    positions=[Pose2(2(i-1),y) for (i,y) in pairs(heights)]
    system=build(LineCableSystem,designs,positions;
        connections=[Dict(:core=>i) for i in eachindex(heights)],line_length=1.,system_id="fem-validation")
    LineParametersProblem(system;temperature=20.,earth_props=homogeneous(;rho,eps_r,mu_r=1.),frequencies=[frequency])
end

function cases(group,settings)
    records=NamedTuple[]
    if group in ("main","all")
        frequencies=get(settings,"frequencies",FREQUENCIES)
        soils=get(settings,"earth",[(rho=r,eps_r=1.) for r in (.1,1.,100.,1000.)])
        radius=Float64(get(settings,"radius",.0425))
        for (placement,heights) in PLACEMENTS,soil in soils,f in frequencies
            rho=Float64(soil[:rho]);eps_r=Float64(soil[:eps_r])
            name="$(placement)_rho$(rho)_epsilon$(eps_r)_f$(Float64(f))"
            push!(records,(;name,problem=make_problem(heights,rho,eps_r,Float64(f);radius),gamma=0im,
                reference=:unified,radius,medium=nothing,reference_distance=nothing))
        end
    end
    if group in ("homogeneous","all")
        f=1e8;epsilon=8.8541878128e-12;mu=4pi*1e-7;k0=2pi*f*sqrt(mu*epsilon)
        for (placement,height) in (("overhead",1.),("buried",-1.)),c in (.5,.999,1.001,1.01,2.)
            push!(records,(name="homogeneous_air_$(placement)_rate$(c)",
                problem=make_problem([height],Inf,1.,f),gamma=im*c*k0,
                reference=:closed,radius=.0425,medium=(sigma=0.,epsilon,mu),
                reference_distance=height>0 ? height : nothing))
        end
        for radius in (.0425,.01)
            push!(records,(name="homogeneous_earth_buried_radius$(radius)",
                problem=make_problem([-1.],1000.,12.,f;radius),gamma=0im,
                reference=:closed,radius,medium=(sigma=.001,epsilon=12epsilon,mu),reference_distance=nothing))
        end
    end
    records
end

function call_native(command,path;timeout=2700.)
    start=time();process=open(path,"w") do io
        run(pipeline(ignorestatus(command);stdout=io,stderr=io);wait=false)
    end
    waitstatus=timedwait(()->process_exited(process),timeout;pollint=.1)
    if waitstatus===:timed_out
        kill(process);wait(process)
    end
    (success=success(process),timed_out=waitstatus===:timed_out,seconds=time()-start,exit_code=process.exitcode)
end

function matrix(path)
    lines=filter(l->!startswith(l,"#"),readlines(path));header=split(first(lines),'\t')
    rows=[split(l,'\t') for l in lines[2:end]]
    a=zeros(ComplexF64,maximum(parse(Int,r[1]) for r in rows),maximum(parse(Int,r[2]) for r in rows))
    ir=findfirst(==("real"),header);ii=findfirst(==("imag"),header)
    for r in rows;a[parse(Int,r[1]),parse(Int,r[2])]=complex(parse(Float64,r[ir]),parse(Float64,r[ii]));end
    a
end

function errors(actual,reference)
    scale=maximum(abs,diag(reference));delta=actual-reference
    percent(x,r)=iszero(r) ? nothing : 100x/r
    (scaled_percent=100maximum(abs,delta)/scale,
        entries=[(row=i,column=j,complex_error_percent=100abs(delta[i,j])/max(abs(reference[i,j]),eps(scale)),
            G_signed_percent=percent(real(delta[i,j]),real(reference[i,j])),
            B_signed_percent=percent(imag(delta[i,j]),imag(reference[i,j])),
            G_sign_agrees=sign(real(actual[i,j]))==sign(real(reference[i,j])))
            for i in axes(actual,1),j in axes(actual,2)],
        reciprocity=maximum(abs,actual-transpose(actual))/maximum(abs,diag(actual)))
end

# Engineering budgets are reference-scaled; tiny conductances remain visible but ungated.
function budget_metrics(actual_y, reference_y, actual_z, reference_z, case, native_mesh)
    scale=maximum(abs,diag(reference_y))
    delta=actual_y-reference_y
    conductance=[(row=i,column=j,relative_error=abs(real(delta[i,j]))/abs(real(reference_y[i,j])))
        for i in axes(reference_y,1),j in axes(reference_y,2)
        if abs(real(reference_y[i,j])) >= 1e-3*scale]
    sign_canaries=[(row=i,column=j,agrees=sign(real(actual_y[i,j]))==sign(real(reference_y[i,j])))
        for i in axes(reference_y,1),j in axes(reference_y,2)
        if abs(real(reference_y[i,j])) >= 1e-6*scale]
    reciprocity(a)=maximum(abs,a-transpose(a))/maximum(abs,diag(a))
    roots=native_mesh["transverse_roots"]
    receiver_sizes=[case.radius*hypot(roots[2m+1],roots[2m+2])
        for m in (case.problem.system.positions[i].y > 0 ? 0 : 1 for i in axes(reference_y,1))]
    air=first(case.problem.earth_props.layers)
    k0=2pi*only(case.problem.frequencies)*sqrt((4pi*1e-7*air.mu_r)*(8.8541878128e-12*air.eps_r))
    scope_reasons=String[]
    !iszero(case.gamma) && imag(case.gamma)<k0 && push!(scope_reasons,"prescribed wave faster than an air plane wave")
    any(>(.1),receiver_sizes) && push!(scope_reasons,"receiver transverse size exceeds mean-field scope")
    any(iszero,hypot(roots[2m+1],roots[2m+2]) for m in 0:1) && push!(scope_reasons,"exact transverse cutoff")
    reference_reciprocity=reciprocity(reference_y)
    excess=reciprocity(actual_y)-reference_reciprocity
    (e_Z=reference_z===nothing ? nothing : maximum(abs,actual_z-reference_z)/maximum(abs,diag(reference_z)),
        e_Y=maximum(abs,delta)/scale,
        worst_Y_entry=Tuple(argmax(abs.(delta))),
        worst_Z_entry=reference_z===nothing ? nothing : Tuple(argmax(abs.(actual_z-reference_z))),
        conductance_gate_pass=all(e.relative_error<=.05 for e in conductance),
        conductance_entries=conductance,conductance_failures=filter(e->e.relative_error>.05,conductance),
        significant_G_signs=sign_canaries,sign_canary_flag=any(!e.agrees for e in sign_canaries),
        ungated_G_sign_agreement=count(e.G_sign_agrees for e in errors(actual_y,reference_y).entries),
        ungated_G_sign_entries=length(reference_y),reference_reciprocity,
        fem_reciprocity=reciprocity(actual_y),reciprocity_excess_percentage_points=100excess,
        reciprocity_canary_flag=excess>.003,in_reference_scope=isempty(scope_reasons),scope_reasons,
        receiver_transverse_sizes=receiver_sizes)
end

function run_case(case,output,options,settings)
    directory=joinpath(output,case.name);mkpath(directory)
    form=Formulation(:LineCableModelsFEM;options=(;REDUCTIONS...,Γ=case.gamma))
    bundle=joinpath(directory,"bundle");entry=joinpath(bundle,"model.pro")
    isfile(entry) || export_data(:onelab,case.problem,form;file_name=entry,options=options)
    # Homogeneous-medium controls edit numerical material facts only, in the bundle.
    if case.reference===:closed && case.medium.sigma!=0
        data=joinpath(bundle,"model_data.pro");text=read(data,String)
        for (name,value) in (("AirSigma",case.medium.sigma),("AirEpsilon",case.medium.epsilon),("AirMu",case.medium.mu))
            text=replace(text,Regex("(?m)^$name = [^;]+;")=>"$name = $(repr(value));")
        end
        write(data,text)
    end
    source_root=joinpath(pkgdir(LineCableModels),"ext","LineCableModelsGmshExt")
    current_sources=join([read(joinpath(dir,file),String) for (dir,_,files) in walkdir(source_root) for file in sort(files) if endswith(file,".jl") || endswith(file,".pro") || endswith(file,".geo")],"\n")
    identity=bytes2hex(sha256(JSON3.write(settings) * read(@__FILE__,String) * read(joinpath(@__DIR__,"closed_forms.jl"),String) * current_sources * join([read(joinpath(dir,file),String) for (dir,_,files) in walkdir(bundle) for file in sort(files) if (endswith(file,".pro") || endswith(file,".geo"))],"\n")))
    recordpath=joinpath(directory,"result.json")
    if isfile(recordpath)
        saved=JSON3.read(read(recordpath,String));saved.source_identity==identity || error("Saved inputs differ: $(case.name)")
        return saved
    end
    impedance_reference=nothing
    reference=if case.reference===:unified
        equation=formula(:unified;options=(Γ=case.gamma,))
        result=compute(case.problem,Formulation(;earth_impedance=equation,earth_admittance=equation,options=REDUCTIONS))
        impedance_reference=Z(result)[:,:,1]
        Y(result)[:,:,1]
    else
        medium=case.medium
        # Read the exported conductor material facts, including temperature evaluation.
        data=read(joinpath(bundle,"model_data.pro"),String)
        scalar(name)=parse(Float64,match(Regex("(?m)^"*name*" = ([^;]+);"),data)[1])
        array(name)=parse(Float64,first(split(match(Regex("(?m)^"*name*"\\(\\) = \\{([^}]+)\\};"),data)[1],',')))
        impedance_reference=reshape([cylinder_impedance(only(case.problem.frequencies),case.radius,
            medium.sigma,medium.epsilon,medium.mu,case.gamma,
            array("Material_1_conductor_Sigma"),array("Material_1_conductor_Epsilon"),scalar("Material_1_conductor_Mu");
            reference_distance=case.reference_distance)],1,1)
        reshape([cylinder_admittance(only(case.problem.frequencies),case.radius,medium.sigma,
            medium.epsilon,medium.mu,case.gamma;reference_distance=case.reference_distance)],1,1)
    end
    open(joinpath(directory,"reference.csv"),"w") do io
        println(io,"row,column,real,imag")
        for i in axes(reference,1),j in axes(reference,2);v=reference[i,j];@printf(io,"%d,%d,%.17g,%.17g\n",i,j,real(v),imag(v));end
    end
    mesh=joinpath(directory,"mesh.msh");valuespath=joinpath(directory,"mesh-values.txt")
    mesh_run=call_native(`$(Gmsh.gmsh_jll.gmsh()) $(replace(entry,r"\.pro$"=>".geo")) -setstring MeshMetadataPath $valuespath -2 -o $mesh -v 3 -nt 1`,joinpath(directory,"mesh.log"))
    status=mesh_run.success ? "meshed" : "mesh failed"
    budget=nothing;observations=nothing;metrics=nothing;impedance_metrics=nothing;solve_run=nothing;preprocess_run=nothing;dofs=nothing
    if mesh_run.success
        controls=FEM.computation_options(FEM.LineCableModelsFEM,ComputationOptions())
        getdp=FEM._getdp_selection(controls).path
        preprocess_run=call_native(`$getdp $entry -msh $mesh -pre LineCableModelsFEM -name $(joinpath(directory,"solver")) -setnumber FrequencyIndex 1 -setnumber GetDPThreads 1 -v 4 -nt 1`,joinpath(directory,"preprocess.log");timeout=Float64(get(settings,"timeout_seconds",2700.)))
        prelog=read(joinpath(directory,"preprocess.log"),String)
        sizes=collect(eachmatch(r"System \d+/\d+: (\d+) Dofs",prelog))
        isempty(sizes) || (dofs=maximum(parse(Int,m[1]) for m in sizes))
        if !preprocess_run.success || dofs===nothing || dofs>Int(get(settings,"max_dofs",4_000_000))
            status=preprocess_run.timed_out ? "preprocess timeout" : !preprocess_run.success ? "preprocess failed" : dofs===nothing ? "missing DOF count" : "skipped: size"
            native_mesh=isfile(valuespath) ? FEM._read_native_mesh_values(valuespath) : nothing
            record=(;name=case.name,status,source_identity=identity,frequency_hz=only(case.problem.frequencies),
                gamma=(real=real(case.gamma),imag=imag(case.gamma)),radius_m=case.radius,
                reference=case.reference,options,mesh_run,preprocess_run,solve_run,dofs,observations,metrics,impedance_metrics,native_mesh,
                budget,material_override=case.reference===:closed && case.medium.sigma!=0 ? case.medium : nothing)
            open(recordpath,"w") do io;JSON3.pretty(io,record);end
            return record
        end
        solve_run=call_native(`$getdp $entry -msh $mesh -cal -name $(joinpath(directory,"solver")) -setnumber FrequencyIndex 1 -setnumber GetDPThreads 1 -v 4 -nt 1`,joinpath(directory,"native.log");timeout=Float64(get(settings,"timeout_seconds",2700.)))
        log=read(joinpath(directory,"native.log"),String);matches=collect(eachmatch(r"N:\s+(\d+)",log));isempty(matches) || (dofs=maximum(parse(Int,m[1]) for m in matches))
        status=solve_run.timed_out ? "timeout" : solve_run.success ? "solved" : "native failed"
        raw=joinpath(bundle,"results","f0001-helmholtz-b0000")
        if solve_run.success
            actual=matrix(joinpath(raw,"matrices","Y.tsv"));metrics=errors(actual,reference)
            impedance=matrix(joinpath(raw,"matrices","Z.tsv"))
            impedance_metrics=impedance_reference===nothing ? nothing : errors(impedance,impedance_reference)
            for (quantity,values) in (("Y",actual),("Z",impedance))
                open(joinpath(directory,"fem-"*quantity*".csv"),"w") do io
                    println(io,"row,column,real,imag")
                    for i in axes(values,1),j in axes(values,2)
                        @printf(io,"%d,%d,%.17g,%.17g\n",i,j,real(values[i,j]),imag(values[i,j]))
                    end
                end
            end
            budget=budget_metrics(actual,reference,impedance,impedance_reference,case,FEM._read_native_mesh_values(valuespath))
            observations=FEM._pml_observation(joinpath(raw,"raw/jobs/pml-f0001.tsv"),1,only(case.problem.frequencies))
            observations===nothing && error("Missing native observations: $(case.name)")
            observations.earth_sizing_ceiling_active && (status="solved, unqualified")
        end
    end
    native_mesh=isfile(valuespath) ? FEM._read_native_mesh_values(valuespath) : nothing
    record=(;name=case.name,status,source_identity=identity,frequency_hz=only(case.problem.frequencies),
        gamma=(real=real(case.gamma),imag=imag(case.gamma)),radius_m=case.radius,
        reference=case.reference,options,mesh_run,preprocess_run,solve_run,dofs,observations,metrics,impedance_metrics,native_mesh,
        budget,material_override=case.reference===:closed && case.medium.sigma!=0 ? case.medium : nothing)
    open(recordpath,"w") do io;JSON3.pretty(io,record);end
    record
end

function main(arguments)
    length(arguments) in (2,3) || error("Usage: run.jl OUTPUT main|homogeneous|all [SETTINGS.json]")
    output=abspath(arguments[1]);group=arguments[2];group in ("main","homogeneous","all") || error("Unknown comparison group")
    settings=length(arguments)==3 ? JSON3.read(read(arguments[3],String)) : Dict{String,Any}()
    options=get(settings,"options",Dict{String,Any}());options=NamedTuple(Symbol(k)=>(k=="overrides" ? NamedTuple(Symbol(n)=>x for (n,x) in pairs(v)) : v) for (k,v) in pairs(options))
    requests=cases(group,settings);mkpath(output)
    if get(settings,"list_only",false)
        open(joinpath(output,"cases.json"),"w") do io;JSON3.pretty(io,[(name=c.name,reference=c.reference,gamma=c.gamma) for c in requests]);end
        println(length(requests)," cases; no native execution")
        return
    end
    selected=get(settings,"case_indices",collect(eachindex(requests)))
    records=Any[]
    for i in selected
        case=requests[Int(i)]
        println("Starting ",i,"/",length(requests),": ",case.name);flush(stdout)
        record=try
            run_case(case,output,options,settings)
        catch exception
            directory=joinpath(output,case.name);mkpath(directory)
            message=sprint(showerror,exception,catch_backtrace())
            write(joinpath(directory,"tool-error.log"),message)
            failed=(name=case.name,status="tool failed",error=message)
            open(joinpath(directory,"result.json"),"w") do io;JSON3.pretty(io,failed);end
            failed
        end
        push!(records,record)
        println("Finished ",case.name,": ",record.status);flush(stdout)
    end
    open(joinpath(output,get(settings,"record_file","results.json")),"w") do io;JSON3.pretty(io,records);end
end
abspath(PROGRAM_FILE)==(@__FILE__) && main(ARGS)
