isdefined(@__MODULE__,:PILOTS) || include("qualify.jl")
include("../fem_conductor_mesh/fixtures.jl")

function fixture_problem(name,f)
    if name=="three"
        src=joinpath(pkgdir(LineCableModels),
            ".linecablemodels/fem/conductor-mesh-qualification/bare-wire-mvp/three-air_1/quasi_fw/problem.json")
        original=import_data(:json,LineParametersProblem;file_name=src)
        return LineParametersProblem(original.system;frequencies=[f],
            temperature=original.temperature,earth_props=original.earth_props)
    end
    design=getproperty(ConductorMeshFixtures,Symbol(name=="screen" ? "screened_cable" : name=="tube" ? "tubular_cable" : "sector_cable"))()
    system=build(LineCableSystem,design,Pose2(0.,name=="sector" ? 1. : -1.);
        connections=Dict(t=>i for (i,t) in enumerate(design.terminal_order)),line_length=1.)
    LineParametersProblem(system;frequencies=[f],temperature=20.,
        earth_props=homogeneous(rho=100.,eps_r=1.,mu_r=1.))
end

function fixture_bundle!(dir,problem,form)
    isfile(joinpath(dir,"sources.toml")) && return
    export_data(:onelab,problem,form;file_name=joinpath(dir,"study.pro"),
        mesh_options=B.MESH,overwrite=isfile(joinpath(dir,"study.pro")))
    for file in readdir(joinpath(dir,"formulations"))
        frozen=joinpath(ROOT,"snapshot/ext/LineCableModelsGmshExt/getdp",file)
        isfile(frozen) && cp(frozen,joinpath(dir,"formulations",file);force=true)
    end
    record(joinpath(dir,"sources.toml"),Dict(f=>digest(joinpath(dir,f)) for f in readlines(joinpath(dir,".onelab-export-files"))))
end

function matrix_n(dir,q,n)
    a=zeros(ComplexF64,n,n)
    for row in split.(readlines(joinpath(dir,"results/f0001-quasi-fw-b0000/matrices/$q.tsv"))[3:end],'\t')
        a[parse(Int,row[1]),parse(Int,row[2])]=complex(parse(Float64,row[5]),parse(Float64,row[6]))
    end
    a
end

function fixtures!(variant; cases=[(name,f) for name in ("three","screen","tube","sector") for f in (.1,1e6)], stop_on_failure=false)
    passed=true
    completed=0
    form=Formulation(:LineCableModelsFEM;options=(;B.REDUCTIONS...,Γ=0.))
    for (name,f) in cases
        problem=fixture_problem(name,f)
        n=length(problem.system.terminal_order)
        root=joinpath(ROOT,"fixtures","$name-f$f")
        base=joinpath(root,"baseline"); dir=joinpath(root,variant)
        fixture_bundle!(base,problem,form); B.mesh!(base,problem,false)
        if !isfile(joinpath(dir,"sources.toml"))
            files=readlines(joinpath(base,".onelab-export-files"))
            for file in [files;".onelab-export-files";"study.msh";"mesh.toml"]
                dest=joinpath(dir,file); mkpath(dirname(dest)); cp(joinpath(base,file),dest;force=true)
            end
            variant in ("harmonic","combined") && harmonic!(joinpath(dir,"formulations/quasi-full.pro"))
            variant in ("physical3","combined") && integration!(joinpath(dir,"formulations/integration.pro"),variant)
            record(joinpath(dir,"sources.toml"),Dict(file=>digest(joinpath(dir,file)) for file in files))
        end
        for target in (base,dir)
            B.solve!(target)
            # solve! executes all native columns; extend its two-column timing summary.
            data=TOML.parsefile(joinpath(target,"solve.toml"))
            times=[parse.(Float64,split(read(joinpath(target,"results/f0001-quasi-fw-b0000/raw/jobs",
                @sprintf("getdp-f0001-b%04d-timing.tsv",b)),String))) for b in 1:n]
            data["assembly_seconds"]=sum(t[4] for t in times)
            data["solve_seconds"]=sum(t[5] for t in times)
            data["source_columns"]=n
            record(joinpath(target,"solve.toml"),data)
        end
        worst=0.; flips=0
        open(joinpath(root,"components-$variant.csv"),"w") do io
            println(io,"quantity,i,j,baseline,candidate,relative_change,sign_match")
            for (q,field,part) in (("R","Z",real),("X","Z",imag),("G","Y",real),("B","Y",imag))
                a=part.(matrix_n(base,field,n)); b=part.(matrix_n(dir,field,n))
                for j in 1:n,i in 1:n
                    relative=iszero(a[i,j]) ? (iszero(b[i,j]) ? 0. : Inf) : abs((b[i,j]-a[i,j])/a[i,j])
                    match=sign(a[i,j])==sign(b[i,j]); worst=max(worst,relative); flips+=!match
                    println(io,join((q,i,j,a[i,j],b[i,j],relative,match),','))
                end
            end
        end
        ok=worst<=.02 && flips==0; passed &= ok
        record(joinpath(root,"comparison-$variant.toml"),Dict("passed"=>ok,"maximum_component_relative_change"=>worst,
            "new_component_sign_changes"=>flips,"source_columns"=>n))
        completed+=1
        say("FIXTURE ",name," f=",f," variant=",variant," passed=",ok," worst=",worst," sign changes=",flips)
        stop_on_failure && !ok && break
    end
    record(joinpath(ROOT,"fixtures-$variant.toml"),Dict("passed"=>passed,"cases"=>completed))
    say("COMPLETE fixture endpoints ",variant," passed=",passed)
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && fixtures!(only(ARGS))
