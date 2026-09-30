# Included only by run_two_bare_wires_fem.jl's explicit comparison/review mode.
# Serial, fixed prescriptions; completed points are data, never acceptance gates.
using TOML, SHA, Dates

const REVIEW_ROOT = joinpath(pkgdir(LineCableModels), ".linecablemodels/fem/native-feature-review")
const RETAINED_NATIVE = joinpath(pkgdir(LineCableModels), ".linecablemodels/fem/native-performance-20260930")
const REVIEW_CHOICES = (
    (name="triangle12", label="Consolidated; triangles, physical/PML 12/12",
        options=(pml_element_family=:triangle,physical_volume_quadrature=12,pml_quadrature=9)),
    (name="triangle3", label="Consolidated; triangles, physical/PML 3/12",
        options=(pml_element_family=:triangle,physical_volume_quadrature=3,pml_quadrature=9)),
    (name="quad9", label="Consolidated; quadrangles, physical/PML 3/9",
        options=(pml_element_family=:quadrangle,physical_volume_quadrature=3,pml_quadrature=9)))
review_say(xs...) = (println(Dates.now(), " ", xs...); flush(stdout))
review_hash(path) = bytes2hex(open(sha256,path))
function review_record(path,data)
    mkpath(dirname(path))
    open(io->TOML.print(io,data),path*".tmp","w")
    mv(path*".tmp",path;force=true)
end
function review_native_matrix(dir,q)
    a=zeros(ComplexF64,2,2)
    for row in split.(readlines(joinpath(dir,"results/f0001-quasi-fw-b0000/matrices/$q.tsv"))[3:end],'\t')
        a[parse(Int,row[1]),parse(Int,row[2])]=complex(parse(Float64,row[5]),parse(Float64,row[6]))
    end
    a
end
function review_matrices(data)
    (reshape(complex.(data["Z_real"],data["Z_imag"]),2,2),
     reshape(complex.(data["Y_real"],data["Y_imag"]),2,2))
end
function review_point(z,y;kwargs...)
    merge(Dict("Z_real"=>vec(real.(z)),"Z_imag"=>vec(imag.(z)),
        "Y_real"=>vec(real.(y)),"Y_imag"=>vec(imag.(y))),Dict(string(k)=>v for (k,v) in kwargs))
end

function review_baseline!(dir,problem,form,mesh,mesh_controls)
    f=only(problem.frequencies)
    retained=joinpath(RETAINED_NATIVE,"mixed-f$f-gamma0.0","baseline")
    if isfile(joinpath(retained,"solve.toml"))
        review_say("REUSE frozen baseline ",retained)
        dir=retained
    else
        entry=joinpath(dir,"study.pro")
        if !isfile(entry)
            export_data(:onelab,problem,form;file_name=entry,mesh_options=mesh_controls)
            # The comparison baseline is the frozen revision, never another
            # formulation branch in the production solver.
            frozen=joinpath(RETAINED_NATIVE,"snapshot/ext/LineCableModelsGmshExt/getdp")
            for name in readdir(frozen)
                endswith(name,".pro") || continue
                cp(joinpath(frozen,name),joinpath(dir,"formulations",name);force=true)
            end
        end
        completed=joinpath(dir,"results/f0001-quasi-fw-b0000/completed.txt")
        if !isfile(completed)
            fem=Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
            getdp=fem._getdp_selection(computation_options(LineCableModelsFEM,ComputationOptions())).path
            cmd=`$getdp $entry -msh $mesh -solve LineCableModelsFEM -setnumber PlotFieldMaps 0 -v 4 -nt 1`
            write(joinpath(dir,"command.txt"),string(cmd)*"\n")
            review_say("BASELINE BEGIN ",f," Hz")
            run(pipeline(pipeline(`/usr/bin/time -f %e,%M -o $(joinpath(dir,"resources.csv")) stdbuf -oL -eL $cmd`;
                stderr=stdout),`tee $(joinpath(dir,"solver.log"))`))
            isfile(completed) || error("Native baseline did not publish matrices: $dir")
        end
    end
    costs=parse.(Float64,split(strip(read(joinpath(dir,"resources.csv"),String)),','))
    times=[parse.(Float64,split(read(joinpath(dir,"results/f0001-quasi-fw-b0000/raw/jobs",
        @sprintf("getdp-f0001-b%04d-timing.tsv",b)),String))) for b in 1:2]
    review_point(review_native_matrix(dir,"Z"),review_native_matrix(dir,"Y");
        directory=dir,frequency_hz=f,source="frozen pre-consolidation revision",
        assembly_seconds=sum(t[4] for t in times),solve_seconds=sum(t[5] for t in times),
        native_seconds=costs[1],managed_seconds=NaN,compile_seconds=NaN)
end

function run_native_review(fs,form,ordinary_options)
    form.options.data.Γ == 0 && temperature == 20. && line_length == 1. &&
        soil_relative_permittivity == 1. && soil_relative_permeability == 1. ||
        error("This saved-reference comparison uses Γ=0, 20 °C, 1 m and earth εr=μr=1; use the ordinary runner for other physical settings.")
    fem=Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    # Match the retained reference, explicitly. The ordinary sweep keeps 5/10.
    controls=merge(ordinary_options,(mesh_size_factor=3.,exterior_mesh_size_factor=8.,
        volume_quadrature=12,frequency_workers=1,solver_threads=1,
        getdp_verbosity=4,timing=true,verbosity=(default=1,)))
    mesh_controls=(; (k=>controls[k] for k in keys(controls) if k in fem.FEM_EXPORT_MESH_OPTIONS)...)
    sources=Dict(relpath(joinpath(d,n),pkgdir(LineCableModels))=>review_hash(joinpath(d,n))
        for (d,_,ns) in walkdir(joinpath(pkgdir(LineCableModels),"ext/LineCableModelsGmshExt")) for n in ns)
    signature=bytes2hex(sha256(repr((sources,controls,form.options.data,fs,REVIEW_CHOICES))))
    manifest=joinpath(REVIEW_ROOT,"settings.toml")
    if isfile(manifest)
        TOML.parsefile(manifest)["signature"]==signature || error(
            "Review settings/sources changed. Preserve these results and select a new REVIEW_ROOT.")
    else
        review_record(manifest,Dict("signature"=>signature,"sources"=>sources,
            "frequencies_hz"=>fs,"controls"=>repr(controls),
            "geometry"=>"rho=0.1 ohm m, r=0.0425 m, x=(0,1) m, y=(-1,1) m, Gamma=0",
            "choices"=>Dict(c.name=>repr(c.options) for c in REVIEW_CHOICES)))
    end
    for (index,f) in enumerate(fs)
        directory=joinpath(REVIEW_ROOT,@sprintf("f%02d",index)); mkpath(directory)
        problem=build_two_bare_wires_problem(.0425,.1,[f];vert=-1.)
        retained=joinpath(RETAINED_NATIVE,"mixed-f$f-gamma0.0","baseline","study.msh")
        triangular_mesh=isfile(retained) ? retained : nothing
        for choice in REVIEW_CHOICES
            target=joinpath(directory,choice.name*".toml")
            if isfile(target)
                data=TOML.parsefile(target)
                choice.name=="triangle12" && (triangular_mesh=joinpath(data["directory"],"mesh/model.msh"))
                review_say("REUSE feature ",choice.name," ",f," Hz")
                continue
            end
            supplied=choice.options.pml_element_family===:triangle ? triangular_mesh : nothing
            selected=merge(controls,choice.options,(mesh_path=supplied,
                mesh_policy=supplied===nothing ? :remesh : :reuse))
            review_say("MANAGED BEGIN ",choice.name," ",f," Hz; ",choice.options)
            measured=@timed compute(problem,form;options=selected)
            value=measured.value
            detail=details(value).data
            timing=detail.timing
            data=review_point(only(eachslice(observe(value,Z);dims=3)),only(eachslice(observe(value,Y);dims=3));
                frequency_hz=f,directory=detail.fem.run.run_directory,source="feature code",
                assembly_seconds=timing.assembly_seconds,solve_seconds=timing.solve_seconds,
                native_seconds=timing.worker_wall_seconds,managed_seconds=measured.time,
                compile_seconds=measured.compile_time)
            review_record(target,data)
            choice.name=="triangle12" && (triangular_mesh=joinpath(data["directory"],"mesh/model.msh"))
            review_say("MANAGED DONE ",choice.name," ",f," Hz; native=",data["native_seconds"],
                " s; full call=",measured.time," s; compilation=",measured.compile_time," s; ",data["directory"])
        end
        baseline=joinpath(directory,"baseline.toml")
        if !isfile(baseline)
            data=review_baseline!(joinpath(directory,"baseline"),problem,form,triangular_mesh,
                merge(mesh_controls,(pml_element_family=:triangle,physical_volume_quadrature=nothing)))
            review_record(baseline,data)
        end
        analytical=joinpath(directory,"analytical.toml")
        if !isfile(analytical)
            reference=compute(problem,Formulation(earth_impedance=formula(:unified;options=(Γ=0.,)),
                earth_admittance=formula(:unified;options=(Γ=0.,));options=(reduce_bundle=false,
                    kron_reduction=false,ideal_transposition=false)))
            review_record(analytical,review_point(observe(reference,Z)[:,:,1],observe(reference,Y)[:,:,1];frequency_hz=f))
        end
        review_say("POINT COMPLETE ",index,"/",length(fs)," f=",f," Hz; all numerical differences retained")
    end
    # One independently runnable native bundle of the selected combined choice.
    entry=joinpath(REVIEW_ROOT,"onelab-quad9","study.pro")
    if !isfile(entry)
        export_data(:onelab,build_two_bare_wires_problem(.0425,.1,fs;vert=-1.),form;
            file_name=entry,mesh_options=merge(mesh_controls,last(REVIEW_CHOICES).options))
    end
    review_say("COMPLETE managed comparisons and detached export ",entry)
end
