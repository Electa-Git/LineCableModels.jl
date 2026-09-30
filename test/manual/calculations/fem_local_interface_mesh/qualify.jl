# Qualification only: native GetDP owns equations and all voltage integration.
# Usage: julia --project=test qualify.jl [pilot|spectrum]
# Completed bundles/meshes/solves are reused; no production method overrides.
using LineCableModels, Gmsh, SHA, TOML, Printf, Dates
const FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
const gmsh = Gmsh.gmsh
const ROOT = joinpath(pkgdir(LineCableModels), ".linecablemodels/fem/local-interface-mesh")
const REDUCTIONS = (reduce_bundle=false, kron_reduction=false, ideal_transposition=false)
const MESH = (domain_skin_depths=24.,
    pml_resolution=(interpolation_cells=72, coefficient_change=.12),
    mesh_size_factor=3., exterior_mesh_size_factor=8., volume_quadrature=12,
    conductor_geometry_tolerance=1e-3, conductor_skin_depth_elements=6.,
    conductor_mesh_growth=sqrt(1.25), conductor_skin_depths=5., conductor_thickness_elements=4)
const GETDP = FEM._getdp_selection(computation_options(LineCableModelsFEM, ComputationOptions())).path
say(xs...) = (println(Dates.now(), " ", xs...); flush(stdout))
digest(path) = bytes2hex(open(sha256, path))
function record(path, data)
    open(io -> TOML.print(io, data), path*".tmp", "w")
    mv(path*".tmp", path; force=true)
end

function fixture(layout, f, gamma_fraction)
    wire = build(CableDesign, "two_bare_wires",
        Stack(Group(:core, Region(:core_metal, Disk(.0425), Material(MaterialsLibrary(add_defaults=true), :copper)))))
    heights = layout == :mixed ? (-1., 1.) : layout == :air ? (1., 1.) : (-1., -1.)
    system = build(LineCableSystem, [wire, wire], [Pose2(0., heights[1]), Pose2(1., heights[2])];
        connections=[Dict(:core=>1), Dict(:core=>2)], system_id="interface-$layout", line_length=1.)
    problem = LineParametersProblem(system; frequencies=[f], temperature=20.,
        earth_props=homogeneous(rho=.1, eps_r=1., mu_r=1.))
    gamma = gamma_fraction*sqrt(im*2pi*f*4pi*1e-7*(10+im*2pi*f*8.8541878128e-12))
    form = Formulation(:LineCableModelsFEM; options=(; REDUCTIONS..., Γ=gamma))
    return problem, form
end

function prepare_bundle(dir, problem, form; live=false)
    isfile(joinpath(dir, "sources.toml")) && return
    mkpath(dir)
    if live
        source = joinpath(pkgdir(LineCableModels), ".linecablemodels/fem/onelab-two-bare-wires")
        for file in readlines(joinpath(source, ".onelab-export-files"))
            mkpath(dirname(joinpath(dir, file)))
            cp(joinpath(source, file), joinpath(dir, file); force=true)
        end
        cp(joinpath(source, ".onelab-export-files"), joinpath(dir, ".onelab-export-files"); force=true)
    else
        export_data(:onelab, problem, form; file_name=joinpath(dir, "study.pro"),
            mesh_options=MESH, overwrite=isfile(joinpath(dir, "study.pro")))
    end
    files = readlines(joinpath(dir, ".onelab-export-files"))
    record(joinpath(dir, "sources.toml"), Dict(file=>digest(joinpath(dir, file)) for file in files))
end

function coordinate_hash(entities)
    points = Set{NTuple{3,Float64}}()
    for (d,t) in entities
        _, xyz, _ = gmsh.model.mesh.get_nodes(d,t,true)
        union!(points, Tuple.(eachcol(reshape(xyz,3,:))))
    end
    bytes2hex(sha256(Vector{UInt8}(codeunits(repr(sort!(collect(points)))))))
end

function localize!(problem; media=:both, decay_footprint=false)
    interface = Set(gmsh.model.get_entities_for_physical_group(1,2002))
    cables = [c for i in eachindex(problem.system.designs)
        for c in gmsh.model.get_entities_for_physical_group(1,5000+i)]
    modified = 0
    changed = 0
    original = gmsh.model.mesh.field.list()
    for distance in original
        gmsh.model.mesh.field.get_type(distance) == "Distance" || continue
        curves = Int.(gmsh.model.mesh.field.get_numbers(distance,"CurvesList"))
        isempty(intersect(curves,interface)) && continue
        modified += 1
        media==:air && modified!=1 && continue
        media==:soil && modified!=2 && continue
        thresholds = [f for f in original if gmsh.model.mesh.field.get_type(f)=="Threshold" &&
            gmsh.model.mesh.field.get_number(f,"InField")==distance]
        if gmsh.model.mesh.field.get_number(only(thresholds),"SizeMin") >=
                gmsh.model.mesh.field.get_number(only(thresholds),"SizeMax")
            say("CONSTANT WAVE FIELD medium=",modified," unchanged")
            continue
        end
        if media==:wave
            k = Int(only(gmsh.parser.get_number("FrequencyIndex")))
            f = gmsh.parser.get_number("Frequencies")[k]
            gamma = complex(gmsh.parser.get_number("GammaReValues")[k],gmsh.parser.get_number("GammaImValues")[k])
            mu = modified==1 ? only(gmsh.parser.get_number("AirMu")) : gmsh.parser.get_number("EarthMu")[k]
            eps = modified==1 ? only(gmsh.parser.get_number("AirEpsilon")) : gmsh.parser.get_number("EarthEpsilon")[k]
            sigma = modified==1 ? 0. : gmsh.parser.get_number("EarthSigma")[k]
            limit = only(gmsh.parser.get_number("MeshSizeFactor"))/(8abs(sqrt(im*2pi*f*mu*(sigma+im*2pi*f*eps)-gamma^2)))
            remote = gmsh.model.mesh.field.get_number(only(thresholds),"SizeMax")
            say("WAVE CAP medium=",modified," remote=",remote," limit=",limit," localize=",remote<=limit)
            remote<=limit || continue
        end
        decay_width = decay_footprint ? gmsh.model.mesh.field.get_number(only(thresholds),"DistMin") : 0.
        gmsh.model.mesh.field.set_numbers(distance,"CurvesList",cables)
        sources = [distance]
        for (design,position) in zip(problem.system.designs,problem.system.positions)
            width = abs(position.y)+LineCableModels.outer_radius(design)+decay_width
            footprint = gmsh.model.mesh.field.add("MathEval")
            gmsh.model.mesh.field.set_string(footprint,"F",
                "Sqrt(y^2+Max(Abs(x-($(position.x)))-($width),0)^2)")
            push!(sources,footprint)
        end
        combined = gmsh.model.mesh.field.add("Min")
        gmsh.model.mesh.field.set_numbers(combined,"FieldsList",sources)
        for threshold in thresholds
            gmsh.model.mesh.field.set_number(threshold,"InField",combined)
        end
        changed += 1
    end
    modified == 2 || error("Expected two interface fields, got $modified")
    return changed
end

function anisotropize!(problem)
    # Preserve the existing normal target; change only the tangential metric.
    original = gmsh.model.mesh.field.list()
    background = maximum(original) # Native writer emits the final Min last.
    gmsh.model.mesh.field.get_type(background)=="Min" || error("Missing native background Min")
    interface = gmsh.model.get_entities_for_physical_group(1,2002)
    cables = [c for i in eachindex(problem.system.designs)
        for c in gmsh.model.get_entities_for_physical_group(1,5000+i)]
    normal_fields = [background]
    for restriction in original
        gmsh.model.mesh.field.get_type(restriction)=="Restrict" || continue
        threshold = round(Int,gmsh.model.mesh.field.get_number(restriction,"InField"))
        gmsh.model.mesh.field.get_type(threshold)=="Threshold" || continue
        distance = round(Int,gmsh.model.mesh.field.get_number(threshold,"InField"))
        gmsh.model.mesh.field.get_type(distance)=="Distance" || continue
        isempty(intersect(interface,gmsh.model.mesh.field.get_numbers(distance,"CurvesList"))) && continue
        copied_distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.set_numbers(copied_distance,"CurvesList",[interface;cables])
        gmsh.model.mesh.field.set_number(copied_distance,"Sampling",200)
        copied_threshold = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.set_number(copied_threshold,"InField",copied_distance)
        for key in ("SizeMin","SizeMax","DistMin","DistMax")
            gmsh.model.mesh.field.set_number(copied_threshold,key,gmsh.model.mesh.field.get_number(threshold,key))
        end
        copied_restriction = gmsh.model.mesh.field.add("Restrict")
        gmsh.model.mesh.field.set_number(copied_restriction,"InField",copied_threshold)
        gmsh.model.mesh.field.set_number(copied_restriction,"IncludeBoundary",1)
        gmsh.model.mesh.field.set_numbers(copied_restriction,"SurfacesList",
            gmsh.model.mesh.field.get_numbers(restriction,"SurfacesList"))
        push!(normal_fields,copied_restriction)
    end
    # localize! scans only the original fields so the preserved normal fields
    # must be added afterwards: temporarily save their source curve lists.
    copied_sources = [(f,gmsh.model.mesh.field.get_numbers(f,"CurvesList")) for f in
        setdiff(gmsh.model.mesh.field.list(),original) if gmsh.model.mesh.field.get_type(f)=="Distance"]
    for (f,_) in copied_sources
        gmsh.model.mesh.field.set_numbers(f,"CurvesList",cables)
    end
    localize!(problem)
    for (f,curves) in copied_sources
        gmsh.model.mesh.field.set_numbers(f,"CurvesList",curves)
    end
    normal = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.set_numbers(normal,"FieldsList",normal_fields)
    metric = gmsh.model.mesh.field.add("MathEvalAniso")
    for (key,field) in (("M11",background),("M22",normal),("M33",background))
        gmsh.model.mesh.field.set_string(metric,key,"1/(F$field^2)")
    end
    gmsh.model.mesh.field.set_as_background_mesh(metric)
    for tag in (1001,1002), surface in gmsh.model.get_entities_for_physical_group(2,tag)
        gmsh.model.mesh.set_algorithm(2,surface,7)
    end
    gmsh.option.set_number("Mesh.AnisoMax",8.)
end

function mesh!(dir, problem, candidate; media=:both)
    marker = joinpath(dir,"mesh.toml")
    if isfile(marker)
        data = TOML.parsefile(marker)
        digest(joinpath(dir,"study.msh")) == data["sha256"] || error("Mesh changed: $dir")
        return data
    end
    session = FEM._start_gmsh(4)
    try
        gmsh.parser.clear(); gmsh.onelab.clear()
        gmsh.parser.set_string("OnelabAction",["audit"])
        gmsh.open(joinpath(dir,"study.geo"))
        if candidate
            if media==:anisotropic
                anisotropize!(problem)
            else
                changed = localize!(problem;media=media==:decay ? :both : media,decay_footprint=media==:decay)
                baseline = joinpath(dirname(dir),"baseline")
                if changed==0 && isfile(joinpath(baseline,"mesh.toml"))
                    TOML.parsefile(joinpath(dir,"sources.toml"))==TOML.parsefile(joinpath(baseline,"sources.toml")) || error("No-op sources differ")
                    data = TOML.parsefile(joinpath(baseline,"mesh.toml"))
                    cp(joinpath(baseline,"study.msh"),joinpath(dir,"study.msh");force=true)
                    data["unchanged_size_fields"] = true
                    data["mesh_seconds"] = 0.
                    record(marker,data)
                    say("NO SIZE FIELD CHANGE; using exact baseline mesh in ",dir)
                    return data
                end
            end
        end
        measured = @timed gmsh.model.mesh.generate(2)
        groups(d,t) = [(d,c) for c in gmsh.model.get_entities_for_physical_group(d,t)]
        count2(tag) = sum(length(block) for (_,s) in groups(2,tag)
            for block in gmsh.model.mesh.get_elements(2,s)[2])
        data = Dict{String,Any}("nodes"=>length(first(gmsh.model.mesh.get_nodes())),
            "air_triangles"=>count2(1001),"soil_triangles"=>count2(1002),
            "pml_triangles"=>count2(1005),"mesh_seconds"=>measured.time,
            "compile_seconds"=>measured.compile_time,
            "pml_coordinates"=>coordinate_hash(groups(2,1005)),
            "path_coordinates"=>coordinate_hash(vcat(groups(1,7001),groups(1,7002))),
            "contour_coordinates"=>coordinate_hash(vcat(groups(1,4001),groups(1,4002))))
        gmsh.write(joinpath(dir,"study.msh"))
        data["sha256"] = digest(joinpath(dir,"study.msh"))
        record(marker,data)
        say("MESH ",basename(dir)," ",data)
        return data
    finally
        FEM._finish_gmsh(session)
    end
end

function solve!(dir; ordering=nothing, extra_args=String[])
    for (file,hash) in TOML.parsefile(joinpath(dir,"sources.toml"))
        digest(joinpath(dir,file)) == hash || error("Frozen source changed: $file")
    end
    output = joinpath(dir,"results/f0001-quasi-fw-b0000")
    marker = joinpath(dir,"solve.toml")
    if !isfile(joinpath(output,"completed.txt"))
        say("SOLVE BEGIN ",dir,"; two source columns; native output follows")
        cmd = `$GETDP $(joinpath(dir,"study.pro")) -msh $(joinpath(dir,"study.msh")) -solve LineCableModelsFEM -setnumber FrequencyIndex 1 -setnumber PlotFieldMaps 0 -v 4 -nt 1 -ksp_diagonal_scale -ksp_diagonal_scale_fix`
        ordering===nothing || (cmd = `$cmd -mat_mumps_icntl_7 $ordering`)
        cmd = `$cmd $extra_args`
        write(joinpath(dir,"command.txt"),string(cmd)*"\n")
        # Line buffering plus tee makes native progress visible in the one live log.
        run(pipeline(pipeline(`/usr/bin/time -f %e,%M -o $(joinpath(dir,"resources.csv")) stdbuf -oL -eL $cmd`;
            stderr=stdout), `tee $(joinpath(dir,"solver.log"))`))
        isfile(joinpath(output,"completed.txt")) || error("Missing native completion: $dir")
    end
    text = read(joinpath(dir,"solver.log"),String)
    costs = parse.(Float64,split(strip(read(joinpath(dir,"resources.csv"),String)),','))
    times = [parse.(Float64,split(strip(read(joinpath(output,"raw/jobs",
        @sprintf("getdp-f0001-b%04d-timing.tsv",basis)),String)))) for basis in 1:2]
    data = Dict("seconds"=>costs[1],"peak_rss_kib"=>round(Int,costs[2]),
        "dofs"=>maximum(parse(Int,m[1]) for m in eachmatch(r"System \d+/\d+: (\d+) Dofs",text)),
        "assembly_seconds"=>sum(t[4] for t in times),"solve_seconds"=>sum(t[5] for t in times))
    record(marker,data)
    say("SOLVE DONE ",dir," ",data)
    return data
end

function native_matrix(dir, quantity)
    a = zeros(ComplexF64,2,2)
    for row in split.(readlines(joinpath(dir,"results/f0001-quasi-fw-b0000/matrices/$quantity.tsv"))[3:end],'\t')
        a[parse(Int,row[1]),parse(Int,row[2])] = complex(parse(Float64,row[5]),parse(Float64,row[6]))
    end
    return a
end

function compare!(dir, problem, form; candidate="localized", tolerance=.01, record_prefix="")
    analytical = compute(problem,Formulation(earth_impedance=formula(:unified;options=(Γ=form.options.data.Γ,)),
        earth_admittance=formula(:unified;options=(Γ=form.options.data.Γ,));options=REDUCTIONS))
    worst = 0.; flips = 0
    suffix = candidate=="localized" ? "" : "-$candidate"
    open(joinpath(dir,"components$suffix.csv"),"w") do io
        println(io,"quantity,i,j,baseline,candidate,analytical,absolute_change,relative_change,sign_match")
        for (q,native,component) in (("R","Z",real),("X","Z",imag),("G","Y",real),("B","Y",imag))
            a = component.(native_matrix(joinpath(dir,"baseline"),native))
            b = component.(native_matrix(joinpath(dir,candidate),native))
            ref = component.(native == "Z" ? Z(analytical)[:,:,1] : Y(analytical)[:,:,1])
            for j in 1:2, i in 1:2
                relative = iszero(a[i,j]) ? (iszero(b[i,j]) ? 0. : Inf) : abs((b[i,j]-a[i,j])/a[i,j])
                match = sign(a[i,j]) == sign(b[i,j])
                worst = max(worst,relative); flips += !match
                println(io,join((q,i,j,a[i,j],b[i,j],ref[i,j],abs(b[i,j]-a[i,j]),relative,match),','))
            end
        end
    end
    passed = worst <= tolerance && flips == 0
    record(joinpath(dir,"$(record_prefix)comparison$suffix.toml"),Dict("passed"=>passed,"relative_tolerance"=>tolerance,"maximum_component_relative_change"=>worst,
        "new_component_sign_changes"=>flips))
    say("COMPARISON ",basename(dir)," / ",candidate," passed=",passed," largest component change=",worst," sign changes=",flips)
    return passed
end

function case!(layout, f, gamma_fraction; media=:both)
    label = "$layout-f$(f)-gamma$(gamma_fraction)"
    dir = mkpath(joinpath(ROOT,label))
    problem, form = fixture(layout,f,gamma_fraction)
    baseline = joinpath(dir,"baseline")
    candidate = media==:both ? "localized" : "localized-$media"
    localized = joinpath(dir,candidate)
    # The live mixed-case export is copied once. All later work uses that snapshot.
    prepare_bundle(baseline,problem,form;live=layout==:mixed && f==.1 && gamma_fraction==0)
    if !isfile(joinpath(localized,"sources.toml"))
        mkpath(localized)
        for file in [readlines(joinpath(baseline,".onelab-export-files"));".onelab-export-files";"sources.toml"]
            mkpath(dirname(joinpath(localized,file)))
            cp(joinpath(baseline,file),joinpath(localized,file);force=true)
        end
    end
    a = mesh!(baseline,problem,false)
    b = mesh!(localized,problem,true;media)
    for key in ("pml_coordinates","path_coordinates","contour_coordinates")
        a[key] == b[key] || error("Preservation failure $label: $key")
    end
    say("PRESERVATION PASS ",label,"; identical PML/path/contour nodes; nodes ",a["nodes"]," -> ",b["nodes"])
    solve!(baseline); solve!(localized)
    compare!(dir,problem,form;candidate) || error("Candidate failed qualification; completed data retained. Inspect $dir/components*.csv")
end

if abspath(PROGRAM_FILE) == @__FILE__
    stage = isempty(ARGS) ? "pilot" : only(ARGS)
    stage in ("pilot","spectrum","decay-spectrum","wave-spectrum") || error("Expected pilot, spectrum, decay-spectrum or wave-spectrum")
    media = stage=="decay-spectrum" ? :decay : stage=="wave-spectrum" ? :wave : :both
    mkpath(ROOT)
    say("BEGIN serial interface localization qualification stage=",stage)
    case!(:mixed,.1,0.;media)
    if stage != "pilot"
        for layout in (:air,:soil,:mixed), f in (.1,1e3,1e6)
            layout==:mixed && f==.1 && continue
            case!(layout,f,0.;media)
        end
        for layout in (:mixed,:air,:soil), f in (.1,1e6)
            case!(layout,f,.99;media)
        end
    end
    say("COMPLETE ",stage," qualification; feature implementation remains separate")
end
