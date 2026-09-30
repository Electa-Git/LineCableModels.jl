# Two independent diagnostics. Current feature sources are never modified.
# Run: julia --project=test diagonals_and_soil.jl
isdefined(@__MODULE__, :B) || include("qualify.jl")
const DIAGNOSTICS = joinpath(ROOT, "diagonals-and-soil")

function opposite_diagonals!(baseline, target)
    marker = joinpath(target, "mesh.toml")
    if isfile(marker)
        data = TOML.parsefile(marker)
        digest(joinpath(target, "study.msh")) == data["sha256"] || error("Changed candidate mesh")
        return data
    end
    copy_bundle!(baseline, target)
    session = FEM._start_gmsh(4)
    try
        gmsh.parser.clear(); gmsh.onelab.clear()
        gmsh.open(joinpath(baseline, "study.msh"))
        gmsh.parser.parse(joinpath(target, "study_data.pro"))
        before = inspect_mesh()
        all_nodes = coordinate_hash(gmsh.model.get_entities())
        nt, coordinates, _ = gmsh.model.mesh.get_nodes(-1,-1,false,false)
        xyz = Dict(Int(t)=>Tuple(p) for (t,p) in zip(nt,eachcol(reshape(coordinates,3,:))))
        area(a,b,c) = ((xyz[b][1]-xyz[a][1])*(xyz[c][2]-xyz[a][2]) -
            (xyz[b][2]-xyz[a][2])*(xyz[c][1]-xyz[a][1]))/2
        changed = 0
        minimum_area = Inf
        max_area_change = 0.
        for surface in reduce(vcat,corner_groups())
            types, blocks, connectivity = gmsh.model.mesh.get_elements(2,surface)
            types == [2] || error("Expected linear triangles")
            tags = Int.(only(blocks))
            triangles = Int.(reshape(only(connectivity),3,:))
            pairs = Dict{Tuple{Int,Int},Vector{Int}}()
            for (i,tri) in enumerate(eachcol(triangles))
                # Native transfinite interpolation leaves tiny coordinate
                # roundoff on nominally horizontal/vertical edges. Identify
                # the diagonal relative to this cell, not global coordinates.
                edges = ((tri[1],tri[2]),(tri[2],tri[3]),(tri[3],tri[1]))
                dx = [abs(xyz[a][1]-xyz[b][1]) for (a,b) in edges]
                dy = [abs(xyz[a][2]-xyz[b][2]) for (a,b) in edges]
                diagonal = [edges[k] for k in 1:3
                    if min(dx[k]/maximum(dx),dy[k]/maximum(dy))>.99]
                a,b = only(diagonal)
                push!(get!(pairs,minmax(a,b),Int[]),i)
            end
            for ((a,b),indices) in pairs
                length(indices)==2 || error("Unpaired rectangle diagonal")
                i,j = indices
                c = only(setdiff(triangles[:,i],[a,b]))
                d = only(setdiff(triangles[:,j],[a,b]))
                old_area = abs(area(triangles[:,i]...))+abs(area(triangles[:,j]...))
                for (index,vertices) in ((i,[c,d,a]),(j,[d,c,b]))
                    if area(vertices...)<0
                        vertices[2],vertices[3] = vertices[3],vertices[2]
                    end
                    minimum_area = min(minimum_area,area(vertices...))
                    triangles[:,index] = vertices
                end
                new_area = area(triangles[:,i]...)+area(triangles[:,j]...)
                max_area_change = max(max_area_change,abs(new_area-old_area)/old_area)
                changed += 1
            end
            gmsh.model.mesh.remove_elements(2,surface,tags)
            gmsh.model.mesh.add_elements_by_type(surface,2,tags,vec(triangles))
        end
        after = inspect_mesh()
        for key in keys(before)
            before[key] == after[key] || error("Fixed-node diagnostic changed $key")
        end
        coordinate_hash(gmsh.model.get_entities()) == all_nodes || error("Node coordinates changed")
        minimum_area>0 && max_area_change<1e-12 || error("Invalid rectangle replacement")
        gmsh.option.set_number("Mesh.Renumber",0)
        gmsh.option.set_number("Mesh.SaveAll",1)
        gmsh.write(joinpath(target,"study.msh"))
        merge!(after,Dict("sha256"=>digest(joinpath(target,"study.msh")),
            "flipped_rectangles"=>changed,"minimum_triangle_area"=>minimum_area,
            "maximum_cell_area_relative_change"=>max_area_change,
            "all_node_coordinates"=>all_nodes,"physical_and_strip_elements_identical"=>true))
        record(marker,after)
        say("OPPOSITE DIAGONALS ",after)
        return after
    finally
        FEM._finish_gmsh(session)
    end
end

function soil_bundles!()
    problem,form = B.fixture(:mixed,.1,0.)
    problem = LineParametersProblem(problem.system;
        frequencies=collect(10. .^ range(-1.,6.;length=10)),temperature=20.,
        earth_props=homogeneous(rho=.1,eps_r=1.,mu_r=1.))
    root = mkpath(joinpath(DIAGNOSTICS,"soil-localization"))
    localized,baseline = joinpath(root,"localized"),joinpath(root,"baseline")
    # This is a fresh detached public export, also usable for interactive inspection.
    B.prepare_bundle(localized,problem,form)
    if !isfile(joinpath(baseline,"sources.toml"))
        copy_bundle!(localized,baseline)
        for file in readdir(joinpath(baseline,"geometry");join=true)
            endswith(file,".geo") || continue
            content = read(file,String)
            length(findall("If(MeshWaveEarth < MeshRemoteEarth)",content))==1 || error("Unexpected exported soil field")
            write(file,replace(content,"If(MeshWaveEarth < MeshRemoteEarth)"=>
                "If(0) // qualification control: retain the full-width soil field"))
        end
        record(joinpath(baseline,"sources.toml"),Dict(file=>digest(joinpath(baseline,file))
            for file in readlines(joinpath(baseline,".onelab-export-files"))))
    end
    a = B.mesh!(baseline,problem,false)
    b = B.mesh!(localized,problem,false)
    for key in ("pml_coordinates","path_coordinates","contour_coordinates")
        a[key]==b[key] || error("Soil diagnostic changed $key")
    end
    say("SOIL LOCALIZATION MESHES; PML/path/conductor coordinates identical; air ",
        a["air_triangles"]," -> ",b["air_triangles"],"; soil ",
        a["soil_triangles"]," -> ",b["soil_triangles"])
    B.solve!(baseline); B.solve!(localized)
    # Only the first frequency is solved in this focused comparison.
    first_problem,first_form = B.fixture(:mixed,.1,0.)
    passed = B.compare!(root,first_problem,first_form;tolerance=.02)
    return root,passed
end

function run_diagonals_and_soil()
    mkpath(DIAGNOSTICS)
    say("BEGIN serial fixed-node diagonal and independent soil-localization diagnostics")
    case = (:air,.1,0.)
    baseline = freeze!(case)
    target = joinpath(dirname(baseline),"opposite-diagonals")
    opposite_diagonals!(baseline,target)
    B.solve!(target)
    problem,form = B.fixture(case...)
    diagonal_passed = B.compare!(dirname(baseline),problem,form;
        candidate="opposite-diagonals",tolerance=.02)
    root,soil_passed = soil_bundles!()
    record(joinpath(DIAGNOSTICS,"complete.toml"),Dict(
        "diagonal_preservation_passed"=>diagonal_passed,"soil_preservation_passed"=>soil_passed,
        "soil_comparison_frequency_hz"=>.1,"fresh_detached_study"=>joinpath(root,"localized/study.pro"),
        "production_changed"=>false,"existing_user_bundle_changed"=>false))
    say("COMPLETE fixed-node diagonal and soil-localization diagnostics")
end

if abspath(PROGRAM_FILE)==@__FILE__
    open(joinpath(ROOT,"live.log"),"a") do log
        redirect_stdio(stdout=log,stderr=log) do
            run_diagonals_and_soil()
        end
    end
end
