# Manual, serial PML-corner qualification. No production method replacement.
module ReferenceTools
include("../fem_local_interface_mesh/qualify.jl")
end
using LineCableModels, Gmsh, SHA, TOML, Dates, Printf, Statistics
const B = ReferenceTools
const gmsh = Gmsh.gmsh
const FEM = B.FEM
const ROOT = joinpath(pkgdir(LineCableModels), ".linecablemodels/fem/pml-corner-mesh")
const PILOTS = ((:air,.1,0.), (:mixed,.1,.99), (:mixed,1e6,.99))
# Fixed before any native candidate solve. No search or automatic refinement.
const CANDIDATES = ((name="metric-1", scale=1.), (name="metric-1.35", scale=1.35))
say(xs...) = B.say(xs...)
record(args...) = B.record(args...)
digest(path) = B.digest(path)
label(c) = "$(c[1])-f$(c[2])-gamma$(c[3])"
value(name) = only(gmsh.parser.get_number(name))
entities(d, group) = [(d,Int(t)) for t in unique(gmsh.model.get_entities_for_physical_group(d,group))]

function freeze!(case)
    target = joinpath(ROOT,label(case),"baseline")
    isfile(joinpath(target,"sources.toml")) && return target
    old = joinpath(B.ROOT,label(case),case[3]==0 ? "localized" : "localized-terminal-shift")
    files = readlines(joinpath(old,".onelab-export-files"))
    for file in [files;".onelab-export-files";"study.msh";"mesh.toml"]
        mkpath(dirname(joinpath(target,file)))
        cp(joinpath(old,file),joinpath(target,file);force=true)
    end
    # Current corrected PDE on both reference and candidate; all physical data
    # and existing mesh coordinates are copied, never regenerated.
    for file in readdir(joinpath(target,"formulations"))
        source = joinpath(pkgdir(LineCableModels),"ext/LineCableModelsGmshExt/getdp",file)
        isfile(source) && cp(source,joinpath(target,"formulations",file);force=true)
    end
    record(joinpath(target,"sources.toml"), Dict(f=>digest(joinpath(target,f)) for f in files))
    record(joinpath(target,"reference.toml"), Dict("original_bundle"=>old,
        "original_mesh_sha256"=>digest(joinpath(old,"study.msh")),
        "frozen_mesh_sha256"=>digest(joinpath(target,"study.msh")),
        "equations"=>"current corrected quasi-fw; identical for both meshes"))
    say("FROZEN ",label(case)," from ",old)
    return target
end

function copy_bundle!(source, target)
    isfile(joinpath(target,"sources.toml")) && return
    for file in [readlines(joinpath(source,".onelab-export-files"));".onelab-export-files";"sources.toml"]
        mkpath(dirname(joinpath(target,file)))
        cp(joinpath(source,file),joinpath(target,file);force=true)
    end
end

function corner_groups()
    D, xc = value("DomainHalfwidths"), value("Xcenter")
    groups = [Int[] for _ in 1:4]
    for (_,s) in entities(2,1005)
        x0,y0,_,x1,y1,_ = gmsh.model.get_bounding_box(2,s)
        x,y = (x0+x1)/2, (y0+y1)/2
        abs(x-xc)>D && abs(y)>D || continue
        push!(groups[(x>xc ? 1 : 0)+(y>0 ? 2 : 0)+1],s)
    end
    return groups
end

function mesh_hash(surfaces)
    # Include coordinates AND connectivity. Node numbering is retained here.
    io = IOBuffer()
    for s in sort(surfaces)
        types,tags,nodes = gmsh.model.mesh.get_elements(2,s)
        print(io,s,types,tags,nodes)
        nt,xyz,_ = gmsh.model.mesh.get_nodes(2,s,true,false)
        order = sortperm(nt)
        print(io,nt[order],reshape(xyz,3,:)[:,order])
    end
    bytes2hex(sha256(take!(io)))
end
count_elements(surfaces) = sum(length(tags) for s in surfaces for tags in gmsh.model.mesh.get_elements(2,s)[2];init=0)

function inspect_mesh()
    corners = reduce(vcat,corner_groups())
    pml = last.(entities(2,1005))
    physical = setdiff(last.(gmsh.model.get_entities(2)),pml)
    strips = setdiff(pml,corners)
    paths = [e for tag in (7001,7002) for e in entities(1,tag)]
    return Dict{String,Any}("nodes"=>length(first(gmsh.model.mesh.get_nodes())),
        "corner_triangles"=>count_elements(corners),"strip_triangles"=>count_elements(strips),
        "pml_triangles"=>count_elements(pml),"physical_triangles"=>count_elements(physical),
        "total_triangles"=>count_elements(last.(gmsh.model.get_entities(2))),
        "physical_hash"=>mesh_hash(physical),"strip_hash"=>mesh_hash(strips),
        "path_hash"=>coordinate_hash(paths))
end

function coordinate_hash(es)
    points = Set{NTuple{3,Float64}}()
    for (d,t) in es
        _,xyz,_ = gmsh.model.mesh.get_nodes(d,t,true,false)
        union!(points,Tuple.(eachcol(reshape(xyz,3,:))))
    end
    bytes2hex(sha256(codeunits(repr(sort!(collect(points))))))
end

function modal_value(Q,u,b,slope)
    z, endz = FEM._pml_coordinate(u,b,slope), FEM._pml_coordinate(1.,b,slope)
    iszero(Q) && return 1-z/endz
    return exp(-Q*z)*(-expm1(-2Q*(endz-z)))/(-expm1(-2Q*endz))
end

function interpolate(x, xs, ys)
    i = clamp(searchsortedlast(xs,x),1,length(xs)-1)
    t = (x-xs[i])/(xs[i+1]-xs[i])
    return (1-t)*ys[i]+t*ys[i+1]
end

function density_samples(direction, scale)
    f = value("Frequencies")
    gamma = complex(value("GammaReValues"),value("GammaImValues"))
    q = map(((value("AirMu"),0.,value("AirEpsilon")),
             (value("EarthMu"),value("EarthSigma"),value("EarthEpsilon")))) do (mu,sigma,eps)
        sqrt(complex(im*2pi*f*mu*(sigma+im*2pi*f*eps)-gamma^2))
    end
    axis = ("Side","Top","Bottom")[direction]
    L,b = value("Pml$(axis)ThicknessValues"),value("Pml$(axis)StrengthValues")
    slope = value("PmlSlopeValues")
    # Conservative clearance from a circle enclosing the two-wire bounding
    # box. The actual Cartesian wall clearance is larger. Boundary nodes
    # retain the original production prescription exactly.
    clearance = value("DomainHalfwidths")-hypot(.5,1.0425)
    grid = FEM._pml_density_grid(q,L,b,clearance,direction,B.MESH.pml_resolution;slope)
    # Reconstruct the prescribed density from its cumulative integral.
    mid = (grid.u[1:end-1]+grid.u[2:end])/2
    density = diff(grid.integral)./diff(grid.u)
    modes = FEM._pml_modes(q...,L,clearance,direction)
    sample_u = sort!(unique([FEM._pml_equal_density_nodes(grid,60);min(1.,b^(-1/3))]))
    # The cross-coordinate multiplier comes from the maximum finite-layer
    # modal amplitude, including the static mode. It is deliberately an
    # envelope, not a fit to the terminal conductance. Squared interpolation
    # contributions give amplitude^(2/3) in the 1D density estimate.
    envelope(u) = maximum(sqrt(m.weight)*abs(modal_value(m.Q,u,b,slope)) for m in modes)
    return (;L,b,slope,u=sample_u,
        density=[interpolate(u,mid,density)/L/scale for u in sample_u],
        amplitude=[envelope(u) for u in sample_u])
end

function corner_metric!(dir, scale)
    axes = ntuple(d->density_samples(d,scale),3)
    D,xc = value("DomainHalfwidths"), value("Xcenter")
    data = Float64[]
    for side in (-1,1), vertical in (-1,1)
        ax,ay = axes[1],axes[vertical>0 ? 2 : 3]
        # Bound the relaxation where a finite-wall mode vanishes. This cap is
        # shared by both candidates; only the fixed metric scale differs.
        function vertex(i,j)
            x,y = xc+side*(D+ax.L*ax.u[i]), vertical*(D+ay.L*ay.u[j])
            wx,wy = max(ax.amplitude[i],.125)^(2/3),max(ay.amplitude[j],.125)^(2/3)
            rx,ry = ax.density[i]*wy,ay.density[j]*wx
            return (x,y,0.),(rx^2,0.,0.,0.,ry^2,0.,0.,0.,max(rx^2,ry^2))
        end
        for j in 1:length(ay.u)-1, i in 1:length(ax.u)-1
            for tri in (((i,j),(i+1,j),(i+1,j+1)),((i,j),(i+1,j+1),(i,j+1)))
                v = [vertex(ij...) for ij in tri]
                # Gmsh list data: x1 x2 x3 y1 y2 y3 z1 z2 z3, then tensors.
                append!(data,[p[1][k] for k in 1:3 for p in v])
                for p in v; append!(data,p[2]); end
            end
        end
    end
    view = gmsh.view.add("prescribed PML corner metric")
    gmsh.view.add_list_data(view,"TT",length(data)÷36,data)
    gmsh.view.write(view,joinpath(dir,"corner_metric.pos"))
    field = gmsh.model.mesh.field.add("PostView")
    gmsh.model.mesh.field.set_number(field,"ViewTag",view)
    gmsh.model.mesh.field.set_as_background_mesh(field)
    gmsh.option.set_number("Mesh.AnisoMax",1e6)
    gmsh.option.set_number("Mesh.SmoothRatio",1.8)
    gmsh.option.set_number("Mesh.MeshSizeFromPoints",0)
    gmsh.option.set_number("Mesh.MeshSizeExtendFromBoundary",0)
    return view
end

function replace_corners!(target,candidate)
    reference_model = gmsh.model.get_current()
    groups = corner_groups()
    old_corners = reduce(vcat,groups)
    memberships = Dict(tag=>Int.(last.(entities(2,tag))) for tag in (1003,1004,1005))
    boundaries = [Int.(last.(gmsh.model.get_boundary([(2,s) for s in ss],true,true,false))) for ss in groups]
    all_edges = unique(abs.(last.(gmsh.model.get_boundary([(2,s) for s in old_corners],false,true,false))))
    internal_edges = setdiff(all_edges,abs.(reduce(vcat,boundaries)))
    all_points = unique(last.(gmsh.model.get_boundary([(1,c) for c in all_edges],false,false,false)))
    border_points = unique(last.(gmsh.model.get_boundary([(1,c) for c in unique(abs.(reduce(vcat,boundaries)))],false,false,false)))
    internal_points = setdiff(all_points,border_points)
    # Mesh corners in separate native Gmsh models. Importing an old mesh into
    # independently exported CAD can change entity classification; this avoids
    # that confound and preserves all existing noncorner elements exactly.
    tags,coords,_ = gmsh.model.mesh.get_nodes(-1,-1,false,false)
    xyz = Dict(Int(t)=>Tuple(v) for (t,v) in zip(tags,eachcol(reshape(coords,3,:))))
    next_node = maximum(tags)+1
    next_element = maximum(reduce(vcat,gmsh.model.mesh.get_elements()[2]))+1
    metric_view = corner_metric!(target,candidate.scale)
    meshes = []
    for (index,boundary) in enumerate(boundaries)
        segments = [(Int(v[1]),Int(v[2])) for c in abs.(boundary)
            for v in eachcol(reshape(only(gmsh.model.mesh.get_elements(1,c)[3]),2,:))]
        boundary_tags = unique(reduce(vcat,([a,b] for (a,b) in segments)))
        points = [xyz[t] for t in boundary_tags]
        xmin,xmax = extrema(first.(points)); ymin,ymax = extrema(p[2] for p in points)
        corners = [(xmin,ymin,0.),(xmax,ymin,0.),(xmax,ymax,0.),(xmin,ymax,0.)]
        corner_tags = [only([t for t in boundary_tags if xyz[t]==p]) for p in corners]
        model_name = "corner-$index"
        gmsh.model.add(model_name)
        p = [gmsh.model.geo.add_point(v...) for v in corners]
        curves = [gmsh.model.geo.add_line(p[i],p[mod1(i+1,4)]) for i in 1:4]
        surface = gmsh.model.geo.add_plane_surface([gmsh.model.geo.add_curve_loop(curves)])
        gmsh.model.geo.synchronize()
        for i in 1:4
            gmsh.model.mesh.add_nodes(0,p[i],[corner_tags[i]],collect(corners[i]))
            gmsh.model.mesh.add_elements_by_type(p[i],15,Int[],[corner_tags[i]])
        end
        for i in 1:4
            a,b = corners[i],corners[mod1(i+1,4)]
            axis = a[1]==b[1] ? 2 : 1
            fixed_axis = 3-axis
            side_nodes = sort([t for t in boundary_tags if xyz[t][fixed_axis]==a[fixed_axis]];
                by=t->(xyz[t][axis]-a[axis])/(b[axis]-a[axis]))
            middle = side_nodes[2:end-1]
            gmsh.model.mesh.add_nodes(1,curves[i],middle,reduce(vcat,(collect(xyz[t]) for t in middle);init=Float64[]),
                [(xyz[t][axis]-a[axis])/(b[axis]-a[axis]) for t in middle])
            gmsh.model.mesh.add_elements_by_type(curves[i],1,Int[],
                [side_nodes[k] for j in 1:length(side_nodes)-1 for k in (j,j+1)])
        end
        gmsh.model.mesh.set_algorithm(2,surface,7)
        field = gmsh.model.mesh.field.add("PostView")
        gmsh.model.mesh.field.set_number(field,"ViewTag",metric_view)
        gmsh.model.mesh.field.set_as_background_mesh(field)
        gmsh.option.set_number("Mesh.MeshOnlyEmpty",1)
        gmsh.option.set_number("Mesh.Renumber",0)
        gmsh.model.mesh.generate(2)
        nt,xc,_ = gmsh.model.mesh.get_nodes(2,surface,false,false)
        triangle_nodes = only(gmsh.model.mesh.get_elements(2,surface)[3])
        mapping = Dict(t=>t for t in boundary_tags)
        new_tags = collect(next_node:next_node+length(nt)-1)
        for (t,new) in zip(nt,new_tags); mapping[Int(t)]=Int(new); end
        next_node += length(nt)
        new_elements = collect(next_element:next_element+length(triangle_nodes)÷3-1)
        next_element += length(new_elements)
        push!(meshes,(new_tags,xc,new_elements,[mapping[Int(t)] for t in triangle_nodes]))
        gmsh.model.remove()
        gmsh.model.set_current(reference_model)
    end
    gmsh.view.remove(metric_view)
    gmsh.model.mesh.clear([(2,s) for s in old_corners])
    gmsh.model.mesh.clear([(1,c) for c in internal_edges])
    gmsh.model.mesh.clear([(0,p) for p in internal_points])
    gmsh.model.remove_entities([(2,s) for s in old_corners],false)
    gmsh.model.remove_entities([(1,c) for c in internal_edges],false)
    gmsh.model.remove_entities([(0,p) for p in internal_points],false)
    new_corners = Int[]
    for (boundary,(nt,xc,et,en)) in zip(boundaries,meshes)
        s = gmsh.model.add_discrete_entity(2,-1,boundary)
        gmsh.model.mesh.add_nodes(2,s,nt,xc)
        gmsh.model.mesh.add_elements_by_type(s,2,et,en)
        push!(new_corners,s)
    end
    gmsh.model.remove_physical_groups([(2,t) for t in keys(memberships)])
    for tag in (1003,1004,1005)
        add = tag==1003 ? new_corners[3:4] : tag==1004 ? new_corners[1:2] : new_corners
        gmsh.model.add_physical_group(2,[setdiff(memberships[tag],old_corners);add],tag)
    end
    return new_corners
end

function mesh_candidate!(baseline, target, candidate)
    marker = joinpath(target,"mesh.toml")
    if isfile(marker)
        d = TOML.parsefile(marker)
        digest(joinpath(target,"study.msh"))==d["sha256"] || error("Changed candidate mesh")
        return d
    end
    copy_bundle!(baseline,target)
    session = FEM._start_gmsh(4)
    try
        gmsh.parser.clear(); gmsh.onelab.clear()
        gmsh.parser.set_string("OnelabAction",["audit"])
        gmsh.open(joinpath(baseline,"study.msh"))
        gmsh.parser.parse(joinpath(target,"study_data.pro"))
        before = inspect_mesh()
        record(joinpath(baseline,"corner-audit.toml"),before)
        say("CORNER AUDIT ",basename(dirname(baseline))," ",before)
        measured = @timed begin
            replace_corners!(target,candidate)
        end
        after = inspect_mesh()
        for key in ("physical_hash","strip_hash","path_hash")
            before[key]==after[key] || error("Mesh isolation failed: $key")
        end
        after["mesh_seconds"] = measured.time
        after["compile_seconds"] = measured.compile_time
        after["metric_scale"] = candidate.scale
        gmsh.write(joinpath(target,"study.msh"))
        after["sha256"] = digest(joinpath(target,"study.msh"))
        record(marker,after)
        say("CANDIDATE MESH ",candidate.name," ",after)
        return after
    finally
        gmsh.option.set_number("Mesh.MeshOnlyEmpty",0)
        FEM._finish_gmsh(session)
    end
end

function run_pilots(;mesh_only=false)
    mkpath(ROOT)
    for case in PILOTS
        baseline = freeze!(case)
        problem,form = B.fixture(case...)
        mesh_only || B.solve!(baseline)
        for candidate in CANDIDATES
            target = joinpath(dirname(baseline),candidate.name)
            say("PREPARE ",label(case)," ",candidate.name)
            mesh_candidate!(baseline,target,candidate)
            if !mesh_only
                B.solve!(target)
                B.compare!(dirname(baseline),problem,form;candidate=candidate.name,tolerance=.02)
            end
        end
    end
    say("COMPLETE corner pilot stage; inspect comparisons before any production change")
end

if abspath(PROGRAM_FILE)==@__FILE__
    run_pilots(;mesh_only="mesh-only" in ARGS)
    if !("mesh-only" in ARGS)
        include("audit_meshes.jl")
        audit_candidates()
        Base.include(Module(:CornerAssessment),joinpath(@__DIR__,"assess.jl"))
    end
end
