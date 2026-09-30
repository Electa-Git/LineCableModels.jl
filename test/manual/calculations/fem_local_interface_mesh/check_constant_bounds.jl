# Inspect the failed finite-Gamma case without meshing or solving it again.
include("qualify.jl")
problem, _ = fixture(:mixed,1e6,.99)
dir = joinpath(ROOT,"mixed-f1.0e6-gamma0.99")
session = FEM._start_gmsh(0)
try
    gmsh.parser.clear(); gmsh.onelab.clear()
    gmsh.parser.set_string("OnelabAction",["audit"])
    gmsh.open(joinpath(dir,"baseline/study.geo"))
    interface = Set(gmsh.model.get_entities_for_physical_group(1,2002))
    distances = [f for f in gmsh.model.mesh.field.list()
        if gmsh.model.mesh.field.get_type(f)=="Distance" &&
            !isempty(intersect(interface,gmsh.model.mesh.field.get_numbers(f,"CurvesList")))]
    records = Dict{String,Any}()
    for (medium,distance,tag) in zip(("air","soil"),distances,(1001,1002))
        threshold = only(f for f in gmsh.model.mesh.field.list()
            if gmsh.model.mesh.field.get_type(f)=="Threshold" &&
                gmsh.model.mesh.field.get_number(f,"InField")==distance)
        hmin = gmsh.model.mesh.field.get_number(threshold,"SizeMin")
        hmax = gmsh.model.mesh.field.get_number(threshold,"SizeMax")
        decay = gmsh.model.mesh.field.get_number(threshold,"DistMin")
        boxes = [gmsh.model.get_bounding_box(2,s) for s in gmsh.model.get_entities_for_physical_group(2,tag)]
        xmin,xmax = minimum(b[1] for b in boxes),maximum(b[4] for b in boxes)
        ymin,ymax = minimum(b[2] for b in boxes),maximum(b[5] for b in boxes)
        # Distance to one segment is convex; its maximum over a rectangle is
        # attained at a corner. The nearest of several segments is no farther.
        bound = minimum(maximum(hypot(y,max(abs(x-position.x)-abs(position.y)-LineCableModels.outer_radius(design),0.))
            for x in (xmin,xmax), y in (ymin,ymax))
            for (design,position) in zip(problem.system.designs,problem.system.positions))
        unchanged = hmin==hmax || bound<=decay
        records[medium] = Dict("size_min"=>hmin,"size_max"=>hmax,"decay_radius"=>decay,
            "compact_distance_upper_bound"=>bound,"unchanged_wave_target_over_physical_region"=>unchanged,
            "bounding_box"=>[xmin,xmax,ymin,ymax])
        say("CONSTANT BOUND ",medium," ",records[medium])
    end
    record(joinpath(dir,"constant-target-bounds.toml"),records)
finally
    FEM._finish_gmsh(session)
end
