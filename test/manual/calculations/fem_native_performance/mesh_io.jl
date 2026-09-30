isdefined(@__MODULE__,:PILOTS) || include("qualify.jl")

function mesh_io!()
    target=mkpath(joinpath(ROOT,"mesh-io"))
    records=Dict{String,Any}()
    for case in PILOTS
        source=joinpath(ROOT,label(case),"baseline/study.msh")
        session=FEM._start_gmsh(0)
        try
            gmsh.open(source)
            before=entity_hash(gmsh.model.get_entities())
            groups=sort(gmsh.model.get_physical_groups())
            membership=[(d,t,sort(gmsh.model.get_entities_for_physical_group(d,t))) for (d,t) in groups]
            for binary in (0,1)
                path=joinpath(target,"$(label(case))-$binary.msh")
                gmsh.option.set_number("Mesh.Binary",binary)
                gmsh.option.set_number("Mesh.MshFileVersion",4.1)
                gmsh.option.set_number("Mesh.SaveAll",1)
                writes=Float64[]; reads=Float64[]
                # Initial call excluded from warm statistics.
                for repeat in 0:5
                    w=@elapsed gmsh.write(path)
                    gmsh.clear()
                    r=@elapsed gmsh.open(path)
                    entity_hash(gmsh.model.get_entities())==before || error("Mesh data changed on round trip")
                    [(d,t,sort(gmsh.model.get_entities_for_physical_group(d,t))) for (d,t) in groups]==membership ||
                        error("Physical group membership changed")
                    repeat==0 || (push!(writes,w);push!(reads,r))
                end
                data=Dict("binary"=>binary,"bytes"=>filesize(path),
                    "median_write_seconds"=>median(writes),"median_read_seconds"=>median(reads),
                    "all_mesh_values_preserved"=>true)
                records["$(label(case))-$binary"]=data
                say("MESH IO ",label(case)," ",data)
            end
        finally
            FEM._finish_gmsh(session)
        end
    end
    record(joinpath(target,"results.toml"),records)
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && mesh_io!()
