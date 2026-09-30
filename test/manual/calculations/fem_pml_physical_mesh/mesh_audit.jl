using Gmsh, TOML, SHA, Dates
const gmsh=Gmsh.gmsh
function inspect_mesh(path)
    gmsh.initialize();gmsh.option.set_number("General.Terminal",0)
    try
        gmsh.open(path)
        ntags,xyz,_=gmsh.model.mesh.get_nodes()
        coordinates=Dict(t=>Tuple(xyz[3i-2:3i]) for (i,t) in enumerate(ntags))
        groups=Dict{String,Any}()
        for (dim,tag) in gmsh.model.get_physical_groups()
            name=gmsh.model.get_physical_name(dim,tag)
            occursin("field_maps",name) && continue
            nodes=Set{UInt64}();elements=0
            for entity in gmsh.model.get_entities_for_physical_group(dim,tag)
                _,etags,enodes=gmsh.model.mesh.get_elements(dim,entity)
                elements+=sum(length,etags;init=0)
                for tags in enodes;union!(nodes,tags);end
            end
            points=sort!([coordinates[t] for t in nodes])
            groups[name]=(;elements,nodes=length(points),hash=bytes2hex(sha256(repr(points))))
        end
        (;nodes=length(ntags),groups)
    finally
        gmsh.finalize()
    end
end
root=joinpath(pwd(),".linecablemodels/fem/pml-physical-mesh")
candidate=inspect_mesh(joinpath(root,"step-a-low/detached/study.msh"))
reference=inspect_mesh(joinpath(dirname(root),"pml-conductance-cost/phase2-balanced144-low/detached/study.msh"))
open(joinpath(root,"mesh-audit.csv"),"w") do io
    println(io,"group,candidate_nodes,reference_nodes,candidate_elements,reference_elements,identical_coordinates")
    for name in sort(collect(keys(candidate.groups)))
        a=candidate.groups[name];haskey(reference.groups,name) || continue;b=reference.groups[name]
        println(io,join((name,a.nodes,b.nodes,a.elements,b.elements,a.hash==b.hash),','))
    end
end
println(Dates.now()," MESH AUDIT: candidate_nodes=",candidate.nodes," reference_nodes=",reference.nodes)
