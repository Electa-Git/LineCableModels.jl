# Post-run mesh integrity checks only; no additional PDE solves.
isdefined(@__MODULE__,:PILOTS) || include("qualify.jl")
function audit_candidates()
    session = FEM._start_gmsh(1)
    try
        for case in PILOTS, candidate in CANDIDATES
            dir = joinpath(ROOT,label(case),candidate.name)
            isfile(joinpath(dir,"mesh.toml")) || continue
            gmsh.clear(); gmsh.parser.clear(); gmsh.onelab.clear()
            # Candidate MSH output can omit ungrouped 1D elements. The frozen
            # reference retains all perimeter segments and is authoritative.
            gmsh.open(joinpath(dirname(dir),"baseline/study.msh"))
            gmsh.parser.parse(joinpath(dir,"study_data.pro"))
            expected_perimeters = map(corner_groups()) do surfaces
                boundary = abs.(last.(gmsh.model.get_boundary([(2,s) for s in surfaces],true,true,false)))
                Set(minmax(v...) for curve in boundary
                    for v in eachcol(reshape(only(gmsh.model.mesh.get_elements(1,curve)[3]),2,:)))
            end
            gmsh.clear(); gmsh.parser.clear(); gmsh.onelab.clear()
            gmsh.open(joinpath(dir,"study.msh"))
            gmsh.parser.parse(joinpath(dir,"study_data.pro"))
            tags,coords,_ = gmsh.model.mesh.get_nodes(-1,-1,false,false)
            xyz = Dict(Int(t)=>Tuple(v) for (t,v) in zip(tags,eachcol(reshape(coords,3,:))))
            areas = Float64[]; corner_count=0
            for (index,surfaces) in enumerate(corner_groups())
                @assert length(surfaces)==1
                s = only(surfaces); corner_count+=1
                perimeter = Dict{Tuple{Int,Int},Int}()
                area = 0.
                triangles = reshape(only(gmsh.model.mesh.get_elements(2,s)[3]),3,:)
                for (a,b,c) in eachcol(triangles)
                    pa,pb,pc = xyz[a],xyz[b],xyz[c]
                    signed = ((pb[1]-pa[1])*(pc[2]-pa[2])-(pb[2]-pa[2])*(pc[1]-pa[1]))/2
                    @assert signed>0
                    area += signed
                    for pair in ((a,b),(b,c),(c,a))
                        edge = minmax(pair...)
                        perimeter[edge] = get(perimeter,edge,0)+1
                    end
                end
                @assert all(n in (1,2) for n in values(perimeter))
                outside = Set(e for (e,n) in perimeter if n==1)
                expected = expected_perimeters[index]
                @assert outside==expected
                border_nodes = unique([t for pair in expected for t in pair])
                x0,x1 = extrema(xyz[t][1] for t in border_nodes)
                y0,y1 = extrema(xyz[t][2] for t in border_nodes)
                expected_area = (x1-x0)*(y1-y0)
                @assert isapprox(area,expected_area;rtol=1e-8)
                push!(areas,area)
            end
            @assert corner_count==4
            record(joinpath(dir,"integrity.toml"),Dict("passed"=>true,
                "positive_signed_areas"=>true,"conforming_boundary_connectivity"=>true,
                "corner_areas"=>areas))
            say("MESH INTEGRITY PASS ",label(case),"/",candidate.name)
        end
    finally
        FEM._finish_gmsh(session)
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    audit_candidates()
end
