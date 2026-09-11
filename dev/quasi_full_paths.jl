# Voltage-path input for the manual quasi-full.pro experiment.
# Uses the existing first-order triangular mesh; no field solve is performed here.
using Gmsh
# Clip a straight segment to each triangle. Lowest-order edge elements have
# constant tangential component on a straight line, so one midpoint is exact.
function quasi_full_path_points(triangles, vertices)
    points = NTuple{4,Float64}[]
    cross2(a,b) = a[1]*b[2]-a[2]*b[1]
    for (p,q) in zip(vertices[1:end-1],vertices[2:end])
        d = q .- p
        intervals = Tuple{Float64,Float64}[]
        for triangle in triangles
            minimum(v[1] for v in triangle) > max(p[1],q[1]) && continue
            maximum(v[1] for v in triangle) < min(p[1],q[1]) && continue
            minimum(v[2] for v in triangle) > max(p[2],q[2]) && continue
            maximum(v[2] for v in triangle) < min(p[2],q[2]) && continue
            lo,hi = 0.0,1.0
            signarea = sign(cross2(triangle[2].-triangle[1],triangle[3].-triangle[1]))
            for k in 1:3
                u,v = triangle[k],triangle[mod1(k+1,3)]
                edge = v.-u
                a,b = signarea*cross2(edge,p.-u),signarea*cross2(edge,d)
                if iszero(b)
                    a < -1e-12*max(sum(abs2,edge),1e-12) && (hi = -1.0)
                elseif b > 0
                    lo = max(lo,-a/b)
                else
                    hi = min(hi,-a/b)
                end
            end
            if hi-lo > 1e-11
                push!(intervals, (lo, hi))
            end
        end
        # The shared tangential trace is unique on a mesh edge. Merge duplicate
        # intervals there, but retain every triangle crossing elsewhere.
        sort!(intervals)
        covered = 0.0
        for (lo, hi) in intervals
            lo <= covered + 1e-8 || error("Voltage path leaves the field domain at $p -> $q")
            lo = max(lo, covered)
            hi-lo > 1e-11 || continue
            mid = p .+ ((lo+hi)/2).*d
            weight = (hi-lo).*d
            push!(points, (mid[1],mid[2],weight[1],weight[2]))
            covered = hi
        end
        abs(covered-1) < 1e-8 || error("Incomplete voltage path on $p -> $q")
    end
    points
end

function write_quasi_full_paths(path, mesh, plan, model, positions, radius; shell_segments=128)
    shell_segments >= 2 || throw(ArgumentError("shell_segments must be at least 2"))
    length(positions) == length(model.terminal_ids) || throw(ArgumentError("one path per terminal is required"))
    gmsh.open(mesh)
    tags, coords, _ = gmsh.model.mesh.get_nodes()
    nodes = Dict(tag => (coords[3i-2],coords[3i-1]) for (i,tag) in enumerate(tags))
    triangles = NTuple{3,NTuple{2,Float64}}[]
    metal_tags = [m.physical_tag for m in model.material_plans if m.kind === :conductor]
    for (_, surface) in gmsh.model.get_entities(2)
        groups = gmsh.model.get_physical_groups_for_entity(2,surface)
        any(tag -> tag in metal_tags, groups) && continue
        types, _, element_nodes = gmsh.model.mesh.get_elements(2,surface)
        for (type, nds) in zip(types,element_nodes)
            type == 2 || error("Path integration expects first-order triangles")
            append!(triangles,[(nodes[nds[i]],nodes[nds[i+1]],nodes[nds[i+2]]) for i in 1:3:length(nds)])
        end
    end
    allpoints = NTuple{4,Float64}[]
    starts, counts = Int[],Int[]
    ri,ro = plan.domain_radius,plan.shell_outer_radius
    cx,cy = model.centre
    for (x,y) in positions
        dx = x-cx
        # Pull the physical vertical ray back through GetDP's spherical shell:
        # r_phys = ri*(ro-ri)/(ro-r_mesh), t = ri/r_phys in [0,1].
        vertices = [(cx+(ro-(ro-ri)*t)*t*dx/ri,
            cy-(ro-(ro-ri)*t)*sqrt(1-(t*dx/ri)^2)) for t in range(0,1;length=shell_segments+1)]
        push!(vertices,(x,y-radius))
        hypot(x-cx, y-radius-cy) < ri || error("Voltage endpoint must lie in the finite field domain")
        pts = quasi_full_path_points(triangles,vertices)
        push!(starts,length(allpoints));push!(counts,length(pts));append!(allpoints,pts)
    end
    open(path,"w") do io
        for (name,values) in (("PathStart",starts),("PathCount",counts),
            ("PathX",first.(allpoints)),("PathY",getindex.(allpoints,2)),
            ("PathDX",getindex.(allpoints,3)),("PathDY",last.(allpoints)))
            println(io,name,"() = {",join(values,","),"};")
        end
    end
    println("Voltage path integration points per terminal: ",counts)
end
