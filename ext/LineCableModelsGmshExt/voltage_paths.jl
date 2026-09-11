# Clip a straight segment to each triangle. Lowest-order edge elements have
# constant tangential component on a straight line, so one midpoint is exact.
function _voltage_path_points(triangles, vertices; active=trues(length(triangles)))
    length(active) == length(triangles) || throw(DimensionMismatch("one activity flag per triangle"))
    points = NTuple{4,Float64}[]
    cross2(a,b) = a[1]*b[2]-a[2]*b[1]
    for (p,q) in zip(vertices[1:end-1],vertices[2:end])
        d = q .- p
        intervals = Tuple{Float64,Float64,Bool}[]
        for (triangle, integrate) in zip(triangles, active)
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
                push!(intervals, (lo, hi, integrate))
            end
        end
        # The shared tangential trace is unique on a mesh edge. Merge duplicate
        # intervals there, but retain every triangle crossing elsewhere.
        sort!(intervals)
        covered = 0.0
        for (lo, hi, integrate) in intervals
            lo <= covered + 1e-8 || error("Voltage path leaves the field domain at $p -> $q")
            lo = max(lo, covered)
            hi-lo > 1e-11 || continue
            mid = p .+ ((lo+hi)/2).*d
            weight = (hi-lo).*d
            integrate && push!(points, (mid[1],mid[2],weight[1],weight[2]))
            covered = hi
        end
        abs(covered-1) < 1e-8 || error("Incomplete voltage path on $p -> $q")
    end
    points
end

# Choose an actual contour node, so endpoints lie on the discretized electrode
# for circular, tubular, stranded, and noncircular terminal geometries alike.
function _voltage_endpoint(model, terminal, nodes)
    curves = gmsh.model.get_entities_for_physical_group(
        1, model.tags.terminal_contour_base + terminal)
    tags = unique!(reduce(vcat, [first(gmsh.model.mesh.get_nodes(1, curve, true))
        for curve in curves]; init=UInt64[]))
    isempty(tags) && error("Terminal $(model.terminal_ids[terminal]) has no contour nodes")
    return nodes[argmin(tag -> (nodes[tag][2], nodes[tag][1]), tags)]
end

function _write_voltage_paths(path, mesh, plan, model; endpoints=nothing, shell_segments=128)
    shell_segments >= 2 || throw(ArgumentError("shell_segments must be at least 2"))
    previous_model = gmsh.model.get_current()
    scratch_model = "LineCableModels-voltage-$(basename(tempname()))"
    gmsh.model.add(scratch_model)
    try
        # Merge into a uniquely named empty model: gmsh.open names its model
        # after the mesh file and can collide with a caller-owned model.
        gmsh.merge(mesh)
        tags, coords, _ = gmsh.model.mesh.get_nodes()
        nodes = Dict(tag => (coords[3i-2], coords[3i-1]) for (i, tag) in enumerate(tags))
        triangles = NTuple{3,NTuple{2,Float64}}[]
        active = Bool[]
        metal_tags = Set(m.physical_tag for m in model.material_plans if m.kind === :conductor)
        for (_, surface) in gmsh.model.get_entities(2)
            groups = gmsh.model.get_physical_groups_for_entity(2, surface)
            integrate = !any(in(metal_tags), groups)
            types, _, element_nodes = gmsh.model.mesh.get_elements(2, surface)
            for (type, nds) in zip(types, element_nodes)
                type == 2 || error("Voltage-path integration requires first-order triangles")
                for i in 1:3:length(nds)
                    push!(triangles, (nodes[nds[i]], nodes[nds[i+1]], nodes[nds[i+2]]))
                    push!(active, integrate)
                end
            end
        end
        endpoints === nothing && (endpoints = [_voltage_endpoint(model, i, nodes)
            for i in eachindex(model.terminal_ids)])
        length(endpoints) == length(model.terminal_ids) ||
            throw(DimensionMismatch("one voltage endpoint per terminal is required"))
        allpoints = NTuple{4,Float64}[]
        starts, counts = Int[], Int[]
        ri, ro = plan.domain_radius, plan.shell_outer_radius
        cx, cy = model.centre
        for (x, y) in endpoints
            hypot(x-cx, y-cy) < ri || error("Voltage endpoint must lie in the finite field domain")
            dx = x-cx
            # Pull back a physical vertical ray from earth infinity through the
            # spherical shell: r_phys = ri*(ro-ri)/(ro-r_mesh), t=ri/r_phys.
            vertices = [(cx+(ro-(ro-ri)*t)*t*dx/ri,
                cy-(ro-(ro-ri)*t)*sqrt(1-(t*dx/ri)^2))
                for t in range(0, 1; length=shell_segments+1)]
            push!(vertices, (x, y))
            # The electric model extends v constantly and bt=0 inside each
            # equipotential metal. Vertical paths can therefore cross shields
            # or other terminals without introducing a transverse voltage there.
            points = _voltage_path_points(triangles, vertices; active)
            push!(starts, length(allpoints))
            push!(counts, length(points))
            append!(allpoints, points)
        end
        open(path, "w") do io
            println(io, "// Physical vertical voltage paths; zero transverse field inside metal.")
            for (name, values) in (("PathStart", starts), ("PathCount", counts),
                ("PathX", first.(allpoints)), ("PathY", getindex.(allpoints, 2)),
                ("PathDX", getindex.(allpoints, 3)), ("PathDY", last.(allpoints)))
                println(io, name, "() = ", _pro_array(values), ";")
            end
        end
        return path
    finally
        gmsh.model.set_current(scratch_model)
        gmsh.model.remove()
        isempty(previous_model) || gmsh.model.set_current(previous_model)
    end
end

_voltage_path_file(run::FEMRun, frequency::Int) =
    joinpath(run.path, "input", @sprintf("paths-f%04d.pro", frequency))

function _prepare_voltage_paths!(run, model, formulation, mesh_paths)
    _quasi_full(formulation.options.physics) || return nothing
    @info "Preparing quasi-full voltage paths" frequencies=length(mesh_paths)
    for (mesh, plan) in zip(mesh_paths, model.mesh_plans)
        _write_voltage_paths(_voltage_path_file(run, plan.frequency_index), mesh, plan, model)
    end
    return nothing
end
