"""
$(TYPEDSIGNATURES)

Write the meshes of `model` at every frequency \\[Hz\\] under `directory`, with at most
`options.workers` native Gmsh processes at once. Keep each process's log under `logs/`.
After every process ends, validate each mesh and add its gauge tree in the caller's
initialized Gmsh session. Return the mesh paths in frequency order.

# Errors

- Throws `ErrorException` for a failed process or mesh check.
- Rethrows `InterruptException` after the pool stops its processes.

The error identifies the frequency, run directory, log path, and last twenty log lines.
A process failure includes its exit code and signal. The run directory stays on failure.
"""
function mesh(directory::AbstractString, model::FieldModel,
        options::ComputationOptions{GmshBackend})
    for name in ("work", "logs", "mesh")
        mkpath(joinpath(directory, name))
    end
    processes = Base.Process[]
    stopping = Ref(false)
    failures = Pair{Int, Exception}[]
    pool = @async asyncmap(eachindex(model.problem.frequencies); ntasks = options.workers) do k
        stopping[] && return
        try
            stem = @sprintf("mesh-%04d", k)
            folder = joinpath(directory, "work", stem)
            mkpath(folder)
            for source in readdir(joinpath(directory, "input", "getdp"); join = true)
                cp(source, joinpath(folder, basename(source)))
            end
            for name in ("model_data.pro", "physical.geo")
                cp(joinpath(directory, "input", name), joinpath(folder, name))
            end
            entry = "Include \"model_data.pro\";\n" *
                    "Include \"parameters.pro\";\n" *
                    "Include \"physical.geo\";\n" *
                    "Include \"geometry.geo\";\n" *
                    "Include \"mesh.geo\";\n"
            write(joinpath(folder, "model.geo"), entry)
            log_path = joinpath(directory, "logs", stem * ".log")
            command = Cmd(
                `$(Gmsh.gmsh_jll.gmsh()) model.geo -setnumber FrequencyIndex $k -setstring MeshMetadataPath mesh-metadata.txt -2 -o model.msh -v $(options.gmsh_verbosity)`;
                dir = folder)
            io = open(log_path, "w")
            try
                process = run(pipeline(command; stdout = io, stderr = io); wait = false)
                push!(processes, process)
                wait(process)
                if !success(process) && !stopping[]
                    error("exit code $(process.exitcode), signal $(process.termsignal)")
                end
            finally
                close(io)
            end
        catch exception
            if !stopping[]
                if exception isa InterruptException
                    push!(failures, k => exception)
                else
                    log_path = joinpath(directory, "logs", @sprintf("mesh-%04d.log", k))
                    tail = isfile(log_path) ? join(last(readlines(log_path), 20), '\n') : ""
                    message = "mesh: frequency $(model.problem.frequencies[k]) Hz, " *
                              "run directory $directory, log $log_path\n" *
                              sprint(showerror, exception) * "\n$tail"
                    push!(failures, k => ErrorException(message))
                end
                stopping[] = true
            end
        end
        return
    end
    try
        while !istaskdone(pool) && !stopping[]
            sleep(0.02)
        end
    finally
        stopping[] = true
        deadline = time() + 2
        terminated = Set{Base.Process}()
        while !istaskdone(pool) || any(process_running, processes)
            for process in processes
                if process_running(process) &&
                   (!(process in terminated) || time() >= deadline)
                    kill(process, time() < deadline ? Base.SIGTERM : Base.SIGKILL)
                    push!(terminated, process)
                end
            end
            sleep(0.02)
        end
        for process in processes
            wait(process)
        end
        wait(pool)
    end
    isempty(failures) || throw(last(first(failures)))
    paths = String[]
    for index in eachindex(model.problem.frequencies)
        stem = @sprintf("mesh-%04d", index)
        folder = joinpath(directory, "work", stem)
        path = joinpath(folder, "model.msh")
        log_path = joinpath(directory, "logs", stem * ".log")
        try
            validate(path, model)
            gauge!(path, model)
        catch exception
            exception isa InterruptException && rethrow()
            tail = join(last(readlines(log_path), 20), '\n')
            error("mesh: frequency $(model.problem.frequencies[index]) Hz, " *
                  "run directory $directory, log $log_path\n" *
                  sprint(showerror, exception) * "\n$tail")
        end
        destination = joinpath(directory, "mesh", @sprintf("frequency-%04d.msh", index))
        mv(path, destination)
        mv(joinpath(folder, "mesh-metadata.txt"),
            joinpath(directory, "logs", stem * "-metadata.txt"))
        rm(folder; recursive = true)
        push!(paths, destination)
    end
    return paths
end

"""
$(TYPEDSIGNATURES)

Validate the two-dimensional mesh at `path` for `model` in the caller's initialized Gmsh
session. Check the physical groups, terminal count, material coverage, surface edge
incidence, material and terminal ownership, and measurement lines with nodes outside the
field support.

For each material, compare the area from element vertices with its regions' area
\\[m²\\], within a relative tolerance of `1e-8` plus the chord allowance

```math
A_{\\mathrm{chord}} = \\frac{h^2}{6} \\sum_{\\mathrm{arcs}} |\\varphi|,
```

where ``h²`` is the longest boundary edge squared \\[m²\\]. The angle ``\\varphi`` is the
signed turn \\[rad\\] of each arc in the regions' contours. An arc of angle ``|\\varphi|``
cut into chords of length at most ``h`` loses at most ``|\\varphi|h²/12`` of area, so the
allowance is twice that loss. Higher-order elements use their Jacobians for the area
integral. Return `path` unchanged, and restore the caller's current Gmsh model.

# Errors

- Throws `ErrorException` with the failed check and mesh path.
"""
function validate(path::AbstractString, model::FieldModel)
    isfile(path) || error("mesh_path: mesh file does not exist: $path")
    current = gmsh.model.getCurrent()
    gmsh.model.add("LineCableModels-validation-$(time_ns())")
    try
        gmsh.merge(path)
        isempty(gmsh.model.getEntities(2)) && error("mesh_dimension: $path has no surfaces")
        isempty(gmsh.model.getEntities(3)) ||
            error("mesh_dimension: $path must be two-dimensional")
        expected = [
            (2, TAGS.air, "LCM/domain/air"),
            (2, TAGS.earth, "LCM/domain/earth"),
            (2, TAGS.air_pml, "LCM/domain/air_pml"),
            (2, TAGS.earth_pml, "LCM/domain/earth_pml"),
            (2, TAGS.pml, "LCM/domain/pml"),
            (1, TAGS.outer_boundary, "LCM/boundary/magnetic_dirichlet"),
            (1, TAGS.outer_air_boundary, "LCM/boundary/electric_reference_air"),
            (1, TAGS.outer_earth_boundary, "LCM/boundary/electric_reference_earth"),
            (1, TAGS.pml_inner_boundary, "LCM/boundary/pml_inner"),
            (1, TAGS.interface, "LCM/interface/air_earth"),
            (2, TAGS.conductors, "LCM/domain/conductors"),
            (2, TAGS.field_maps, "LCM/domain/field_maps")
        ]
        for (i, label) in enumerate(model.label)
            name = @sprintf("LCM/measurement_line/%04d", i)
            push!(expected, (1, TAGS.measurement_line_base + i, name))
            name = @sprintf("LCM/terminal/%04d/%s", i, label)
            push!(expected, (2, TAGS.terminal_base + i, name))
            name = @sprintf("LCM/terminal_contour/%04d", i)
            push!(expected, (1, TAGS.terminal_contour_base + i, name))
        end
        for material in model.materials
            push!(expected, (2, material.tag, material.name))
        end
        available = Set((dim, tag, gmsh.model.getPhysicalName(dim, tag))
        for (dim, tag) in gmsh.model.getPhysicalGroups())
        for group in expected
            group in available || error("physical_groups: $path is missing $group")
            dim, tag, name = group
            entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
            isempty(entities) && error("physical_groups: $path has empty group $name")
            for entity in entities
                element_tags = gmsh.model.mesh.getElements(dim, entity)[2]
                any(!isempty, element_tags) || error(
                    "physical_groups: $path has no elements on entity $entity in $name")
            end
        end
        terminal_count = count(gmsh.model.getPhysicalGroups(2)) do (dim, tag)
            return TAGS.terminal_base < tag <= TAGS.terminal_base + length(model.label)
        end
        terminal_count == length(model.label) ||
            error("terminal_count: $path has $terminal_count terminals")
        node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
        points = Dict(tag => (coordinates[3i - 2], coordinates[3i - 1])
        for (i, tag) in enumerate(node_tags))
        element_properties = Dict(kind => gmsh.model.mesh.getElementProperties(kind)
        for kind in first(gmsh.model.mesh.getElements()))
        curve_edges = Dict{Int, Vector{Tuple{UInt64, UInt64}}}()
        for (_, curve) in gmsh.model.getEntities(1)
            edges = Tuple{UInt64, UInt64}[]
            kinds, _, blocks = gmsh.model.mesh.getElements(1, curve)
            for (kind, vertices) in zip(kinds, blocks)
                _, _, _, n, _, primary = element_properties[kind]
                primary == 2 ||
                    error("material_partition: $path has unsupported boundary element $kind")
                for i in 1:n:length(vertices)
                    push!(edges, minmax(vertices[i], vertices[i + 1]))
                end
            end
            curve_edges[curve] = edges
        end
        areas = Dict{Int, Float64}()
        sizes = Dict{Int, Float64}()
        field_nodes = Set{UInt64}()
        for (_, surface) in gmsh.model.getEntities(2)
            incidence = Dict{Tuple{UInt64, UInt64}, Int}()
            surface_area = 0.0
            kinds, _, blocks = gmsh.model.mesh.getElements(2, surface)
            for (kind, vertices) in zip(kinds, blocks)
                union!(field_nodes, vertices)
                _, _, order, n, _, primary = element_properties[kind]
                primary in (3, 4) ||
                    error("material_partition: $path has unsupported surface element $kind")
                for i in 1:n:length(vertices)
                    for corner in 0:(primary - 1)
                        edge = minmax(vertices[i + corner], vertices[i + mod(corner + 1, primary)])
                        incidence[edge] = get(incidence, edge, 0) + 1
                    end
                    if order == 1
                        a = points[vertices[i]]
                        for corner in 1:(primary - 2)
                            b,
                            c = points[vertices[i + corner]],
                            points[vertices[i + corner + 1]]
                            surface_area += abs((b[1] - a[1]) * (c[2] - a[2]) -
                                                (b[2] - a[2]) * (c[1] - a[1])) / 2
                        end
                    end
                end
                if order > 1
                    local_points,
                    weights = gmsh.model.mesh.getIntegrationPoints(kind, "Gauss$(2order)")
                    _, determinants,
                    _ = gmsh.model.mesh.getJacobians(kind, local_points, surface)
                    surface_area += sum(abs(value) * weights[mod1(i, length(weights))]
                    for (i, value) in enumerate(determinants))
                end
            end
            surface_area > 0 ||
                error("material_partition: $path surface $surface has no positive-area elements")
            boundary = Set{Tuple{UInt64, UInt64}}()
            longest = 0.0
            for (_, curve) in gmsh.model.getBoundary([(2, surface)], false, false, false)
                edges = curve_edges[curve]
                isempty(edges) &&
                    error("material_partition: $path curve $curve has no boundary elements")
                for edge in edges
                    push!(boundary, edge)
                    get(incidence, edge, 0) == 1 || error(
                        "material_partition: $path surface $surface does not cover boundary curve $curve")
                    a, b = points[edge[1]], points[edge[2]]
                    longest = max(longest, (a[1] - b[1])^2 + (a[2] - b[2])^2)
                end
            end
            all(n == (edge in boundary ? 1 : 2) for (edge, n) in incidence) || error(
                "material_partition: $path surface $surface has missing, overlapping or nonconformal elements")
            areas[surface] = surface_area
            sizes[surface] = longest
        end
        assigned = Set{Int}()
        contours = [Contour(region.shape) for region in model.problem.system.regions]
        for (g, material) in enumerate(model.materials)
            surfaces = gmsh.model.getEntitiesForPhysicalGroup(2, material.tag)
            any(surface -> surface in assigned, surfaces) && error(
                "material_partition: $path has a surface with multiple material laws")
            union!(assigned, surfaces)
            indices = findall(==(g), model.material)
            expected_area = sum(area(model.problem.system.regions[r].shape)
            for r in indices)
            actual = sum(areas[surface] for surface in surfaces)
            h2 = maximum(sizes[surface] for surface in surfaces)
            total_turn = sum(
                abs(curve.span) for r in indices
                for contour in (contours[r],)
                for boundary in (contour.outer, contour.holes...)
                for curve in boundary if curve isa Arc;
                init = 0.0)
            allowance = total_turn * h2 / 6
            abs(actual - expected_area) <= 1e-8 * expected_area + allowance || error(
                "material_partition: $path material $(material.name) covers $actual m², expected $expected_area m²")
        end
        conductor_surfaces = Set{Int}()
        for material in model.materials
            material.conductor || continue
            union!(conductor_surfaces, gmsh.model.getEntitiesForPhysicalGroup(2, material.tag))
        end
        terminals = Set{Int}()
        for i in eachindex(model.label)
            surfaces = gmsh.model.getEntitiesForPhysicalGroup(2, TAGS.terminal_base + i)
            any(surface -> surface in terminals || surface ∉ conductor_surfaces, surfaces) &&
                error(
                    "material_partition: $path terminal $i has duplicate or nonconductive surface ownership")
            union!(terminals, surfaces)
        end
        terminals == conductor_surfaces ||
            error("material_partition: $path has a conductor surface with no terminal owner")
        for i in eachindex(model.label)
            for curve in
                gmsh.model.getEntitiesForPhysicalGroup(1, TAGS.measurement_line_base + i)
                for block in last(gmsh.model.mesh.getElements(1, curve))
                    isdisjoint(field_nodes, block) ||
                        error("measurement_line: $path line $i shares field nodes")
                end
            end
        end
    finally
        gmsh.model.remove()
        isempty(current) || gmsh.model.setCurrent(current)
    end
    return path
end

"""
$(TYPEDSIGNATURES)

Add the gauge tree, physical curve 9001, and voltage sampling nodes, physical point 9002,
to the mesh at `path` when those groups are absent. Process boundary elements before media
elements, in GetDP 3.5 element and edge order, with union by size and path compression.
Select the vertices of media elements whose coordinate ranges cross a measurement line.
The model's geometry and mesh coordinates are in \\[m\\].

Write the mesh in binary version 4.1 with `Mesh.SaveAll` to keep unlabelled curves.
Return `path` and restore the caller's current Gmsh model and mesh output options.
"""
function gauge!(path::AbstractString, model::FieldModel)
    current = gmsh.model.getCurrent()
    saved = Dict(name => gmsh.option.getNumber(name)
    for name in ("Mesh.Binary", "Mesh.SaveAll", "Mesh.MshFileVersion"))
    gmsh.model.add("LineCableModels-gauge-$(time_ns())")
    try
        gmsh.merge(path)
        tree_ready = (1, TAGS.gauge_tree) in gmsh.model.getPhysicalGroups(1)
        sampling_ready = (0, TAGS.sampling_nodes) in gmsh.model.getPhysicalGroups(0)
        tree_ready && sampling_ready && return path
        boundaries = Set([TAGS.outer_boundary;
                          TAGS.terminal_contour_base .+ eachindex(model.label)])
        media = Set([TAGS.air, TAGS.earth, TAGS.air_pml, TAGS.earth_pml])
        union!(media, (material.tag for material in model.materials if !material.conductor))
        _, all_tags, _ = gmsh.model.mesh.getElements()
        last_tag = maximum(maximum, all_tags)
        boundary = @NamedTuple{tag::Int, nodes::Vector{Int}}[]
        domain = @NamedTuple{tag::Int, nodes::Vector{Int}}[]
        for (dim, entity) in gmsh.model.getEntities()
            kinds, tags, nodes = gmsh.model.mesh.getElements(dim, entity)
            physical = gmsh.model.getPhysicalGroupsForEntity(dim, entity)
            for (p, tag) in enumerate(physical), block in eachindex(kinds)

                n_elements = length(tags[block])
                if (dim == 1 && tag in boundaries) || (dim == 2 && tag in media)
                    selected = dim == 1 ? boundary : domain
                    element_properties = gmsh.model.mesh.getElementProperties(kinds[block])
                    _, _, _, n, _, primary = element_properties
                    for i in 1:n_elements
                        id = p == 1 ? Int(tags[block][i]) : last_tag + i
                        connectivity = Int[nodes[block][(i - 1) * n + j] for j in 1:primary]
                        push!(selected, (tag = id, nodes = connectivity))
                    end
                end
                p > 1 && (last_tag += n_elements)
            end
        end
        sort!(boundary; by = cell -> cell.tag)
        sort!(domain; by = cell -> cell.tag)
        if !tree_ready
            largest = maximum(maximum(cell.nodes) for cells in (boundary, domain)
            for cell in cells)
            parents, sizes = collect(1:largest), ones(Int, largest)
            tree = NTuple{2, Int}[]
            for cells in (boundary, domain), cell in cells

                n = length(cell.nodes)
                for i in 1:n, j in (i + 1):n

                    (j == i + 1 || (i == 1 && j == n)) || continue
                    a, b = cell.nodes[i], cell.nodes[j]
                    ra, rb = a, b
                    while parents[ra] != ra
                        parents[ra] = parents[parents[ra]]
                        ra = parents[ra]
                    end
                    while parents[rb] != rb
                        parents[rb] = parents[parents[rb]]
                        rb = parents[rb]
                    end
                    ra == rb && continue
                    push!(tree, minmax(a, b))
                    sizes[ra] < sizes[rb] && ((ra, rb) = (rb, ra))
                    parents[rb] = ra
                    sizes[ra] += sizes[rb]
                end
            end
            curve = gmsh.model.addDiscreteEntity(1)
            edges = collect(Iterators.flatten(tree))
            gmsh.model.mesh.addElementsByType(curve, 1, Int[], edges)
            name = "LCM/gauge_tree/GetDP-3.5"
            gmsh.model.addPhysicalGroup(1, [curve], TAGS.gauge_tree, name)
        end
        if !sampling_ready
            tags, coordinates, _ = gmsh.model.mesh.getNodes()
            x, y = zeros(maximum(tags)), zeros(maximum(tags))
            for (i, tag) in enumerate(tags)
                x[tag], y[tag] = coordinates[3i - 2], coordinates[3i - 1]
            end
            paths = NTuple{3, Float64}[]
            for i in eachindex(model.label)
                tag = TAGS.measurement_line_base + i
                _, xyz = gmsh.model.mesh.getNodesForPhysicalGroup(1, tag)
                low = minimum(@view(xyz[2:3:end]))
                high = maximum(@view(xyz[2:3:end]))
                push!(paths, (xyz[1], low, high))
            end
            selected_nodes = Set{Int}()
            for cell in domain
                xmin = minimum(x[node] for node in cell.nodes)
                xmax = maximum(x[node] for node in cell.nodes)
                ymin = minimum(y[node] for node in cell.nodes)
                ymax = maximum(y[node] for node in cell.nodes)
                any(xmin <= px <= xmax && ymin <= hi && ymax >= lo
                for (px, lo, hi) in paths) || continue
                union!(selected_nodes, cell.nodes)
            end
            points = gmsh.model.addDiscreteEntity(0)
            sampling_vertices = sort!(collect(selected_nodes))
            gmsh.model.mesh.addElementsByType(points, 15, Int[], sampling_vertices)
            name = "LCM/voltage_sampling_nodes"
            gmsh.model.addPhysicalGroup(0, [points], TAGS.sampling_nodes, name)
        end
        gmsh.option.setNumber("Mesh.Binary", 1)
        gmsh.option.setNumber("Mesh.SaveAll", 1)
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.write(path)
    finally
        gmsh.model.remove()
        isempty(current) || gmsh.model.setCurrent(current)
        for (name, value) in saved
            gmsh.option.setNumber(name, value)
        end
    end
    return path
end
