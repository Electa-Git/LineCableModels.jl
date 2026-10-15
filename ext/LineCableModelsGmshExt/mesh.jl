const _GMSH_GEOMETRY_TOLERANCE = 1.0e-8
const _GMSH_BOOLEAN_TOLERANCE = 0.0

function _write_json_atomic(path::String, value)
    mkpath(dirname(path))
    temporary = tempname(dirname(path))
    open(temporary, "w") do io
        JSON3.pretty(io, value)
        write(io, '\n')
    end
    mv(temporary, path; force = true)
    return path
end

function _gmsh_version()
    try
        return gmsh.option.get_string("General.Version")
    catch
        return "unknown"
    end
end

function _mesh_fingerprint(model::FEMResolvedModel, gmsh_version::String,
        execution::ComputationOptions, frequency_index::Int=length(model.problem.frequencies);
        physical_geometry=nothing)
    evidence = (
        native_mesh_version=2,
        system=ImportExport.serialize_value(model.problem.system),
        regions=model.region_plans, cable_boundaries=model.cable_boundaries,
        cable_hosts=model.cable_hosts, terminal_ids=model.terminal_ids, tags=model.tags,
        frequency=model.problem.frequencies[frequency_index],
        gamma=model.prescribed_gamma[frequency_index],
        materials=[(m.mu_r,m.admittivity[frequency_index]) for m in model.material_plans],
        earth=(rho=ImportExport.serialize_value(model.earth_materials[frequency_index].rho),
            eps_r=model.earth_materials[frequency_index].eps_r,
            mu_r=model.earth_materials[frequency_index].mu_r),
        air=(sigma=inv(model.problem.earth_props.layers[1].rho),
            epsilon_r=model.problem.earth_props.layers[1].eps_r,
            mu_r=model.problem.earth_props.layers[1].mu_r),
        controls=(; (key=>getproperty(execution.data,key) for key in FEM_MESH_CONTROL_FIELDS)...),
        native_sources=(; (key=>bytes2hex(sha256(Vector{UInt8}(codeunits(getproperty(FEM_GETDP_SOURCES,key)))))
            for key in (:parameters,:geometry,:mesh))...),
        physical_geometry, gmsh_version,
        geometry_tolerance=_GMSH_GEOMETRY_TOLERANCE,
        boolean_tolerance=_GMSH_BOOLEAN_TOLERANCE)
    io=IOBuffer(); JSON3.write(io,evidence)
    return bytes2hex(sha256(take!(io)))
end

function _expected_physical_groups(model::FEMResolvedModel)
    groups = Tuple{Int, Int, String}[
        (
            2, model.tags.air, "LCM/domain/air"),
        (
            2, model.tags.earth, "LCM/domain/earth"),
        (
            2, model.tags.air_pml, "LCM/domain/air_pml"),
        (
            2, model.tags.earth_pml, "LCM/domain/earth_pml"),
        (
            2, model.tags.pml, "LCM/domain/pml"),
        (
            1,
            model.tags.outer_boundary,
            "LCM/boundary/magnetic_dirichlet"),
        (
            1,
            model.tags.outer_air_boundary,
            "LCM/boundary/electric_reference_air"),
        (
            1,
            model.tags.outer_earth_boundary,
            "LCM/boundary/electric_reference_earth"),
        (
            1,
            model.tags.pml_inner_boundary,
            "LCM/boundary/pml_inner"
        ),
        (
            1, model.tags.interface, "LCM/interface/air_earth"),
        (
            2, 6_001, "LCM/domain/conductors"),
        (
            2, 6_003, "LCM/domain/field_maps")
    ]
    for (index, name) in enumerate(model.terminal_names)
        push!(groups, (1, model.tags.measurement_line_base + index,
            @sprintf("LCM/measurement_line/%04d", index)))
        push!(groups, (2, model.tags.terminal_base + index, name))
        push!(groups,
            (
                1,
                model.tags.terminal_contour_base + index,
                @sprintf("LCM/terminal_contour/%04d", index)
            ))
    end
    for material in model.material_plans
        push!(groups, (2, material.physical_tag, material.physical_name))
    end
    return groups
end

function _inspect_loaded_mesh(model::FEMResolvedModel, mesh_path::String)
    isempty(gmsh.model.get_entities(2)) && _fem_error(
        :mesh,
        model.problem.system.system_id,
        :mesh_dimension,
        "mesh has no two-dimensional elements"
    )
    !isempty(gmsh.model.get_entities(3)) && _fem_error(
        :mesh,
        model.problem.system.system_id,
        :mesh_dimension,
        "mesh must be two-dimensional"
    )
    available = Set(
        (Int(dim), Int(tag), gmsh.model.get_physical_name(dim, tag))
    for (dim, tag) in gmsh.model.get_physical_groups()
    )
    for expected in _expected_physical_groups(model)
        expected in available || _fem_error(
            :mesh,
            model.problem.system.system_id,
            :physical_groups,
            "mesh $(mesh_path) is missing physical group $(expected)"
        )
        entities = gmsh.model.get_entities_for_physical_group(expected[1], expected[2])
        isempty(entities) && _fem_error(
            :mesh,
            model.problem.system.system_id,
            :physical_groups,
            "mesh $(mesh_path) has an empty physical group $(expected)"
        )
        for entity in entities
            _, element_tags, _ = gmsh.model.mesh.get_elements(expected[1], entity)
            any(!isempty, element_tags) || _fem_error(
                :mesh,
                model.problem.system.system_id,
                :physical_groups,
                "mesh $(mesh_path) has no elements on entity $(entity) " *
                "in physical group $(expected)"
            )
        end
    end
    terminal_groups = count(
        group -> begin
            dim, tag = group
            dim == 2 &&
                model.tags.terminal_base < tag <=
                model.tags.terminal_base + length(model.terminal_ids)
        end,
        gmsh.model.get_physical_groups())
    terminal_groups == length(model.terminal_ids) || _fem_error(
        :mesh,
        model.problem.system.system_id,
        :terminal_count,
        "mesh terminal count $terminal_groups differs from " *
        "$(length(model.terminal_ids))"
    )
    _inspect_material_coverage(model, mesh_path)
    # Sampling curves must have their own nodes; they are absent from the
    # field support and cannot alter the tree gauge or electrode contours.
    field_nodes = Set{UInt64}()
    for (_,surface) in gmsh.model.get_entities(2)
        for block in last(gmsh.model.mesh.get_elements(2,surface))
            union!(field_nodes,block)
        end
    end
    for i in eachindex(model.terminal_ids)
        for curve in gmsh.model.get_entities_for_physical_group(1,model.tags.measurement_line_base+i)
            for block in last(gmsh.model.mesh.get_elements(1,curve))
                isempty(intersect(field_nodes,block)) || _fem_error(
                    :mesh,model.problem.system.system_id,:measurement_line,
                    "measurement line $i shares field nodes; regenerate the mesh")
            end
        end
    end
    return _measurement_line_ratios(model)
end

# Read-only measurements of the emitted mesh, not a second sizing law.
function _measurement_line_ratios(model)
    coordinates = Dict{UInt64,Vector{Float64}}()
    xyz(node) = get!(coordinates,UInt64(node)) do
        first(gmsh.model.mesh.get_node(node))
    end
    metal = Set(surface for material in model.material_plans if material.kind === :conductor
        for surface in gmsh.model.get_entities_for_physical_group(2,material.physical_tag))
    maximums = Float64[]
    for terminal in eachindex(model.terminal_ids)
        maximum_ratio = 0.0
        for curve in gmsh.model.get_entities_for_physical_group(1,model.tags.measurement_line_base+terminal)
            kinds,_,blocks = gmsh.model.mesh.get_elements(1,curve)
            for (kind,nodes) in zip(kinds,blocks)
                @assert kind == 1
                for i in 1:2:length(nodes)
                    a,b = xyz(nodes[i]),xyz(nodes[i+1])
                    line_length = norm(b-a)
                    for t in (.2113248654051871,.5,.7886751345948129)
                        point = (1-t)*a+t*b
                        element,_,_,_,_ = gmsh.model.mesh.get_element_by_coordinates(point...,2)
                        _,vertices,_,surface = gmsh.model.mesh.get_element(element)
                        surface in metal && continue
                        longest = maximum(norm(xyz(vertices[j])-xyz(vertices[mod1(j+1,length(vertices))])) for j in eachindex(vertices))
                        maximum_ratio = max(maximum_ratio,line_length/longest)
                    end
                end
            end
        end
        push!(maximums,maximum_ratio)
    end
    return maximums
end

# A conservative chord-area allowance, proportional to the actual boundary
# segment length squared. Polygon polylines use the same absolute-turn bound
# as their native spacing rule; rectangles retain exact straight-line corners.
_mesh_area_allowance(shape::DataModel.Polygon, h2) =
    h2 * _polygon_total_turn(_normalize_polygon_vertices(_shape_points(shape))) / 8
_mesh_area_allowance(::DataModel.Rectangle, h2) = 0.0
_mesh_area_allowance(::DataModel.Disk, h2) = π * h2 / 3
_mesh_area_allowance(shape::DataModel.Annulus, h2) = (iszero(shape.ri) ? π : 2π) * h2 / 3
_mesh_area_allowance(shape::DataModel.Ellipse, h2) = π * (max(shape.a, shape.b) / min(shape.a, shape.b))^2 * h2 / 3
# Offset ellipses use the adapter's existing polygonal boundary, not native
# arcs. Its fixed tessellation error must not shrink when mesh edges subdivide.
_mesh_area_allowance(shape::DataModel.EllipseOffset, h2) =
    abs(LineCableModels.area(shape) - abs(_signed_area(_shape_points(shape))))
_mesh_area_allowance(::DataModel.SectorShape, h2) = 4π * h2 / 3
_mesh_area_allowance(::DataModel.BentStrip, h2) = 2π * h2 / 3
_mesh_area_allowance(shape::DataModel.ShellShape, h2) =
    _mesh_area_allowance(shape.outer, h2) + _mesh_area_allowance(shape.inner, h2)
_mesh_area_allowance(shape::DataModel.DifferenceShape, h2) =
    _mesh_area_allowance(shape.outer, h2) + sum(hole -> _mesh_area_allowance(hole, h2), shape.holes; init=0.0)
_mesh_area_allowance(shape::DataModel.AssemblyShape, h2) =
    sum(member -> _mesh_area_allowance(member, h2), shape.members; init=0.0)

function _inspect_material_coverage(model::FEMResolvedModel, mesh_path::String)
    fail(message) = _fem_error(:mesh, model.problem.system.system_id,
        :material_partition, "mesh $mesh_path: $message")
    node_tags, coordinates, _ = gmsh.model.mesh.get_nodes()
    points = Dict(UInt64(tag) => (coordinates[3i-2], coordinates[3i-1])
                  for (i, tag) in enumerate(node_tags))
    properties = Dict{Int, Tuple{Int, Int, Int}}()
    function element_properties(kind)
        get!(properties, kind) do
            _, _, order, count, _, primary = gmsh.model.mesh.get_element_properties(kind)
            (Int(order), Int(count), Int(primary))
        end
    end
    curve_edges = Dict{Int, Vector{Tuple{UInt64, UInt64}}}()
    function boundary_edges(curve)
        get!(curve_edges, curve) do
            edges = Tuple{UInt64, UInt64}[]
            kinds, _, blocks = gmsh.model.mesh.get_elements(1, curve)
            for (kind, vertices) in zip(kinds, blocks)
                _, count, primary = element_properties(kind)
                primary == 2 || fail("unsupported boundary element $kind")
                for index in 1:count:length(vertices)
                    push!(edges, minmax(vertices[index], vertices[index+1]))
                end
            end
            isempty(edges) && fail("curve $curve has no boundary elements")
            edges
        end
    end
    areas = Dict{Int, Float64}()
    sizes = Dict{Int, Float64}()
    for (_, surface) in gmsh.model.get_entities(2)
        incidence = Dict{Tuple{UInt64, UInt64}, Int}()
        surface_area = 0.0
        kinds, _, blocks = gmsh.model.mesh.get_elements(2, surface)
        for (kind, vertices) in zip(kinds, blocks)
            order, count, primary = element_properties(kind)
            primary in (3, 4) || fail("unsupported surface element $kind")
            for index in 1:count:length(vertices)
                for corner in 0:primary-1
                    edge = minmax(vertices[index+corner], vertices[index+mod(corner+1, primary)])
                    incidence[edge] = get(incidence, edge, 0) + 1
                end
                if order == 1
                    a = points[vertices[index]]
                    for corner in 1:primary-2
                        b, c = points[vertices[index+corner]], points[vertices[index+corner+1]]
                        surface_area += abs((b[1]-a[1])*(c[2]-a[2]) -
                                            (b[2]-a[2])*(c[1]-a[1])) / 2
                    end
                end
            end
            if order > 1
                local_points, weights = gmsh.model.mesh.get_integration_points(kind, "Gauss$(2order)")
                _, determinants, _ = gmsh.model.mesh.get_jacobians(kind, local_points, surface)
                surface_area += sum(abs(value) * weights[mod1(index, length(weights))]
                                    for (index, value) in enumerate(determinants))
            end
        end
        surface_area > 0 || fail("surface $surface has no positive-area elements")
        boundary = Set{Tuple{UInt64, UInt64}}()
        maximum_edge_squared = 0.0
        for (_, curve) in gmsh.model.get_boundary([(2, surface)], false, false, false)
            for edge in boundary_edges(Int(curve))
                push!(boundary, edge)
                get(incidence, edge, 0) == 1 || fail(
                    "surface $surface does not cover boundary curve $curve")
                a, b = points[edge[1]], points[edge[2]]
                maximum_edge_squared = max(maximum_edge_squared,
                    (a[1]-b[1])^2 + (a[2]-b[2])^2)
            end
        end
        all(count == (edge in boundary ? 1 : 2) for (edge, count) in incidence) ||
            fail("surface $surface has missing, overlapping or nonconformal elements")
        areas[surface] = surface_area
        sizes[surface] = maximum_edge_squared
    end
    assigned = Set{Int}()
    for (index, material) in enumerate(model.material_plans)
        surfaces = Int.(gmsh.model.get_entities_for_physical_group(2, material.physical_tag))
        any(surface -> surface in assigned, surfaces) && fail("a surface owns multiple material laws")
        union!(assigned, surfaces)
        regions = filter(region -> region.material_index == index, model.region_plans)
        expected = sum(region -> Float64(LineCableModels.area(region.shape)), regions)
        actual = sum(surface -> areas[surface], surfaces)
        h2 = maximum(surface -> sizes[surface], surfaces)
        chord_allowance = sum(region -> _mesh_area_allowance(region.shape, h2), regions)
        abs(actual - expected) <= 1e-8 * expected + chord_allowance || fail(
            "material $(material.physical_name) covers $actual m², expected $expected m²")
    end
    conductor_surfaces = Set{Int}()
    for material in model.material_plans
        material.kind === :conductor || continue
        union!(conductor_surfaces, gmsh.model.get_entities_for_physical_group(2, material.physical_tag))
    end
    terminals = Set{Int}()
    for index in eachindex(model.terminal_ids)
        surfaces = Int.(gmsh.model.get_entities_for_physical_group(2, model.tags.terminal_base + index))
        any(surface -> surface in terminals || surface ∉ conductor_surfaces, surfaces) &&
            fail("terminal $index has duplicate or nonconductive surface ownership")
        union!(terminals, surfaces)
    end
    terminals == conductor_surfaces || fail("a conductor surface has no terminal owner")
    return nothing
end

function _inspect_mesh_file(mesh_path::String, model::FEMResolvedModel)
    isfile(mesh_path) || _fem_error(
        :mesh,
        model.problem.system.system_id,
        :mesh_path,
        "mesh file does not exist: $mesh_path"
    )
    current = try
        gmsh.model.get_current()
    catch
        ""
    end
    validation_name = "LineCableModelsFEM-validation-$(time_ns())"
    gmsh.model.add(validation_name)
    try
        gmsh.merge(mesh_path)
        return _inspect_loaded_mesh(model, mesh_path)
    catch exception
        exception isa LineCableModelsFEMError && rethrow()
        _fem_error(
            :mesh,
            model.problem.system.system_id,
            :mesh_path,
            "failed to read mesh $mesh_path: $(sprint(showerror, exception))"
        )
    finally
        try
            gmsh.model.remove()
        catch
        end
        isempty(current) || try
            gmsh.model.set_current(current)
        catch
        end
    end
    return nothing
end

function validate(mesh_path::String, model::FEMResolvedModel)
    _inspect_mesh_file(mesh_path, model)
    return mesh_path
end

function _read_native_mesh_values(path)
    values=Dict{String,Any}()
    repeated=("region_mesh_sizes_m","cable_outer_mesh_sizes_m","cable_interface_mesh_sizes_m","conductor_targets")
    for line in eachline(path)
        fields=split(line); isempty(fields) && continue
        name=String(first(fields))
        value=if name == "conductor_targets"
            fields[2] == "null" ? nothing : NamedTuple{(:delta,:first_size,:extent,:bulk)}(Tuple(parse.(Float64,fields[2:end])))
        else
            numbers=parse.(Float64,fields[2:end])
            length(numbers) == 1 ? only(numbers) : numbers
        end
        if name in repeated
            push!(get!(Vector{Any},values,name),value)
        else
            haskey(values,name) && error("duplicate native mesh value $name")
            values[name]=value
        end
    end
    return values
end

function _mesh_metadata(model, fingerprint, gmsh_version, source, execution,
        frequency_index, native_values)
    metadata = Dict{String,Any}(String(key)=>value for (key,value) in pairs(native_values))
    merge!(metadata,Dict(
        "schema"=>"LineCableModels.FEMMesh", "version"=>3, "fingerprint"=>fingerprint,
        "gmsh_version"=>gmsh_version, "source"=>String(source), "mesh_dimension"=>2,
        "terminal_count"=>length(model.terminal_ids), "terminal_ids"=>model.terminal_ids,
        "physical_groups"=>[(dimension=dim,tag,name) for (dim,tag,name) in _expected_physical_groups(model)],
        "frequency_index"=>frequency_index, "frequency_hz"=>model.problem.frequencies[frequency_index],
        "earth_inputs"=>ImportExport.serialize_value((rho=model.earth_materials[frequency_index].rho,
            eps_r=model.earth_materials[frequency_index].eps_r, mu_r=model.earth_materials[frequency_index].mu_r)),
        "gamma"=>ImportExport.serialize_value(model.prescribed_gamma[frequency_index]),
        "pml_element_family"=>execution.data.pml_element_family == 1 ? "quadrangle" : "triangle",
        "pml_layers"=>execution.data.pml_layers, "pml_grading"=>execution.data.pml_grading))
    return metadata
end

function _copy_or_link_mesh(source::String, destination::String)
    try
        Base.Filesystem.hardlink(source, destination)
    catch exception
        @debug "Falling back to a physical FEM mesh copy" source destination exception
        cp(source, destination; force = true)
    end
    return destination
end

function _copy_mesh_snapshot!(
        source::String,
        run_mesh::String,
        metadata_path::String,
        metadata
)
    mkpath(dirname(run_mesh))
    if !isfile(run_mesh) || !Base.Filesystem.samefile(source,run_mesh)
        temporary = tempname(dirname(run_mesh))
        try
            _copy_or_link_mesh(source, temporary)
            mv(temporary,run_mesh;force=true)
        finally
            rm(temporary;force=true)
        end
    end
    _write_json_atomic(metadata_path, metadata)
    return run_mesh
end

function _resolve_mesh(run::FEMRun,model::FEMResolvedModel,
        execution::ComputationOptions,frequency_index::Int)
    gmsh_version = _gmsh_version()
    cad = joinpath(run.path,"input","physical.geo")
    fingerprint = _mesh_fingerprint(model,gmsh_version,execution,frequency_index;
        physical_geometry=bytes2hex(open(sha256,cad)))
    reference_mesh = frequency_index == length(model.problem.frequencies)
    stem = @sprintf("frequency_%04d",frequency_index)
    run_mesh = joinpath(run.path,"mesh",stem*".msh")
    run_metadata = joinpath(run.path,"mesh",stem*".json")
    mkpath(dirname(run_mesh))
    selected=nothing; source=:generated; native_values=Dict{String,Any}()
    if execution.data.mesh_policy === :reuse
        for (candidate,record,origin) in ((run_mesh,run_metadata,:resume),)
            isfile(candidate) && isfile(record) || continue
            try
                metadata=JSON3.read(read(record,String))
                String(metadata.fingerprint)==fingerprint || error("fingerprint mismatch")
                validate(candidate,model)
                selected=candidate; source=origin; native_values=metadata
                break
            catch exception
                @warn "Ignoring an invalid retained FEM mesh" candidate exception
            end
        end
        if selected === nothing && execution.data.mesh_path !== nothing
            path=abspath(execution.data.mesh_path)
            if isdir(path)
                for record in readdir(joinpath(path,"mesh");join=true)
                    occursin(r"^frequency_\d{4,}\.json$",basename(record)) || continue
                    candidate=replace(record,r"\.json$"=>".msh")
                    try
                        metadata=JSON3.read(read(record,String))
                        String(metadata.fingerprint)==fingerprint || continue
                        validate(candidate,model)
                        selected=candidate; source=:previous_run
                        native_values=Dict{String,Any}(String(key)=>value for (key,value) in pairs(metadata))
                        native_values["source_mesh"]=candidate
                        break
                    catch exception
                        @warn "Ignoring an invalid retained FEM mesh" candidate exception
                    end
                end
            elseif reference_mesh
                selected=path; source=:explicit
                native_values["measurement_line_max_size_ratios"]=_inspect_mesh_file(selected,model)
            end
        end
    end
    return (; frequency_index, gmsh_version, fingerprint, reference_mesh, run_mesh, run_metadata, stem, selected, source, native_values)
end

function _launch_mesh(run, plan, execution, temporary_mesh, temporary_native)
    directory=dirname(temporary_mesh)
    for path in _getdp_assets(joinpath(run.path,"input","getdp"))
        cp(path,joinpath(directory,basename(path)))
    end
    for name in ("model_data.pro","physical.geo")
        cp(joinpath(run.path,"input",name),joinpath(directory,name))
    end
    _write_native_mesh_entry(joinpath(directory,"model.geo"),"model_data.pro","physical.geo",".")
    log_path=joinpath(run.path,"logs",plan.stem*"-gmsh.log")
    command=Cmd(`$(Gmsh.gmsh_jll.gmsh()) model.geo -setnumber FrequencyIndex $(plan.frequency_index) -setstring MeshMetadataPath $(basename(temporary_native)) -2 -o $(basename(temporary_mesh)) -v $(max(1,execution.data.gmsh_verbosity))`;dir=directory)
    for (index, argument) in enumerate(command)
        name = index > 2 && command[index-2] == "-setstring" ? command[index-1] :
            "Gmsh argument $index"
        _check_native_path(name,argument;run_directory=run.path,solver="Gmsh",category=:mesh)
    end
    return open(log_path,"w") do log
        println(log,JSON3.write((;directory,command=collect(command))))
        flush(log)
        Base.run(pipeline(ignorestatus(command);stdout=log,stderr=log);wait=false)
    end
end

function _mesh_launch_error(run, model, plan)
    log_path=joinpath(run.path,"logs",plan.stem*"-gmsh.log")
    _fem_error(:mesh,model.problem.system.system_id,:native_geometry,
        "native Gmsh failed; see $log_path";run_directory=run.path)
end

const FEM_GAUGE_TREE_TAG = 9001
const FEM_GAUGE_TREE_NAME = "LCM/gauge_tree/GetDP-3.5"
const FEM_VOLTAGE_SAMPLING_NODE_TAG = 9002

# GetDP 3.5 GeoEntity.h edge order, including its non-Gmsh triangle ordering.
# Preserve GetDP 3.5 GeoData.cpp's sorted element order, including the tags
# assigned to secondary physical groups. Union by size/path compression changes
# only component bookkeeping, not the selected tree or the boundary forest.
_getdp_gauge_cell_edges(t) = t == 1 ? ((1,2),) :
    t == 2 ? ((1,2),(1,3),(2,3)) :
    t == 3 ? ((1,2),(1,4),(2,3),(3,4)) :
    throw(ArgumentError("unsupported GetDP gauge element type $t"))

function _getdp_gauge_elements(boundaries, media)
    _, all_tags, _ = gmsh.model.mesh.get_elements()
    last_tag = maximum(maximum, all_tags)
    record = Tuple{Int,Int,NTuple{4,Int}}
    boundary, domain = record[], record[]
    for (dim, entity) in gmsh.model.get_entities()
        types, tags, nodes = gmsh.model.mesh.get_elements(dim, entity)
        physical = gmsh.model.get_physical_groups_for_entity(dim, entity)
        for (p, tag) in enumerate(physical), block in eachindex(types)
            count = length(tags[block])
            selected = dim == 1 && tag in boundaries ? boundary :
                       dim == 2 && tag in media ? domain : nothing
            if selected !== nothing
                n = Int(gmsh.model.mesh.get_element_properties(types[block])[4])
                for i in 1:count
                    id = p == 1 ? Int(tags[block][i]) : last_tag + i
                    connectivity = ntuple(k -> k <= n ? Int(nodes[block][(i-1)*n+k]) : 0, 4)
                    push!(selected, (id, Int(types[block]), connectivity))
                end
            end
            p > 1 && (last_tag += count)
        end
    end
    sort!(boundary; by=first)
    sort!(domain; by=first)
    return boundary, domain
end

function _gauge_root!(parent, a)
    while parent[a] != a
        parent[a] = parent[parent[a]]
        a = parent[a]
    end
    return a
end

function _getdp_gauge_tree(boundary, domain)
    largest = maximum(maximum(last(cell)) for cells in (boundary,domain) for cell in cells)
    parent, sizes = collect(1:largest), ones(Int,largest)
    tree = NTuple{2,Int}[]
    for cells in (boundary,domain), (_,t,nodes) in cells, (i,j) in _getdp_gauge_cell_edges(t)
        a, b = nodes[i], nodes[j]
        ra, rb = _gauge_root!(parent,a), _gauge_root!(parent,b)
        ra == rb && continue
        push!(tree, minmax(a,b))
        sizes[ra] < sizes[rb] && ((ra,rb) = (rb,ra))
        parent[rb] = ra
        sizes[ra] += sizes[rb]
    end
    return tree
end

# Keep original element interpolation: select vertices of every media cell
# whose bounding box crosses a voltage path. GetDP expands these vertices to
# a one-ring band, including neighbouring cells at interfaces and endpoints.
function _measurement_field_nodes(model, domain)
    tags, coordinates, _ = gmsh.model.mesh.get_nodes()
    x, y = zeros(maximum(tags)), zeros(maximum(tags))
    for (i,tag) in pairs(tags)
        x[tag], y[tag] = coordinates[3i-2], coordinates[3i-1]
    end
    paths = map(eachindex(model.terminal_ids)) do terminal
        _, coordinates = gmsh.model.mesh.get_nodes_for_physical_group(
            1,model.tags.measurement_line_base+terminal)
        (coordinates[1],minimum(@view(coordinates[2:3:end])),
            maximum(@view(coordinates[2:3:end])))
    end
    selected = Set{Int}()
    for (_,type,nodes) in domain
        n = type == 2 ? 3 : 4
        xmin = xmax = x[nodes[1]]
        ymin = ymax = y[nodes[1]]
        for i in 2:n
            xmin, xmax = min(xmin,x[nodes[i]]), max(xmax,x[nodes[i]])
            ymin, ymax = min(ymin,y[nodes[i]]), max(ymax,y[nodes[i]])
        end
        any(paths) do (px,lo,hi)
            xmin <= px <= xmax && ymin <= hi && ymax >= lo
        end || continue
        for i in 1:n
            push!(selected,nodes[i])
        end
    end
    return sort!(collect(selected))
end

function _prepare_gauge_mesh!(mesh, model)
    current = gmsh.model.get_current()
    saved = Dict(name=>gmsh.option.get_number(name)
        for name in ("Mesh.Binary","Mesh.SaveAll","Mesh.MshFileVersion"))
    temporary = tempname(dirname(mesh))*".msh"
    gmsh.model.add("LineCableModelsFEM-gauge-$(time_ns())")
    try
        gmsh.merge(mesh)
        groups = gmsh.model.get_physical_groups(1)
        tree_ready = (1,FEM_GAUGE_TREE_TAG) in groups &&
            gmsh.model.get_physical_name(1,FEM_GAUGE_TREE_TAG)==FEM_GAUGE_TREE_NAME
        sampling_ready = (0,FEM_VOLTAGE_SAMPLING_NODE_TAG) in gmsh.model.get_physical_groups(0)
        tree_ready && sampling_ready && return mesh
        boundaries = Set([model.tags.outer_boundary;
            model.tags.terminal_contour_base .+ eachindex(model.terminal_ids)])
        media = Set([model.tags.air,model.tags.earth,model.tags.air_pml,model.tags.earth_pml,
            [m.physical_tag for m in model.material_plans if m.kind !== :conductor]...])
        boundary,domain = _getdp_gauge_elements(boundaries,media)
        if !tree_ready
            tree = _getdp_gauge_tree(boundary,domain)
            curve = gmsh.model.add_discrete_entity(1)
            gmsh.model.mesh.add_elements_by_type(curve,1,Int[],collect(Iterators.flatten(tree)))
            gmsh.model.add_physical_group(1,[curve],FEM_GAUGE_TREE_TAG,FEM_GAUGE_TREE_NAME)
        end
        if !sampling_ready
            nodes = _measurement_field_nodes(model,domain)
            points = gmsh.model.add_discrete_entity(0)
            gmsh.model.mesh.add_elements_by_type(points,15,Int[],nodes)
            gmsh.model.add_physical_group(0,[points],FEM_VOLTAGE_SAMPLING_NODE_TAG,
                "LCM/voltage_sampling_nodes")
        end
        gmsh.option.set_number("Mesh.Binary",1)
        # Retain unlabelled internal curves used by mesh coverage validation.
        gmsh.option.set_number("Mesh.SaveAll",1)
        gmsh.option.set_number("Mesh.MshFileVersion",4.1)
        gmsh.write(temporary)
        # Reused meshes can be hard links: replace rather than overwrite.
        mv(temporary,mesh;force=true)
    finally
        rm(temporary;force=true)
        gmsh.model.remove()
        isempty(current) || gmsh.model.set_current(current)
        for (name,value) in saved
            gmsh.option.set_number(name,value)
        end
    end
    return mesh
end

function _publish_mesh!(run, model, execution, plan; temporary_mesh=nothing,
        temporary_native=nothing, owner=nothing)
    (;frequency_index,gmsh_version,fingerprint,reference_mesh,
        run_mesh,run_metadata,selected,source,native_values)=plan
    if owner !== nothing
        selected=owner.run_mesh
        source=:shared
        native_values=JSON3.read(read(owner.run_metadata,String))
    end
    if selected === nothing
        ratios=_inspect_mesh_file(temporary_mesh,model)
        native_values=_read_native_mesh_values(temporary_native)
        native_values["measurement_line_max_size_ratios"]=ratios
        metadata=_mesh_metadata(model,fingerprint,gmsh_version,:generated,execution,frequency_index,native_values)
        _write_json_atomic(run_metadata,metadata)
        mv(temporary_mesh,run_mesh;force=true)
    elseif selected != run_mesh
        metadata=_mesh_metadata(model,fingerprint,gmsh_version,source,execution,frequency_index,native_values)
        _copy_mesh_snapshot!(selected,run_mesh,run_metadata,metadata)
    end
    if selected == run_mesh
        metadata=_mesh_metadata(model,fingerprint,gmsh_version,source,execution,frequency_index,native_values)
        _write_json_atomic(run_metadata,metadata)
    end
    if reference_mesh
        run.mesh_fingerprint=fingerprint
        run.mesh_source=source
    end
    _prepare_gauge_mesh!(run_mesh,model)
    return run_mesh
end

function _stop_mesh_processes!(processes)
    for process in processes
        process === nothing || !process_running(process) || kill(process)
    end
    for process in processes
        process === nothing && continue
        if timedwait(() -> process_exited(process),2.;pollint=.02) === :timed_out
            kill(process,9)
        end
        wait(process)
    end
end

function _select_meshes!(run::FEMRun,model::FEMResolvedModel,
        execution::ComputationOptions)
    # Resolve shared-session Gmsh operations before any external process starts.
    plans=[_resolve_mesh(run,model,execution,index)
        for index in eachindex(model.problem.frequencies)]
    owners=Dict{String,Int}(); owner_indices=zeros(Int,length(plans)); pending=Int[]
    # A retained sibling can supply a missing mesh during resume.
    for (index,plan) in pairs(plans)
        plan.selected==plan.run_mesh && get!(owners,plan.fingerprint,index)
    end
    for (index,plan) in pairs(plans)
        plan.selected === nothing || continue
        owner_indices[index]=get!(owners,plan.fingerprint,index)
        owner_indices[index]==index && push!(pending,index)
    end
    work=joinpath(run.path,"work"); mkpath(work)
    mesh_directories=[mktempdir(work;prefix="mesh-") for plan in plans]
    temporary_meshes=[joinpath(directory,"model.msh") for directory in mesh_directories]
    temporary_native=[joinpath(directory,"mesh-metadata.txt") for directory in mesh_directories]
    processes=Vector{Union{Nothing,Base.Process}}(nothing,length(plans))
    launch_error=Ref{Any}(nothing); aborted=Ref(false)
    try
        asyncmap(pending;ntasks=min(execution.data.frequency_workers,length(pending))) do index
            aborted[] && return
            try
                process=_launch_mesh(run,plans[index],execution,
                    temporary_meshes[index],temporary_native[index])
                processes[index]=process
                aborted[] && process_running(process) && kill(process)
                wait(process)
                if !success(process) && !aborted[]
                    _mesh_launch_error(run,model,plans[index])
                end
            catch exception
                if !aborted[]
                    launch_error[]=(index,exception)
                    aborted[]=true
                    _stop_mesh_processes!(processes)
                end
            end
        end
        if launch_error[] !== nothing
            index,exception=launch_error[]
            # Normal nonzero exits are failures; killed siblings are not.
            first_failed=findfirst(process -> process !== nothing &&
                !success(process) && process.termsignal==0,processes)
            if first_failed !== nothing && first_failed<index
                _mesh_launch_error(run,model,plans[first_failed])
            end
            throw(exception)
        end
        # Publication and validation use the caller's Gmsh session sequentially.
        return [ _publish_mesh!(run,model,execution,plan;
            temporary_mesh=temporary_meshes[index],temporary_native=temporary_native[index],
            owner=owner_indices[index] in (0,index) ? nothing : plans[owner_indices[index]])
            for (index,plan) in pairs(plans)]
    finally
        aborted[]=true
        _stop_mesh_processes!(processes)
        foreach(path -> rm(path;recursive=true,force=true),mesh_directories)
    end
end
