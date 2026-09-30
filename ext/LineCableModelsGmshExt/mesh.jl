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

function _mesh_fingerprint(
        model::FEMResolvedModel,
        gmsh_version::String,
    mesh_plan::FEMMeshPlan = last(model.mesh_plans)
)
    evidence = (
        material_partition_version = 18,
        conductor_mesh = model.conductor_mesh,
        conductor_materials = [(material.mu_r, material.admittivity[mesh_plan.frequency_index])
            for material in model.material_plans if material.kind === :conductor],
        system = ImportExport.serialize_value(model.problem.system),
        temperature = ImportExport.serialize_value(model.problem.temperature),
        earth_props = ImportExport.serialize_value(model.problem.earth_props),
        # Declarations alone cannot identify geometry produced by a changed
        # resolver. Fingerprint the exact domains and ownership sent to Gmsh.
        regions = model.region_plans,
        cable_boundaries = model.cable_boundaries,
        cable_hosts = model.cable_hosts,
        centre = model.centre,
        terminal_ids = model.terminal_ids,
        material_tags = getproperty.(model.material_plans, :physical_tag),
        material_names = getproperty.(model.material_plans, :physical_name),
        tags = model.tags,
        region_mesh_sizes = getproperty.(model.region_plans, :mesh_size),
        cable_outer_mesh_sizes = model.cable_outer_mesh_sizes,
        mesh_growth_factor = model.mesh_growth_factor,
        interface_refinement_factor = model.interface_refinement_factor,
        mesh_frequency = mesh_plan.frequency,
        domain_halfwidth = mesh_plan.domain_halfwidth,
        pml_thickness = mesh_plan.pml_thickness,
        pml_layers = mesh_plan.pml_layers,
        pml_grading = mesh_plan.pml_grading,
        pml_strips = mesh_plan.pml_strips,
        pml_element_family = mesh_plan.pml_element_family,
        domain_mesh_size = mesh_plan.domain_mesh_size,
        exterior_mesh_sizes = mesh_plan.exterior_mesh_sizes,
        exterior_start_radius = mesh_plan.exterior_start_radius,
        interface_mesh_size = mesh_plan.interface_mesh_size,
        cable_interface_mesh_sizes = mesh_plan.cable_interface_mesh_sizes,
        wave_mesh_sizes = mesh_plan.wave_mesh_sizes,
        wave_size_limits = mesh_plan.wave_size_limits,
        wave_decay_radii = mesh_plan.wave_decay_radii,
        gmsh_version,
        geometry_tolerance = _GMSH_GEOMETRY_TOLERANCE,
        boolean_tolerance = _GMSH_BOOLEAN_TOLERANCE
    )
    # Hash owned bytes: SHA's string/CodeUnits path can repeatedly hash the
    # entire immutable string during copyto! alias checks on Julia 1.12.
    io = IOBuffer()
    JSON3.write(io, evidence)
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
        push!(groups, (1, model.tags.voltage_path_base + index,
            @sprintf("LCM/voltage_path/%04d", index)))
        push!(groups, (0, model.tags.voltage_reference_base + index,
            @sprintf("LCM/voltage_reference/%04d", index)))
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
    # A detached line would introduce unrelated BF_Edge coefficients instead
    # of the trace of the solved field. Reject legacy and supplied meshes with
    # such paths before GetDP assembles the system.
    metal = Set(s for material in model.material_plans if material.kind === :conductor
        for s in gmsh.model.get_entities_for_physical_group(2, material.physical_tag))
    field_edges = Set{Tuple{UInt64, UInt64}}()
    for (_, surface) in gmsh.model.get_entities(2)
        surface in metal && continue
        for kind in first(gmsh.model.mesh.get_elements(2, surface))
            nodes = gmsh.model.mesh.get_element_edge_nodes(kind, surface, true)
            for k in 1:2:length(nodes)
                push!(field_edges, minmax(nodes[k], nodes[k+1]))
            end
        end
    end
    for i in eachindex(model.terminal_ids)
        for curve in gmsh.model.get_entities_for_physical_group(1, model.tags.voltage_path_base+i)
            kinds, _, blocks = gmsh.model.mesh.get_elements(1, curve)
            for (kind, nodes) in zip(kinds, blocks)
                count = gmsh.model.mesh.get_element_properties(kind)[4]
                for k in 1:count:length(nodes)
                    minmax(nodes[k], nodes[k+1]) in field_edges || _fem_error(
                        :mesh, model.problem.system.system_id, :voltage_path,
                        "voltage path $i is not conforming to the electric field mesh; regenerate the mesh")
                end
            end
        end
    end
    return nothing
end

# A conservative chord-area allowance, proportional to the actual boundary
# segment length squared. Polygon boundaries have no discretization allowance.
_mesh_area_allowance(::Union{DataModel.Polygon, DataModel.Rectangle}, h2) = 0.0
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

function _validate_mesh_file(model::FEMResolvedModel, mesh_path::String)
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
        _inspect_loaded_mesh(model, mesh_path)
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

_interface_footprint(design, position) =
    (position.x, abs(position.y) + LineCableModels.outer_radius(design))

function _configure_mesh!(
        model::FEMResolvedModel,
        geometry::FEMGeometry,
        mesh_plan::FEMMeshPlan = last(model.mesh_plans)
)
    gmsh.option.set_number("Mesh.MshFileVersion", 4.1)
    gmsh.option.set_number("Mesh.Binary", 1)
    # Preserve boundary elements even on internal same-material seams. MSH 4
    # retains physical groups with SaveAll, enabling coverage checks on reload.
    gmsh.option.set_number("Mesh.SaveAll", 1)
    gmsh.option.set_number("Mesh.MeshSizeMin",
        Float64(min(model.fine_mesh_size, minimum(mesh_plan.wave_mesh_sizes))))
    exterior_size = max(mesh_plan.domain_mesh_size, maximum(mesh_plan.exterior_mesh_sizes))
    gmsh.option.set_number("Mesh.MeshSizeMax", Float64(exterior_size))
    gmsh.option.set_number("Mesh.MeshSizeFromPoints", 1)
    gmsh.option.set_number("Mesh.MeshSizeExtendFromBoundary", 0)
    transition_fields = Int[]
    for cable_index in eachindex(geometry.cable_curves)
        cable_curves = sort!(unique(geometry.cable_curves[cable_index]))
        isempty(cable_curves) && continue
        distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.set_numbers(distance, "CurvesList", cable_curves)
        gmsh.model.mesh.field.set_number(distance, "Sampling", 100)
        threshold = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.set_number(threshold, "InField", distance)
        cable_size = model.cable_outer_mesh_sizes[cable_index]
        gmsh.model.mesh.field.set_number(
            threshold, "SizeMin", Float64(cable_size)
        )
        gmsh.model.mesh.field.set_number(
            threshold, "SizeMax", Float64(exterior_size)
        )
        gmsh.model.mesh.field.set_number(threshold, "DistMin", 0.0)
        transition_distance = max(
            cable_size,
            (exterior_size - cable_size) /
            max(model.mesh_growth_factor - one(model.mesh_growth_factor), 1e-12)
        )
        gmsh.model.mesh.field.set_number(
            threshold, "DistMax", Float64(transition_distance)
        )
        push!(transition_fields, threshold)
    end
    if exterior_size > mesh_plan.domain_mesh_size
        # The previous global cap also bounded cable interiors. Retain it on
        # every cable material, including cores whose local target is larger.
        constant = gmsh.model.mesh.field.add("MathEval")
        gmsh.model.mesh.field.set_string(constant, "F", string(mesh_plan.domain_mesh_size))
        restricted = gmsh.model.mesh.field.add("Restrict")
        gmsh.model.mesh.field.set_number(restricted, "InField", constant)
        gmsh.model.mesh.field.set_numbers(restricted, "SurfacesList",
            reduce(vcat, geometry.material_surfaces; init=Int[]))
        gmsh.model.mesh.field.set_number(restricted, "IncludeBoundary", 1)
        push!(transition_fields, restricted)
    end
    # Retain the central bulk cap, then grow gradually into the remote buffer.
    # This does not multiply conductor sizes or change the local growth slope.
    for (surfaces, size) in zip((geometry.air_surfaces, geometry.earth_surfaces),
            mesh_plan.exterior_mesh_sizes)
        size == mesh_plan.domain_mesh_size && exterior_size == size && continue
        radial = gmsh.model.mesh.field.add("MathEval")
        cx = model.centre[1]
        expression = "Min($(size), $(mesh_plan.domain_mesh_size) + " *
            "$(model.mesh_growth_factor-1) * Max(0, " *
            "Sqrt((x-($(cx)))^2+y^2)-$(mesh_plan.exterior_start_radius)))"
        gmsh.model.mesh.field.set_string(radial, "F", expression)
        restricted = gmsh.model.mesh.field.add("Restrict")
        gmsh.model.mesh.field.set_number(restricted, "InField", radial)
        gmsh.model.mesh.field.set_numbers(restricted, "SurfacesList", surfaces)
        gmsh.model.mesh.field.set_number(restricted, "IncludeBoundary", 1)
        push!(transition_fields, restricted)
    end
    # Localize the interface forcing to projected cable footprints in each
    # medium. Keep the prescribed wave sizes, decay distances and remote caps.
    for (surfaces, wave_size, decay_radius, remote_size) in zip(
        (geometry.air_surfaces, geometry.earth_surfaces),
        mesh_plan.wave_mesh_sizes,
        mesh_plan.wave_decay_radii, mesh_plan.exterior_mesh_sizes)
        wave_size < remote_size || continue
        cable_curves = reduce(vcat, geometry.cable_curves; init=Int[])
        distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.set_numbers(distance, "CurvesList", unique(cable_curves))
        gmsh.model.mesh.field.set_number(distance, "Sampling", 200)
        sources = [distance]
        for (design,position) in zip(model.problem.system.designs,model.problem.system.positions)
            x,width = _interface_footprint(design,position)
            footprint = gmsh.model.mesh.field.add("MathEval")
            gmsh.model.mesh.field.set_string(footprint,"F",
                "Sqrt(y^2+Max(Abs(x-($x))-($(model.interface_refinement_factor*width)),0)^2)")
            push!(sources,footprint)
        end
        distance = gmsh.model.mesh.field.add("Min")
        gmsh.model.mesh.field.set_numbers(distance,"FieldsList",sources)
        threshold = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.set_number(threshold, "InField", distance)
        gmsh.model.mesh.field.set_number(threshold, "SizeMin", Float64(wave_size))
        gmsh.model.mesh.field.set_number(threshold, "SizeMax", Float64(remote_size))
        gmsh.model.mesh.field.set_number(threshold, "DistMin", Float64(decay_radius))
        gmsh.model.mesh.field.set_number(threshold, "DistMax", Float64(2decay_radius))
        restricted = gmsh.model.mesh.field.add("Restrict")
        gmsh.model.mesh.field.set_number(restricted, "InField", threshold)
        gmsh.model.mesh.field.set_numbers(restricted, "SurfacesList", surfaces)
        gmsh.model.mesh.field.set_number(restricted, "IncludeBoundary", 1)
        push!(transition_fields, restricted)
    end
    _configure_conductor_mesh!(model, geometry, mesh_plan, transition_fields)
    # Extend the actual boundary-edge sizes into cable insulation. Keep this
    # local: global boundary-size extension would also refine the remote domain.
    for cable_index in eachindex(model.cable_boundaries)
        surfaces = reduce(vcat, (geometry.region_surfaces[index]
            for (index,region) in enumerate(model.region_plans)
            if region.cable_index == cable_index &&
                model.material_plans[region.material_index].kind !== :conductor); init=Int[])
        isempty(surfaces) && continue
        extension = gmsh.model.mesh.field.add("Extend")
        gmsh.model.mesh.field.set_numbers(extension, "CurvesList", _entity_boundary(surfaces))
        size = model.cable_outer_mesh_sizes[cable_index]
        gmsh.model.mesh.field.set_number(extension, "SizeMax", size)
        gmsh.model.mesh.field.set_number(extension, "DistMax", size/(model.mesh_growth_factor-1))
        restricted = gmsh.model.mesh.field.add("Restrict")
        gmsh.model.mesh.field.set_number(restricted, "InField", extension)
        gmsh.model.mesh.field.set_numbers(restricted, "SurfacesList", surfaces)
        push!(transition_fields, restricted)
    end
    isempty(transition_fields) && return nothing
    background = if length(transition_fields) == 1
        only(transition_fields)
    else
        combined = gmsh.model.mesh.field.add("Min")
        gmsh.model.mesh.field.set_numbers(
            combined, "FieldsList", transition_fields
        )
        combined
    end
    gmsh.model.mesh.field.set_as_background_mesh(background)
    return nothing
end

function _configure_conductor_mesh!(model, geometry, plan, background_fields)
    empty!(geometry.conductor_fields)
    controls = model.conductor_mesh
    segments = _conductor_circle_segments(controls.geometry_tolerance)
    all_surfaces = last.(gmsh.model.get_entities(2))
    for (index, region) in enumerate(model.region_plans)
        sizes = _conductor_mesh_sizes(model, region, plan)
        sizes === nothing && continue
        surfaces = geometry.region_surfaces[index]
        curves = _entity_boundary(surfaces)
        partition = get(geometry.sector_partitions,index,nothing)
        layer = 0
        if partition !== nothing
            # Sectors use twice the round tangential resolution, retaining
            # the qualified straight-side cap as well as exact arc curvature.
            sector_segments = 2segments
            radius = region.shape.primitive.r_back
            for (outer,inner,length,turn) in partition.curve_pairs
                count = max(2,ceil(Int,turn*sector_segments/(2π)),
                    ceil(Int,length/(radius/5*96/sector_segments))) + 1
                for curve in (outer,inner)
                    gmsh.model.mesh.set_transfinite_curve(curve,count)
                    geometry.transfinite_curves[curve] = (count,1.0)
                end
            end
            intervals = controls.growth == 1 ? ceil(Int,partition.depth/sizes.first_size) :
                ceil(Int,log1p((controls.growth-1)*partition.depth/sizes.first_size)/log(controls.growth))
            for spoke in partition.spokes
                ratio = spoke > 0 ? controls.growth : inv(controls.growth)
                gmsh.model.mesh.set_transfinite_curve(abs(spoke),intervals+1,"Progression",ratio)
                geometry.transfinite_curves[abs(spoke)] = (intervals+1,ratio)
            end
        else
            for curve in curves
                arc = _conductor_curve_geometry(region.shape, curve)
                intervals = segments * arc.fraction
                # Angular fidelity is a lower bound on resolution, not a reason
                # to discard the user's existing local characteristic length.
                count = max(2, ceil(Int, intervals - 64eps(intervals)) + 1,
                    ceil(Int, arc.length/region.mesh_size) + 1)
                gmsh.model.mesh.set_transfinite_curve(curve, count)
                geometry.transfinite_curves[curve] = (count, 1.0)
            end
            layer = gmsh.model.mesh.field.add("BoundaryLayer")
            gmsh.model.mesh.field.set_numbers(layer, "CurvesList", curves)
            # Gmsh retains activated boundary-layer IDs after field removal. An
            # inactive replacement must exclude every surface, including when a
            # frequency scan reuses an ID that was previously active.
            gmsh.model.mesh.field.set_numbers(layer, "ExcludedSurfacesList",
                sizes.active ? setdiff(all_surfaces, surfaces) : all_surfaces)
            gmsh.model.mesh.field.set_number(layer, "Size", sizes.first_size)
            gmsh.model.mesh.field.set_number(layer, "Ratio", controls.growth)
            gmsh.model.mesh.field.set_number(layer, "Thickness", sizes.extent * (1+1e-8))
            gmsh.model.mesh.field.set_number(layer, "Quads", 0)
            sizes.active && gmsh.model.mesh.field.set_as_boundary_layer(layer)
        end
        bulk = gmsh.model.mesh.field.add("MathEval")
        gmsh.model.mesh.field.set_string(bulk, "F", string(sizes.bulk))
        restricted = gmsh.model.mesh.field.add("Restrict")
        gmsh.model.mesh.field.set_number(restricted, "InField", bulk)
        gmsh.model.mesh.field.set_numbers(restricted, "SurfacesList",
            partition === nothing ? surfaces : [partition.core])
        push!(background_fields, restricted)
        geometry.conductor_fields[index] = (layer, bulk)
    end
    # Local conductor fields own the small sizes. The former global floor
    # would erase them, including on the conductor side of shared interfaces.
    isempty(geometry.conductor_fields) || gmsh.option.set_number("Mesh.MeshSizeMin", 0.0)
    return nothing
end

function _mesh_metadata(
        model::FEMResolvedModel,
        fingerprint::String,
        gmsh_version::String,
        source::Symbol,
        mesh_plan::FEMMeshPlan = last(model.mesh_plans)
)
    return (
        schema = "LineCableModels.FEMMesh",
        version = 1,
        fingerprint,
        gmsh_version,
        source = String(source),
        mesh_dimension = 2,
        terminal_count = length(model.terminal_ids),
        terminal_ids = model.terminal_ids,
        physical_groups = [(dimension = dim, tag, name)
                           for (dim, tag, name) in _expected_physical_groups(model)],
        frequency_index = mesh_plan.frequency_index,
        frequency_hz = mesh_plan.frequency,
        physical_domain_halfwidth_m = mesh_plan.domain_halfwidth,
        pml_thickness_m = mesh_plan.pml_thickness,
        pml_element_family = mesh_plan.pml_element_family,
        pml_layers = mesh_plan.pml_layers,
        pml_grading = mesh_plan.pml_grading,
        pml_strips = mesh_plan.pml_strips,
        conductor_mesh = model.conductor_mesh,
        minimum_mesh_size_m = any(r -> _conductor_mesh_sizes(model,r,mesh_plan) !== nothing,
            model.region_plans) ? 0.0 : min(model.fine_mesh_size, minimum(mesh_plan.wave_mesh_sizes)),
        region_mesh_sizes_m = getproperty.(model.region_plans, :mesh_size),
        cable_outer_mesh_sizes_m = model.cable_outer_mesh_sizes,
        interface_mesh_size_m = mesh_plan.interface_mesh_size,
        cable_interface_mesh_sizes_m = mesh_plan.cable_interface_mesh_sizes,
        domain_mesh_size_m = mesh_plan.domain_mesh_size,
        exterior_mesh_sizes_m = mesh_plan.exterior_mesh_sizes,
        exterior_start_radius_m = mesh_plan.exterior_start_radius,
        wave_mesh_sizes_m = mesh_plan.wave_mesh_sizes,
        wave_decay_radii_m = mesh_plan.wave_decay_radii,
        adjacent_growth_factor = model.mesh_growth_factor
    )
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
    _copy_or_link_mesh(source, run_mesh)
    _write_json_atomic(metadata_path, metadata)
    return run_mesh
end

function _cache_mesh!(
        source::String,
        cache_mesh::String,
        cache_metadata::String,
        metadata
)
    mkpath(dirname(cache_mesh))
    temporary_mesh = tempname(dirname(cache_mesh))
    _copy_or_link_mesh(source, temporary_mesh)
    mv(temporary_mesh, cache_mesh; force = true)
    _write_json_atomic(cache_metadata, metadata)
    return cache_mesh
end

function _select_mesh!(
        run::FEMRun,
        model::FEMResolvedModel,
        geometry::FEMGeometry,
        execution::ComputationOptions,
        runtime_root::String,
        mesh_plan::FEMMeshPlan = last(model.mesh_plans)
)
    gmsh_version = _gmsh_version()
    fingerprint = _mesh_fingerprint(model, gmsh_version, mesh_plan)
    reference_mesh = mesh_plan.frequency_index == length(model.mesh_plans)
    reference_mesh && (run.mesh_fingerprint = fingerprint)
    cache_directory = joinpath(runtime_root, "meshes", fingerprint)
    cache_mesh = joinpath(cache_directory, "model.msh")
    cache_metadata = joinpath(cache_directory, "mesh.json")
    stem = reference_mesh ? "model" : @sprintf(
        "frequency_%04d", mesh_plan.frequency_index
    )
    run_mesh = joinpath(run.path, "mesh", "$stem.msh")
    run_metadata = joinpath(run.path, "mesh", "$stem.json")

    selected = nothing
    source = :generated
    if execution.data.mesh_policy === :reuse && isfile(run_mesh) && isfile(run_metadata)
        try
            metadata = JSON3.read(read(run_metadata, String))
            String(metadata.fingerprint) == fingerprint || error("fingerprint mismatch")
            _validate_mesh_file(model, run_mesh)
            selected = run_mesh
            source = :resume
        catch exception
            @warn "Ignoring an invalid retained FEM mesh" run_mesh exception
        end
    end
    if selected === nothing && reference_mesh && execution.data.mesh_policy === :reuse &&
       execution.data.mesh_path !== nothing
        explicit = abspath(execution.data.mesh_path)
        _validate_mesh_file(model, explicit)
        selected = explicit
        source = :explicit
    elseif selected === nothing && execution.data.mesh_policy === :reuse &&
           isfile(cache_mesh) && isfile(cache_metadata)
        try
            metadata = JSON3.read(read(cache_metadata, String))
            String(metadata.fingerprint) == fingerprint || error("fingerprint mismatch")
            _validate_mesh_file(model, cache_mesh)
            selected = cache_mesh
            source = :cache
        catch exception
            @warn "Ignoring an invalid cached FEM mesh" cache_mesh exception
        end
    end

    if selected === nothing
        gmsh.model.set_current(geometry.model_name)
        _configure_mesh!(model, geometry, mesh_plan)
        gmsh.model.mesh.generate(2)
        mkpath(dirname(run_mesh))
        gmsh.write(run_mesh)
        _inspect_loaded_mesh(model, run_mesh)
        metadata = _mesh_metadata(
            model, fingerprint, gmsh_version, :generated, mesh_plan
        )
        _write_json_atomic(run_metadata, metadata)
        _cache_mesh!(run_mesh, cache_mesh, cache_metadata, metadata)
        source = :generated
    elseif selected != run_mesh
        metadata = _mesh_metadata(
            model, fingerprint, gmsh_version, source, mesh_plan
        )
        _copy_mesh_snapshot!(selected, run_mesh, run_metadata, metadata)
    end
    reference_mesh && (run.mesh_source = source)
    return run_mesh
end

function _update_exterior_mesh!(model, geometry, plan)
    gmsh.model.set_current(geometry.model_name)
    gmsh.model.mesh.clear()
    for tag in gmsh.model.mesh.field.list()
        gmsh.model.mesh.field.remove(tag)
    end
    # Rebuild only the small exterior domain. GEO coordinate transforms can
    # invalidate cached arcs elsewhere in the model; cable entities stay fixed.
    gmsh.model.remove_physical_groups()
    gmsh.model.mesh.remove_embedded(
        [(2, s) for s in Iterators.flatten(geometry.material_surfaces)], 1)
    surfaces = unique([geometry.air_surfaces; geometry.earth_surfaces])
    curves = geometry.exterior_curves
    gmsh.model.geo.remove([(2, tag) for tag in surfaces], false)
    gmsh.model.geo.remove([(1, tag) for tag in curves], false)
    gmsh.model.geo.remove([(0, tag) for tag in geometry.exterior_points], false)
    gmsh.model.geo.synchronize()
    rebuilt = _build_geometry!(model, geometry.model_name, plan; reuse=geometry)
    for field in (:terminal_curves, :air_surfaces, :earth_surfaces, :pml_surfaces,
                  :outer_curves, :outer_air_curves, :outer_earth_curves,
                  :pml_inner_curves, :interface_curves, :exterior_curves, :exterior_points)
        target = getproperty(geometry, field)
        empty!(target)
        append!(target, getproperty(rebuilt, field))
    end
    for field in (:transfinite_curves,:transfinite_surfaces)
        target = getproperty(geometry,field)
        empty!(target)
        merge!(target,getproperty(rebuilt,field))
    end
    return nothing
end

function _select_meshes!(
        run::FEMRun,
        model::FEMResolvedModel,
        display_geometry::FEMGeometry,
        execution::ComputationOptions,
        runtime_root::String
)
    mesh_paths = Vector{String}(undef, length(model.mesh_plans))
    for mesh_plan in model.mesh_plans
        _update_exterior_mesh!(model, display_geometry, mesh_plan)
        mesh_paths[mesh_plan.frequency_index] = _select_mesh!(
            run, model, display_geometry, execution, runtime_root, mesh_plan
        )
    end
    return mesh_paths
end
