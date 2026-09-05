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
        system = ImportExport.serialize_value(model.problem.system),
        temperature = ImportExport.serialize_value(model.problem.temperature),
        earth_props = ImportExport.serialize_value(model.problem.earth_props),
        terminal_ids = model.terminal_ids,
        material_tags = getproperty.(model.material_plans, :physical_tag),
        material_names = getproperty.(model.material_plans, :physical_name),
        tags = model.tags,
        region_mesh_sizes = getproperty.(model.region_plans, :mesh_size),
        cable_outer_mesh_sizes = model.cable_outer_mesh_sizes,
        mesh_growth_factor = model.mesh_growth_factor,
        mesh_frequency = mesh_plan.frequency,
        domain_radius = mesh_plan.domain_radius,
        shell_outer_radius = mesh_plan.shell_outer_radius,
        domain_mesh_size = mesh_plan.domain_mesh_size,
        infinite_mesh_size = mesh_plan.infinite_mesh_size,
        interface_mesh_size = mesh_plan.interface_mesh_size,
        cable_interface_mesh_sizes = mesh_plan.cable_interface_mesh_sizes,
        gmsh_version
    )
    return bytes2hex(sha256(JSON3.write(evidence)))
end

function _expected_physical_groups(model::FEMResolvedModel)
    groups = Tuple{Int, Int, String}[
        (
            2, model.tags.air, "LCM/domain/air"),
        (
            2, model.tags.earth, "LCM/domain/earth"),
        (
            2, model.tags.air_infinite, "LCM/domain/air_infinite"),
        (
            2, model.tags.earth_infinite, "LCM/domain/earth_infinite"),
        (
            2, model.tags.infinite_domain, "LCM/domain/infinite_shell"),
        (
            1,
            model.tags.outer_boundary,
            "LCM/boundary/magnetic_dirichlet"),
        (
            1,
            model.tags.outer_air_boundary,
            "LCM/boundary/electric_insulation_air"),
        (
            1,
            model.tags.outer_earth_boundary,
            "LCM/boundary/electric_reference_earth"),
        (
            1,
            model.tags.inner_shell_boundary,
            "LCM/boundary/inner_infinite_shell"
        ),
        (
            1, model.tags.interface, "LCM/interface/air_earth"),
        (
            2, 6_001, "LCM/domain/conductors"),
        (
            2, 6_003, "LCM/domain/field_maps")
    ]
    for (index, name) in enumerate(model.terminal_names)
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

function _configure_mesh!(
        model::FEMResolvedModel,
        geometry::FEMGeometry,
        mesh_plan::FEMMeshPlan = last(model.mesh_plans)
)
    gmsh.option.set_number("Mesh.MshFileVersion", 4.1)
    gmsh.option.set_number("Mesh.MeshSizeMin", Float64(model.fine_mesh_size))
    gmsh.option.set_number(
        "Mesh.MeshSizeMax", Float64(mesh_plan.infinite_mesh_size)
    )
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
            threshold, "SizeMax", Float64(mesh_plan.domain_mesh_size)
        )
        gmsh.model.mesh.field.set_number(threshold, "DistMin", 0.0)
        transition_distance = max(
            cable_size,
            (mesh_plan.domain_mesh_size - cable_size) /
            max(model.mesh_growth_factor - one(model.mesh_growth_factor), 1e-12)
        )
        gmsh.model.mesh.field.set_number(
            threshold, "DistMax", Float64(transition_distance)
        )
        push!(transition_fields, threshold)
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
        inner_shell_radius_m = mesh_plan.domain_radius,
        outer_shell_radius_m = mesh_plan.shell_outer_radius,
        minimum_mesh_size_m = model.fine_mesh_size,
        region_mesh_sizes_m = getproperty.(model.region_plans, :mesh_size),
        cable_outer_mesh_sizes_m = model.cable_outer_mesh_sizes,
        interface_mesh_size_m = mesh_plan.interface_mesh_size,
        cable_interface_mesh_sizes_m = mesh_plan.cable_interface_mesh_sizes,
        domain_mesh_size_m = mesh_plan.domain_mesh_size,
        infinite_mesh_size_m = mesh_plan.infinite_mesh_size,
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

function _load_mesh_for_display!(mesh_path::String)
    name = "LineCableModelsFEM-mesh-$(time_ns())"
    gmsh.model.add(name)
    gmsh.merge(mesh_path)
    return name
end

function _select_mesh!(
        run::FEMRun,
        model::FEMResolvedModel,
        geometry::FEMGeometry,
        formulation::LineCableModelsFEM,
        runtime_root::String,
        mesh_plan::FEMMeshPlan = last(model.mesh_plans)
)
    gmsh_version = _gmsh_version()
    fingerprint = _mesh_fingerprint(model, gmsh_version, mesh_plan)
    displayed = mesh_plan.frequency_index == length(model.mesh_plans)
    displayed && (run.mesh_fingerprint = fingerprint)
    cache_directory = joinpath(runtime_root, "meshes", fingerprint)
    cache_mesh = joinpath(cache_directory, "model.msh")
    cache_metadata = joinpath(cache_directory, "mesh.json")
    stem = displayed ? "model" : @sprintf(
        "frequency_%04d", mesh_plan.frequency_index
    )
    run_mesh = joinpath(run.path, "mesh", "$stem.msh")
    run_metadata = joinpath(run.path, "mesh", "$stem.json")
    execution = formulation.execution

    selected = nothing
    source = :generated
    if execution.mesh_policy === :reuse && isfile(run_mesh) && isfile(run_metadata)
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
    if selected === nothing && displayed && execution.mesh_policy === :reuse &&
       execution.mesh_path !== nothing
        explicit = abspath(execution.mesh_path)
        _validate_mesh_file(model, explicit)
        selected = explicit
        source = :explicit
    elseif selected === nothing && execution.mesh_policy === :reuse &&
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
        formulation.execution.ui && displayed && _load_mesh_for_display!(run_mesh)
    elseif formulation.execution.ui && displayed
        _load_mesh_for_display!(run_mesh)
    end
    displayed && (run.mesh_source = source)
    return run_mesh
end

function _select_meshes!(
        run::FEMRun,
        model::FEMResolvedModel,
        display_geometry::FEMGeometry,
        formulation::LineCableModelsFEM,
        runtime_root::String
)
    mesh_paths = Vector{String}(undef, length(model.mesh_plans))
    for mesh_plan in model.mesh_plans
        geometry = if mesh_plan.frequency_index == length(model.mesh_plans)
            display_geometry
        else
            _build_geometry!(
                model,
                @sprintf(
                    "LineCableModelsFEM-%s-f%04d",
                    basename(run.path),
                    mesh_plan.frequency_index
                ),
                mesh_plan
            )
        end
        mesh_paths[mesh_plan.frequency_index] = _select_mesh!(
            run,
            model,
            geometry,
            formulation,
            runtime_root,
            mesh_plan
        )
    end
    return mesh_paths
end
