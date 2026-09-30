# Native file I/O shares the Gmsh session owner, but never prepares a solve.
function _read_fem_file(read_data, filename::AbstractString)
    path = abspath(filename)
    isfile(path) || throw(ArgumentError("FEM file does not exist: $path"))
    return lock(FEM_SESSION_LOCK) do
        session = _start_gmsh(0)
        try
            gmsh.model.add("LineCableModels-import-$(time_ns())")
            gmsh.merge(path)
            read_data(path, session)
        finally
            _finish_gmsh(session)
        end
    end
end

function _fem_coordinate_scale(value)
    value isa Real && !(value isa Bool) && isfinite(value) && value > 0 &&
    0 < Float64(value) < Inf ||
        throw(ArgumentError("coordinate_scale must be finite and positive"))
    return Float64(value)
end

function _read_fem_mesh(path, coordinate_scale)
    tags, xyz, _ = gmsh.model.mesh.get_nodes()
    isempty(tags) && throw(ArgumentError("mesh contains no nodes: $path"))
    allunique(tags) || throw(ArgumentError("mesh contains duplicate node tags: $path"))
    all(isfinite, xyz) || throw(ArgumentError("mesh coordinates must be finite: $path"))
    indices = Dict(tag => i for (i, tag) in enumerate(tags))
    blocks = Engine.FEMElementBlock[]
    for (dim, entity) in gmsh.model.get_entities()
        physical = Int.(gmsh.model.get_physical_groups_for_entity(dim, entity))
        kinds, elements, nodes = gmsh.model.mesh.get_elements(dim, entity)
        for (kind, element_tags, node_tags) in zip(kinds, elements, nodes)
            _, dimension, order, count,
            _, primary = gmsh.model.mesh.get_element_properties(kind)
            length(node_tags) == count * length(element_tags) ||
                throw(ArgumentError("invalid element connectivity in $path"))
            connectivity = reshape([indices[tag] for tag in node_tags], count, :)
            push!(blocks,
                Engine.FEMElementBlock(kind, dimension, order, primary,
                    entity, copy(physical), UInt64.(element_tags), connectivity))
        end
    end
    isempty(blocks) && throw(ArgumentError("mesh contains no elements: $path"))
    names = Dict((Int(dim), Int(tag)) => gmsh.model.get_physical_name(dim, tag)
    for (dim, tag) in gmsh.model.get_physical_groups())
    coordinates = reshape(xyz .* coordinate_scale, 3, :)
    all(isfinite, coordinates) ||
        throw(ArgumentError("scaled mesh coordinates must be finite"))
    return Engine.FEMMesh(path, UInt64.(tags), coordinates, blocks, names)
end

"""
    import_data(:msh, path; coordinate_scale=1)

Read a native Gmsh mesh into a detached `FEMMesh`. Retain all element blocks,
node tags, coordinates, and physical groups. `coordinate_scale` converts file
lengths to meters (use `1e-3` for millimeter coordinates). LineCableModels files
already use meters. No mesh generation or solver is invoked.
Load Gmsh before calling this method.
"""
function ImportExport.import_data(::Val{:msh}, path::AbstractString; coordinate_scale = 1)
    scale = _fem_coordinate_scale(coordinate_scale)
    return _read_fem_file(path) do filename, _
        _read_fem_mesh(filename, scale)
    end
end

# POS element codes are an external file-format contract, not a solver registry.
const _POS_ELEMENT_TYPES = Dict('P' => 15, 'L' => 1, 'T' => 2, 'Q' => 3,
    'S' => 4, 'H' => 5, 'I' => 6, 'Y' => 7)

function _read_fem_view(path, tag, representation, coordinate_scale)
    label = gmsh.view.option.get_string(tag, "Name")
    kinds, counts, data = gmsh.view.get_list_data(tag)
    isempty(kinds) && throw(ArgumentError(
        "view $(repr(label)) has no list-based field data; export a native POS view"))
    steps = Int(gmsh.view.option.get_number(tag, "NbTimeStep"))
    steps > 0 || throw(ArgumentError("field view has no output steps: $path"))
    times = map(0:(steps - 1)) do step
        gmsh.view.option.set_number(tag, "TimeStep", step)
        gmsh.view.option.get_number(tag, "Time")
    end
    encoding = representation === :auto ?
               (occursin("; phasor=real,imag", label) ? :complex : :real) : representation
    encoding === :complex && steps != 2 &&
        throw(ArgumentError(
            "complex field representation requires exactly two real/imaginary steps"))
    blocks = Engine.FEMFieldBlock[]
    for (kind, count, raw) in zip(kinds, counts, data)
        count == 0 && continue
        length(kind) == 2 && haskey(_POS_ELEMENT_TYPES, kind[2]) ||
            throw(ArgumentError("unsupported POS element encoding $kind"))
        components = kind[1] == 'S' ? 1 :
                     kind[1] == 'V' ? 3 :
                     kind[1] == 'T' ? 9 :
                     throw(ArgumentError("unsupported POS field encoding $kind"))
        element = _POS_ELEMENT_TYPES[kind[2]]
        _, _, _, nodes, _, _ = gmsh.model.mesh.get_element_properties(element)
        width = nodes * (3 + components * steps)
        length(raw) == count * width || throw(ArgumentError(
            "view $(repr(label)) uses unsupported higher-order interpolation"))
        packed = reshape(raw, width, count)
        coordinates = Array{Float64}(undef, 3, nodes, count)
        values = Array{Float64}(undef, components, nodes, steps, count)
        for i in 1:count
            coordinates[:, :, i] = transpose(reshape(packed[1:3nodes, i], nodes, 3))
            values[:, :, :, i] = reshape(packed[(3nodes + 1):end, i], components, nodes, steps)
        end
        coordinates .*= coordinate_scale
        all(isfinite, coordinates) ||
            throw(ArgumentError("field coordinates must be finite"))
        push!(blocks, Engine.FEMFieldBlock(element, coordinates, values))
    end
    isempty(blocks) &&
        throw(ArgumentError("field view contains no sampled elements: $path"))
    return Engine.FEMFieldMap(path, label, blocks, times, encoding)
end

"""
    import_data(:pos, path; view=nothing, representation=:auto, coordinate_scale=1)

Read native list-based Gmsh field views into detached `FEMFieldMap` values.
Return one map for a single view, or a vector for multiple views. `view` selects
a one-based view index. `coordinate_scale` converts file lengths to meters.
Preserve element-local values, native output-step times, and the original label,
including any physical units. No GetDP process or graphical interface is used.

`representation=:auto` recognizes the explicit `phasor=real,imag` label emitted
by this FEM backend; other views retain independent real steps. Pass
`representation=:complex` for older harmonic maps with two known real/imaginary
steps, or `:real` for independent steps. Merely having two steps never identifies
a complex field. Load Gmsh before importing.
"""
function ImportExport.import_data(::Val{:pos}, path::AbstractString;
        view = nothing, representation::Symbol = :auto, coordinate_scale = 1)
    scale = _fem_coordinate_scale(coordinate_scale)
    representation in (:auto, :real, :complex) || throw(ArgumentError(
        "field representation must be :auto, :real, or :complex"))
    view === nothing || (view isa Integer && !(view isa Bool) && view > 0) ||
        throw(ArgumentError("view must be a positive one-based index"))
    return _read_fem_file(path) do filename, session
        tags = sort!(collect(setdiff(Set(Int.(gmsh.view.get_tags())), session.initial_views)))
        isempty(tags) && throw(ArgumentError("file contains no field views: $filename"))
        view === nothing || view <= length(tags) ||
            throw(ArgumentError("field view index is out of range"))
        selected = view === nothing ? tags : [tags[view]]
        maps = [_read_fem_view(filename, tag, representation, scale) for tag in selected]
        length(maps) == 1 ? only(maps) : maps
    end
end
