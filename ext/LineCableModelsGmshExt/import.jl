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
    blocks = FEMElementBlock[]
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
                FEMElementBlock(kind, dimension, order, primary,
                    entity, copy(physical), UInt64.(element_tags), connectivity))
        end
    end
    isempty(blocks) && throw(ArgumentError("mesh contains no elements: $path"))
    names = Dict((Int(dim), Int(tag)) => gmsh.model.get_physical_name(dim, tag)
    for (dim, tag) in gmsh.model.get_physical_groups())
    coordinates = reshape(xyz .* coordinate_scale, 3, :)
    all(isfinite, coordinates) ||
        throw(ArgumentError("scaled mesh coordinates must be finite"))
    return FEMMesh(path, UInt64.(tags), coordinates, blocks, names)
end

"""
$(TYPEDSIGNATURES)

Read a native Gmsh mesh into a detached `FEMMesh`. Retain all element blocks,
node tags, coordinates, and physical groups. Load Gmsh before calling this method.

# Arguments

- `path`: A mesh file or a saved FEM run, including its `mesh` directory.

# Keywords

- `frequency_index=nothing`: Select a one-based frequency index from a saved
  FEM run. A file path locates that file's mesh directory. Selection uses the
  runner's `model.json` metadata and requires a retained mesh for the index.
  With `nothing`, read the supplied file, or `model.msh` for a directory input.
  In a saved run, `model.msh` corresponds to the last (highest) frequency.
- `coordinate_scale=1`: Convert file lengths to meters. Use `1e-3` for
  millimeter coordinates; LineCableModels files already use meters.

# Returns

- A `FEMMesh` whose `source` is the absolute path of the selected mesh.

# Notes

Frequency selection reads saved files. It does not generate a mesh or run a
solver. A standalone mesh is independent of sidecar metadata when the index is omitted.

# Errors

- `ArgumentError`: Invalid frequency index, missing or invalid saved frequency
  metadata, missing mesh, or invalid coordinate scale.

# Examples

```julia
mesh = import_data(:msh, run_directory; frequency_index=3)
```
"""
function ImportExport.import_data(::Val{:msh}, path::AbstractString;
        frequency_index = nothing, coordinate_scale = 1)
    scale = _fem_coordinate_scale(coordinate_scale)
    frequency_index === nothing ||
        (frequency_index isa Integer && !(frequency_index isa Bool) && frequency_index > 0) ||
        throw(ArgumentError("frequency_index must be a positive one-based integer or nothing"))
    filename = abspath(path)
    if isdir(filename)
        directory = isdir(joinpath(filename, "mesh")) ? joinpath(filename, "mesh") : filename
        filename = joinpath(directory, "model.msh")
    end
    if frequency_index !== nothing
        directory = dirname(filename)
        metadata_path = joinpath(directory, "model.json")
        isfile(metadata_path) || throw(ArgumentError(
            "frequency_index requires saved FEM mesh metadata: $metadata_path"))
        metadata = JSON3.read(read(metadata_path, String))
        reference_index = metadata isa AbstractDict &&
            get(metadata, :schema, nothing) == "LineCableModels.FEMMesh" ?
            get(metadata, :frequency_index, nothing) : nothing
        reference_index isa Integer && !(reference_index isa Bool) && reference_index > 0 ||
            throw(ArgumentError("invalid saved FEM mesh frequency index in $metadata_path"))
        frequency_index <= reference_index || throw(ArgumentError(
            "frequency_index must be in 1:$reference_index for this saved FEM run"))
        stem = frequency_index == reference_index ? "model" :
            "frequency_$(lpad(string(frequency_index), 4, '0'))"
        filename = joinpath(directory, "$stem.msh")
    end
    return _read_fem_file(filename) do selected, _
        _read_fem_mesh(selected, scale)
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
    blocks = FEMFieldBlock[]
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
        push!(blocks, FEMFieldBlock(element, coordinates, values))
    end
    isempty(blocks) &&
        throw(ArgumentError("field view contains no sampled elements: $path"))
    return FEMFieldMap(path, label, blocks, times, encoding)
end

"""
    import_data(:pos, path; view=nothing, representation=:auto, coordinate_scale=1)

Read native list-based Gmsh field views into detached `FEMFieldMap` values.
Return one map for a single view, or a vector for multiple views. `view` selects
a one-based view index. `coordinate_scale` converts file lengths to meters.
Preserve element-local values, native output-step times, and the original label,
including any physical units. No GetDP process or graphical interface is used.

`representation=:auto` recognizes the explicit `phasor=real,imag` label emitted
by this FEM backend. other views retain independent real steps. Pass
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
