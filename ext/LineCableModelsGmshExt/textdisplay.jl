TextDisplay.@showfields LineCableModelsFEMError "LineCableModelsFEMError" error -> (
    category = error.category,
    object = error.object_id,
    field = error.field,
    message = error.message,
    run_directory = error.run_directory
)

TextDisplay.@showfields FEMElementBlock "FEMElementBlock" block -> (
    element_type = block.element_type,
    dimension = block.dimension,
    order = block.order,
    entity = block.entity,
    elements = length(block.element_tags),
    physical_tags = block.physical_tags
)

TextDisplay.@showfields FEMMesh "FEMMesh" mesh -> (
    source = mesh.source,
    nodes = length(mesh.node_tags),
    elements = sum(block -> length(block.element_tags), mesh.blocks; init = 0),
    blocks = length(mesh.blocks),
    physical_groups = length(mesh.physical_names)
)

TextDisplay.@showfields FEMFieldBlock "FEMFieldBlock" block -> (
    element_type = block.element_type,
    components = size(block.values, 1),
    nodes_per_element = size(block.values, 2),
    steps = size(block.values, 3),
    elements = size(block.values, 4)
)

TextDisplay.@showfields FEMFieldMap "FEMFieldMap" field -> (
    source = field.source,
    label = field.label,
    blocks = length(field.blocks),
    steps = length(field.times),
    representation = field.representation
)

TextDisplay.name(::Type{<:LineCableModelsFEM}) = "LineCableModelsFEM"
Base.summary(io::IO, ::LineCableModelsFEM) = print(io, "LineCableModels FEM backend")
function Base.show(io::IO, backend::LineCableModelsFEM)
    print(io, "LineCableModelsFEM(", length(backend.methods), " material laws)")
end
function Base.show(io::IO, ::MIME"text/plain", backend::LineCableModelsFEM)
    get(io, :compact, false) && return show(io, backend)
    selections = map(backend.methods) do selected
        selected === nothing && return nothing
        description(selected;compact=true)
    end
    return TextDisplay.fields(
        io,
        "LineCableModels FEM backend",
        (; selections..., options = backend.options);
        multiline = true
    )
end

