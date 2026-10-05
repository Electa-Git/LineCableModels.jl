"""
$(TYPEDEF)

Retain a homogeneous block of native finite elements and its physical groups.
Connectivity indexes columns of the enclosing [`FEMMesh`](@ref) coordinates.
Native node and element tags remain separate from these dense indices.

$(TYPEDFIELDS)
"""
struct FEMElementBlock
    "Native Gmsh element type."
    element_type::Int
    "Topological dimension."
    dimension::Int
    "Polynomial geometry order."
    order::Int
    "Number of corner nodes per element."
    primary_nodes::Int
    "Native entity tag."
    entity::Int
    "Every physical group containing this entity."
    physical_tags::Vector{Int}
    "Native element tags."
    element_tags::Vector{UInt64}
    "Node indices, with one element per column."
    connectivity::Matrix{Int}
end

"""
$(TYPEDEF)

Retain a mesh independently of a Gmsh session or a solver run. Import through
`import_data(:msh, path)` with Gmsh loaded. Coordinates preserve the file's
Cartesian frame. `coordinate_scale` at import converts file lengths to meters.

$(TYPEDFIELDS)
"""
struct FEMMesh
    "Absolute source filename."
    source::String
    "Native node tags corresponding to coordinate columns."
    node_tags::Vector{UInt64}
    "Cartesian coordinates \\[m\\], three rows and one column per node."
    coordinates::Matrix{Float64}
    "Element blocks, retaining all dimensions and orders."
    blocks::Vector{FEMElementBlock}
    "Physical names keyed by (dimension, physical tag)."
    physical_names::Dict{Tuple{Int, Int}, String}
    "Passive saved-run case metadata, or `nothing` for a standalone mesh."
    provenance::Union{Nothing,NamedTuple}
end

"""
$(TYPEDEF)

Retain element-local field samples without averaging coincident vertices.
Coordinates are stored in meters. sample units are those recorded by the
source field map.

$(TYPEDFIELDS)
"""
struct FEMFieldBlock
    "Native Gmsh element type for the sampled geometry."
    element_type::Int
    "Cartesian coordinates \\[m\\] indexed by (coordinate, node, element)."
    coordinates::Array{Float64, 3}
    "Samples indexed by (component, node, output step, element)."
    values::Array{Float64, 4}
end

"""
$(TYPEDEF)

Retain one native field view independently of its reader session. The label
preserves recorded units, excitation, frequency, and any PML qualification.
`representation=:complex` identifies exactly two output steps as the real and
imaginary parts of a phasor. `:real` keeps steps independent.

$(TYPEDFIELDS)
"""
struct FEMFieldMap
    "Absolute source filename."
    source::String
    "Unmodified native view label, including recorded physical units."
    label::String
    "Element-local geometry and samples."
    blocks::Vector{FEMFieldBlock}
    "Native output-step times. these need not denote physical time or frequency."
    times::Vector{Float64}
    "Either independent real steps or a real and imaginary phasor pair."
    representation::Symbol
end
