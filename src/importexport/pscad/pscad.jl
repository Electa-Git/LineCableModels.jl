include("schema.jl")
include("project.jl")
include("import.jl")

"""
$(TYPEDSIGNATURES)

Export a [`LineCableSystem`](@ref) as a minimal PSCAD project.

The generated project preserves PSCAD's `master:Line_FrePhase_Options`,
`master:Cable_Coax`, and `master:Line_Ground` component bindings. Cable
geometry, material properties, dielectric losses, phase eliminations, line
length, base frequency, and static earth properties are emitted as component
parameters. PSCAD may normalise the deterministic placeholder identifiers when
it opens the project.

# Arguments

- `system`: Materialised line and cable geometry.
- `earth`: Static earth properties. The final layer supplies ground values.
- `base_freq`: Base frequency in hertz.
- `file_name`: Destination `.pscx` file. The system identifier is prepended to
  an explicitly supplied basename.

# Returns

The written path. Filesystem errors are propagated to the caller.
"""
function export_data(
        ::Val{:pscad},
        system::LineCableSystem,
        earth::EarthModel;
        base_freq::Real = 50.0,
        file_name::Union{AbstractString, Nothing} = nothing
)
    base_freq > zero(base_freq) || throw(DomainError(
        base_freq, "PSCAD base frequency must be positive"
    ))
    path = _pscad_output_path(system, file_name)
    isdir(path) && throw(EzXML.XMLError(
        8,
        0,
        "PSCAD output path is a directory: $path",
        2,
        0
    ))
    document = _pscad_project(system, earth, base_freq)
    write(path, document)
    return path
end
