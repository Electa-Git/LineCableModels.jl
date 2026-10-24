include("schema.jl")
include("project.jl")
include("import.jl")

"""
$(TYPEDSIGNATURES)

Write a [`LineCableSystem`](@ref) as a minimal PSCAD project, the document of
[`pscad_project`](@ref LineCableModels.ImportExport.pscad_project).

# Arguments

- `system`: materialized line and cable geometry.
- `earth`: physical air and one infinite homogeneous soil half-space.
- `base_freq`: base frequency in hertz.
- `file_name`: destination `.pscx` file. The system identifier is prepended to
  an explicitly supplied basename.
- `formulation`: selected line-parameter or cable-constant formulation. Its insulation,
  semicon and temperature laws evaluate the materials. The default routes to the
  lossless dielectric relation. Request `:lossy` explicitly to include the supplied
  material losses.
- `temperature=nothing`: optional operating temperature \\[°C\\]. Correction is
  applied here during export, not in geometric flattening. `nothing` retains
  the material reference temperatures.

# Returns

The written path. Filesystem errors are propagated to the caller.
"""
function export_data(
        ::Val{:pscad},
        system::LineCableSystem,
        earth::EarthModel;
        formulation::Union{Engine.LineParametersFormulation,
            Engine.CableConstantsFormulation} = Engine.Formulation(),
        base_freq::Real = 50.0,
        temperature::Union{Nothing, Real} = nothing,
        file_name::Union{AbstractString, Nothing} = nothing
)
    path = _pscad_output_path(system, file_name)
    #! explicit-imports: off
    # EzXML does not mark XMLError public, but this exporter preserves the
    # established IOError behavior for invalid XML output destinations.
    isdir(path) && throw(EzXML.XMLError(
        8,
        0,
        "PSCAD output path is a directory: $path",
        2,
        0
    ))
    #! explicit-imports: on
    write(path, pscad_project(system, earth, formulation.formulas; base_freq, temperature))
    return path
end
