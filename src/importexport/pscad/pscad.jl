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
- `formulation`: Selected line-parameter or cable-constant formulation.
  The default selects lossless dielectric relations. Request `:Ametani2004`
  explicitly to include the supplied material losses.
- `temperature=nothing`: Optional operating temperature \\[°C\\]. Correction is
  applied here during export, not in geometric flattening; `nothing` retains
  the material reference temperatures.

!!! note
    PSCAD uses a reference-frequency equivalent loss tangent, bounded at ten.
    It does not reproduce an arbitrary broadband constitutive law. Exporting
    the selected relation matches its radial admittance at `base_freq` before
    that cap, not necessarily at every frequency in a later PSCAD scan.

# Returns

The written path. Filesystem errors are propagated to the caller.
"""
function export_data(
        ::Val{:pscad},
        system::LineCableSystem,
        earth::EarthModel;
        base_freq::Real = 50.0,
        formulation::Union{Engine.LineParametersFormulation,
            Engine.CableConstantsFormulation} = Engine.Formulation(),
        temperature::Union{Nothing, Real} = nothing,
        file_name::Union{AbstractString, Nothing} = nothing
)
    isfinite(base_freq) && base_freq > zero(base_freq) || throw(DomainError(
        base_freq, "PSCAD base frequency must be positive and finite"
    ))
    path = _pscad_output_path(system, file_name)
    #! explicit-imports: off
    # EzXML does not mark XMLError public, but this exporter preserves the
    # established exception contract for invalid XML output destinations.
    isdir(path) && throw(EzXML.XMLError(
        8,
        0,
        "PSCAD output path is a directory: $path",
        2,
        0
    ))
    #! explicit-imports: on
    document = _pscad_project(system, earth, base_freq; formulation, temperature)
    write(path, document)
    return path
end
