include("schema.jl")
include("project.jl")
include("import.jl")

"""
$(TYPEDSIGNATURES)

Export a cable transient problem as a PSCAD project with sources, terminations
and voltage and current channels. Unlisted terminals remain open. A `Matched()`
end connects to an open extension of length `1.1 * c * tmax / 2` \\[m\\].
Here `c` is the speed of light \\[m/s\\], and `tmax` is the duration \\[s\\].

# Arguments

- `problem`: Transient analysis of an uncomputed line-parameter problem.
- `formulation`: Settings for PSCAD transient export.

# Keywords

- `file_name`: Destination `.pscx` file. Its basename receives the system identifier.

# Returns

- Written project path.

# Notes

The window sets the duration and integration step \\[s\\]. Voltage and current
channels use terminal labels with units \\[V\\] and \\[A\\], respectively.

# Errors

Throw `ArgumentError` for unsupported geometry, earth, abscissa or settings.
"""
function export_data(::Val{:pscad},
        problem::TransientAnalysisProblem{<:LineParametersProblem},
        formulation::PSCADTransientFormulation; file_name)
    validate(formulation)
    validate(problem, formulation)
    line = problem.line
    system = line.system
    blueprints = _pscad_blueprints(system)
    settings = (duration = problem.window.tmax, step = problem.window.step,
        sample = something(formulation.sample, problem.window.step),
        output_stem = formulation.output_stem)
    document = _pscad_project(system, line.earth_props, formulation.frequency, blueprints,
        formulation.line, line.temperature; settings,
        model = formulation.model, fit = formulation.fit)
    for parameter in findall("//User[@defn='master:Cable_Coax']/paramlist/param[starts-with(@name,'elim')]", document)
        parameter["value"] = "0"
    end
    schematic = only(findall("/project/definitions/Definition[@name='Main']/schematic", document))
    canvas = PSCADCanvas(schematic)
    canvas.next_cell = 6 # Reserve the first row for the cable wrapper.
    _pscad_transient_end!(canvas, problem.sending, problem, "sending", formulation)
    _pscad_transient_end!(canvas, problem.receiving, problem, "receiving", formulation)
    path = _pscad_output_path(system, file_name)
    write(path, document)
    return path
end

"""
$(TYPEDSIGNATURES)

Export a [`LineCableSystem`](@ref) as a minimal PSCAD project.

The generated project preserves PSCAD's `master:Line_FrePhase_Options`,
`master:Cable_Coax`, and `master:Line_Ground` component bindings. Cable
geometry, material properties, dielectric losses, phase eliminations, line
length, base frequency, and static earth properties are emitted as component
parameters. PSCAD may normalize the deterministic placeholder identifiers when
it opens the project.

# Arguments

- `system`: materialized line and cable geometry.
- `earth`: physical air and one infinite homogeneous soil half-space.
- `base_freq`: base frequency in hertz.
- `file_name`: destination `.pscx` file. The system identifier is prepended to
  an explicitly supplied basename.
- `formulation`: selected line-parameter or cable-constant formulation.
  The default routes to the lossless dielectric relation. Request `:lossy`
  explicitly to include the supplied material losses.
- `native_settings=(;)`: optional validated native `ground` and `frequency`
  field records supplied by the PSCAD formula adapter. Each field includes its
  `value` and expected `readback`. The same record accompanies native execution.
- `temperature=nothing`: optional operating temperature \\[°C\\]. Correction is
  applied here during export, not in geometric flattening. `nothing` retains
  the material reference temperatures.

Each coaxial group of each design becomes one `Cable_Coax` at the group's position.
Each conductor writes its resistivity at `temperature`. The exporter evaluates its
dielectric group at `base_freq` and `temperature`. `EPS` is the group's equivalent
relative permittivity, and `LT` its ratio G/(ωC) of conductance to susceptance.

!!! note
    PSCAD reads the conductance of a group as ωC·LT, and the exporter caps `LT` at
    ten. Where G/(ωC) exceeds ten, the exporter writes ten, and PSCAD receives
    10ωC: a fraction 10ωC/G of the group's conductance. This happens only where
    conduction dominates the group, as in an interval of semicon alone or at very low
    frequency. Below the bound, the exported group matches the selected laws at
    `base_freq`, not necessarily at every frequency of a later PSCAD scan.

# Returns

The written path. Filesystem errors are propagated to the caller.
"""
function export_data(
        ::Val{:pscad},
        system::LineCableSystem,
        earth::EarthModel;
        formulation::Union{Engine.LineParametersFormulation,
            Engine.CableConstantsFormulation, PSCADFormulation} = Engine.Formulation(),
        base_freq::Real = formulation isa PSCADFormulation ?
                          formulation.options.data.base_frequency : 50.0,
        temperature::Union{Nothing, Real} = nothing,
        file_name::Union{AbstractString, Nothing} = nothing,
        native_settings::NamedTuple = (;)
)
    _pscad_deterministic(eltype(system), eltype(earth), typeof(base_freq), typeof(temperature))
    isfinite(base_freq) && base_freq > zero(base_freq) || throw(DomainError(
        base_freq, "PSCAD base frequency must be positive and finite"
    ))
    validate(earth, PSCADFormulation)
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
    blueprints = _pscad_blueprints(system)
    document = _pscad_project(system, earth, base_freq, blueprints, formulation,
        temperature; native_settings)
    write(path, document)
    return path
end
