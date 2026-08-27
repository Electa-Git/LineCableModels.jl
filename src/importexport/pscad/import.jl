const _PSCAD_EPSILON_0 = 8.8541878128e-12

const _PSCAD_PART_FIELDS = (
    (
        inner = "R1", outer = "R2", insulation = "R3",
        rho = "RHOC", conductor_mu = "PERMC",
        eps = "EPS1", insulation_mu = "PERM1", loss = "LT1"
    ),
    (
        inner = nothing, outer = "R4", insulation = "R5",
        rho = "RHOS", conductor_mu = "PERMS",
        eps = "EPS2", insulation_mu = "PERM2", loss = "LT2"
    ),
    (
        inner = nothing, outer = "R6", insulation = "R7",
        rho = "RHOA", conductor_mu = "PERMA",
        eps = "EPS3", insulation_mu = "PERM3", loss = "LT3"
    ),
    (
        inner = nothing, outer = "R8", insulation = "R9",
        rho = "RHOO", conductor_mu = "PERMO",
        eps = "EPS4", insulation_mu = "PERM4", loss = "LT4"
    )
)

function _pscad_parameters(node)
    values = Dict{String, String}()
    for parameter in findall("./paramlist/param", node)
        name = parameter["name"]
        haskey(values, name) && throw(ArgumentError(
            "duplicate PSCAD parameter '$name'",
        ))
        values[name] = parameter["value"]
    end
    return values
end

function _pscad_parameter(values, name::AbstractString)
    haskey(values, name) || throw(KeyError(name))
    return values[name]
end

function _pscad_number(values, name::AbstractString)
    raw = _pscad_parameter(values, name)
    scalar = strip(first(split(raw, '['; limit = 2)))
    scalar = replace(scalar, 'D' => 'E', 'd' => 'e')
    value = tryparse(Float64, scalar)
    value === nothing && throw(ArgumentError(
        "PSCAD parameter '$name' is not numeric: '$raw'",
    ))
    return value
end

function _pscad_integer(values, name::AbstractString)
    value = _pscad_number(values, name)
    isinteger(value) || throw(ArgumentError(
        "PSCAD parameter '$name' must be an integer; got $value",
    ))
    return round(Int, value)
end

function _pscad_binding(node)
    return haskey(node, "defn") ? node["defn"] : ""
end

function _pscad_output_enabled(node)
    values = _pscad_parameters(node)
    haskey(values, "Output") || return nothing
    return uppercase(strip(values["Output"])) in ("1", "YES", "ENABLED", "TRUE")
end

function _pscad_row_definition(
        project,
        requested::Union{Nothing, AbstractString} = nothing
)
    candidates = NamedTuple[]
    for definition in findall("./definitions/Definition", project)
        users = findall("./schematic/User", definition)
        frequency = filter(
            node -> _pscad_binding(node) == _PSCAD_FREQUENCY_BINDING,
            users
        )
        ground = filter(
            node -> _pscad_binding(node) == _PSCAD_GROUND_BINDING,
            users
        )
        isempty(frequency) && continue
        isempty(ground) && continue
        length(frequency) == 1 || throw(ArgumentError(
            "PSCAD definition '$(definition["name"])' has multiple frequency options",
        ))
        length(ground) == 1 || throw(ArgumentError(
            "PSCAD definition '$(definition["name"])' has multiple ground definitions",
        ))
        push!(candidates, (;
            definition,
            users,
            frequency = only(frequency),
            ground = only(ground)
        ))
    end
    isempty(candidates) && throw(ArgumentError(
        "PSCAD project contains no supported frequency-dependent row definition",
    ))
    if requested !== nothing
        qualified = String(requested)
        name = last(split(qualified, ':'; limit = 2))
        selected = filter(candidate -> candidate.definition["name"] == name, candidates)
        length(selected) == 1 || throw(ArgumentError(
            isempty(selected) ?
            "PSCAD project contains no supported row definition '$qualified'" :
            "PSCAD project contains multiple supported row definitions named '$name'",
        ))
        return only(selected)
    end
    enabled = filter(candidate -> _pscad_output_enabled(candidate.frequency) === true,
        candidates)
    if length(enabled) == 1
        return only(enabled)
    end
    isempty(enabled) && length(candidates) == 1 && return only(candidates)
    throw(ArgumentError(
        isempty(enabled) ?
        "PSCAD project has multiple row definitions and none is selected for output" :
        "PSCAD project has multiple row definitions selected for output",
    ))
end

function _pscad_active_cables(row)
    unsupported = String[]
    cables = typeof(row.definition)[]
    for node in row.users
        binding = _pscad_binding(node)
        if binding in (_PSCAD_CABLE_BINDING, "master:Cable_CoaxSimpl")
            values = _pscad_parameters(node)
            _pscad_integer(values, "LL") >= 0 && push!(cables, node)
        elseif startswith(binding, "master:Line_Tower_") ||
               binding == "master:Cable_PipeType"
            push!(unsupported, binding)
        end
    end
    isempty(unsupported) || throw(ArgumentError(
        "PSCAD row definition uses unsupported physical components: " *
        join(unique(unsupported), ", "),
    ))
    isempty(cables) && throw(ArgumentError(
        "PSCAD row definition contains no active master:Cable_Coax component",
    ))
    sort!(cables; by = node -> _pscad_integer(_pscad_parameters(node), "CABNUM"))
    return cables
end

function _pscad_instance(document, project, row)
    namespace = project["name"]
    target = "$namespace:$(row.definition["name"])"
    candidates = filter(findall("//User", document)) do node
        _pscad_binding(node) == target || return false
        values = _pscad_parameters(node)
        return haskey(values, "Length") && haskey(values, "Name")
    end
    length(candidates) == 1 || throw(ArgumentError(
        "PSCAD project must contain one instance of '$target' with line parameters",
    ))
    return only(candidates)
end

function _pscad_material(kind::Symbol, rho, eps_r, mu_r)
    return Material(kind, rho, eps_r, mu_r, 20.0, 0.0)
end

function _pscad_dielectric(values, fields, frequency)
    eps_r = _pscad_number(values, fields.eps)
    mu_r = _pscad_number(values, fields.insulation_mu)
    loss = _pscad_number(values, fields.loss)
    0 <= loss <= 10 || throw(DomainError(
        loss, "PSCAD loss tangent must be between zero and ten"
    ))
    conductivity = 2π * frequency * _PSCAD_EPSILON_0 * eps_r * loss
    rho = iszero(conductivity) ? Inf : inv(conductivity)
    return _pscad_material(:insulator, rho, eps_r, mu_r)
end

function _pscad_design(values, cable_number::Int, frequency)
    line_layers = _pscad_integer(values, "LL")
    if iszero(line_layers)
        conductor_inner = _pscad_number(values, "R1")
        conductor_outer = _pscad_number(values, "R2")
        conductor = ConductorGroup(Tubular(
            conductor_inner,
            conductor_outer,
            _pscad_material(
                :conductor,
                _pscad_number(values, "RHOC"),
                0.0,
                _pscad_number(values, "PERMC")
            )
        ))
        air_shell = max(abs(conductor_outer) * 1e-6, 1e-9)
        insulation = InsulatorGroup(
            Insulator(
                conductor_outer,
                conductor_outer + air_shell,
                _pscad_material(:insulator, Inf, 1.0, 1.0)
            );
            reference_frequency = frequency
        )
        cable_name = strip(get(values, "Name", ""))
        isempty(cable_name) && (cable_name = "cable$cable_number")
        return CableDesign(
            cable_name,
            CableComponent("conductor", conductor, insulation)
        )
    end
    isodd(line_layers) || throw(ArgumentError(
        "PSCAD LL must describe alternating conductor and insulation layers",
    ))
    component_count = (line_layers + 1) ÷ 2
    component_count in eachindex(_PSCAD_PART_FIELDS) || throw(ArgumentError(
        "PSCAD Cable_Coax supports one to four concentric components",
    ))

    components = CableComponent[]
    previous_outer = 0.0
    for index in 1:component_count
        fields = _PSCAD_PART_FIELDS[index]
        conductor_inner = fields.inner === nothing ? previous_outer :
                          _pscad_number(values, fields.inner)
        conductor_outer = _pscad_number(values, fields.outer)
        insulation_outer = _pscad_number(values, fields.insulation)
        conductor = ConductorGroup(Tubular(
            conductor_inner,
            conductor_outer,
            _pscad_material(
                :conductor,
                _pscad_number(values, fields.rho),
                0.0,
                _pscad_number(values, fields.conductor_mu)
            )
        ))
        insulation = InsulatorGroup(
            Insulator(
                conductor_outer,
                insulation_outer,
                _pscad_dielectric(values, fields, frequency)
            );
            reference_frequency = frequency
        )
        name = strip(_pscad_parameter(values, "CONNAM$index"))
        isempty(name) && (name = "component$index")
        lowercase(name) == "none" && throw(ArgumentError(
            "PSCAD active component $index cannot be named 'none'",
        ))
        push!(components, CableComponent(lowercasefirst(name), conductor, insulation))
        previous_outer = insulation_outer
    end
    cable_name = strip(get(values, "Name", ""))
    isempty(cable_name) && (cable_name = "cable$cable_number")
    return CableDesign(cable_name, components)
end

const _PSCAD_SIMPLIFIED_FIELDS = (
    (
        inner = "R1", outer = "R2", insulation = "R3",
        rho = "RHOC", resistance = "DCRC", material_mode = "DTC",
        mu = "PERMC", eps = "EPS1", capacitance = "CI1",
        dielectric_mode = "DTI1", insulation_mu = "mu_r1", loss = "LT1",
        name = "core"
    ),
    (
        inner = "R3", outer = "R4", insulation = "R5",
        rho = "RHOS", resistance = "DCRS", material_mode = "DTS",
        mu = "PERMS", eps = "EPS2", capacitance = "CI2",
        dielectric_mode = "DTI2", insulation_mu = "mu_r2", loss = "LT2",
        name = "sheath"
    ),
    (
        inner = "R5", outer = "R6", insulation = "R7",
        rho = "RHOA", resistance = "DCRA", material_mode = "DTA",
        mu = "PERMA", eps = "EPS3", capacitance = "CI3",
        dielectric_mode = "DTI3", insulation_mu = "mu_r3", loss = "LT3",
        name = "armor"
    )
)

function _pscad_simplified_rho(values, fields, r_in, r_ex)
    mode = _pscad_integer(values, fields.material_mode)
    mode in (0, 1) || throw(ArgumentError(
        "PSCAD simplified-cable material selector must be zero or one",
    ))
    mode == 1 && return _pscad_number(values, fields.rho)
    resistance = _pscad_number(values, fields.resistance) / 1000
    return resistance * π * (r_ex^2 - r_in^2)
end

function _pscad_simplified_eps(values, fields, r_in, r_ex)
    mode = _pscad_integer(values, fields.dielectric_mode)
    mode in (0, 1) || throw(ArgumentError(
        "PSCAD simplified-cable dielectric selector must be zero or one",
    ))
    mode == 1 && return _pscad_number(values, fields.eps)
    capacitance = _pscad_number(values, fields.capacitance) * 1e-9
    return capacitance * log(r_ex / r_in) / (2π * _PSCAD_EPSILON_0)
end

function _pscad_simplified_design(values, cable_number::Int, frequency)
    _pscad_integer(values, "RorT") == 0 || throw(ArgumentError(
        "PSCAD simplified-cable import currently requires radius input",
    ))
    _pscad_integer(values, "SemiCL") == 0 || throw(ArgumentError(
        "PSCAD simplified-cable import does not support semiconductive layers",
    ))
    line_layers = _pscad_integer(values, "LL")
    isodd(line_layers) || throw(ArgumentError(
        "PSCAD LL must describe alternating conductor and insulation layers",
    ))
    component_count = (line_layers + 1) ÷ 2
    component_count in eachindex(_PSCAD_SIMPLIFIED_FIELDS) || throw(ArgumentError(
        "PSCAD Cable_CoaxSimpl supports one to three concentric components",
    ))
    components = CableComponent[]
    for index in 1:component_count
        fields = _PSCAD_SIMPLIFIED_FIELDS[index]
        conductor_inner = _pscad_number(values, fields.inner)
        conductor_outer = _pscad_number(values, fields.outer)
        insulation_outer = _pscad_number(values, fields.insulation)
        conductor = ConductorGroup(Tubular(
            conductor_inner,
            conductor_outer,
            _pscad_material(
                :conductor,
                _pscad_simplified_rho(
                    values, fields, conductor_inner, conductor_outer
                ),
                0.0,
                _pscad_number(values, fields.mu)
            )
        ))
        eps_r = _pscad_simplified_eps(
            values, fields, conductor_outer, insulation_outer
        )
        loss = _pscad_number(values, fields.loss)
        conductivity = 2π * frequency * _PSCAD_EPSILON_0 * eps_r * loss
        insulation = InsulatorGroup(
            Insulator(
                conductor_outer,
                insulation_outer,
                _pscad_material(
                    :insulator,
                    iszero(conductivity) ? Inf : inv(conductivity),
                    eps_r,
                    _pscad_number(values, fields.insulation_mu)
                )
            );
            reference_frequency = frequency
        )
        push!(components, CableComponent(fields.name, conductor, insulation))
    end
    cable_name = strip(get(values, "Name", ""))
    isempty(cable_name) && (cable_name = "cable$cable_number")
    return CableDesign(cable_name, components)
end

function _pscad_position(values, cable_number::Int, next_phase::Ref{Int})
    frequency = _pscad_number(values, "FLT")
    frequency > 0 || throw(DomainError(
        frequency, "PSCAD cable reference frequency must be positive"
    ))
    design = _pscad_design(values, cable_number, frequency)
    connections = Dict{String, Int}()
    for (index, component) in enumerate(design.components)
        eliminated = index > 1 && _pscad_integer(values, "elim$(index - 1)") != 0
        connections[component.id] = eliminated ? 0 : next_phase[]
        eliminated || (next_phase[] += 1)
    end

    horizontal = _pscad_number(values, "X")
    overhead = _pscad_integer(values, "OHC") != 0
    vertical = overhead ? _pscad_number(values, "Y2") :
               -abs(_pscad_number(values, "Y"))
    return CablePosition(design, horizontal, vertical, connections)
end

function _pscad_simplified_positions(values, next_phase::Ref{Int})
    circuit_count = _pscad_integer(values, "NC")
    circuit_count > 0 || throw(DomainError(
        circuit_count, "PSCAD simplified-cable circuit count must be positive"
    ))
    first_number = _pscad_integer(values, "CABNUM")
    spacing = _pscad_number(values, "D")
    horizontal = _pscad_number(values, "X")
    vertical = -abs(_pscad_number(values, "Y"))
    frequency = _pscad_number(values, "FLT")
    positions = CablePosition[]
    for offset in 0:(3circuit_count - 1)
        cable_number = first_number + offset
        design = _pscad_simplified_design(values, cable_number, frequency)
        connections = Dict(
            component.id => (next_phase[] += 1; next_phase[] - 1)
        for component in design.components)
        push!(positions, CablePosition(
            design,
            horizontal + offset * spacing,
            vertical,
            connections
        ))
    end
    return positions
end

function _pscad_earth(values)
    return EarthModel(
        _pscad_number(values, "GRRES"),
        _pscad_number(values, "GRP"),
        _pscad_number(values, "GPERM")
    )
end

"""
$(TYPEDSIGNATURES)

Import a self-contained PSCAD project that follows the supported
frequency-dependent coaxial-cable schema.

# Arguments

- `file_name`: Existing PSCAD project with a `.pscx` extension.

# Keywords

- `definition`: Qualified or local name of the frequency-dependent row to import. Leave it unset when the project contains one eligible row or one row selected for output.

# Returns

- `earth`: Materialised homogeneous [`EarthModel`](@ref).
- `system`: Materialised [`LineCableSystem`](@ref).

# Notes

The PSCAD schema stores equivalent concentric layers but does not retain the
source cable's detailed strand geometry, material reference temperature, or
temperature coefficient. Imported layers therefore use the emitted equivalent
properties at 20 °C with zero temperature coefficient. Non-eliminated PSCAD
conductors receive consecutive phase indices in cable and layer order.

PSCAD bounds dielectric loss tangent at ten. A project produced from a more
conductive equivalent dielectric therefore contains the bounded PSCAD value,
and import materialises that value rather than the pre-export resistivity.

The importer accepts detailed and simplified coaxial cables. Simplified cables
must use radius input and may not contain semiconductive layers. Tower and
pipe-type definitions are rejected.

# Errors

Throws `ArgumentError`, `DomainError`, `DimensionMismatch`, or `KeyError` when
the file, schema, geometry, or physical parameters are invalid.
"""
function import_data(
        ::Val{:pscad},
        file_name::AbstractString;
        definition::Union{Nothing, AbstractString} = nothing
)
    extension = lowercase(splitext(file_name)[2])
    extension == ".pscx" || throw(ArgumentError(
        "PSCAD import requires a .pscx file",
    ))
    isfile(file_name) || throw(ArgumentError(
        "PSCAD project does not exist: $(_display_path(file_name))",
    ))
    document = readxml(file_name)
    project = root(document)
    nodename(project) == "project" || throw(ArgumentError(
        "PSCAD input root must be <project>",
    ))
    haskey(project, "name") || throw(ArgumentError(
        "PSCAD project has no name",
    ))

    row = _pscad_row_definition(project, definition)
    cables = _pscad_active_cables(row)
    instance = _pscad_parameters(_pscad_instance(document, project, row))
    system_id = strip(_pscad_parameter(instance, "Name"))
    isempty(system_id) && throw(ArgumentError("PSCAD line instance has no name"))
    line_length = 1000 * _pscad_number(instance, "Length")

    next_phase = Ref(1)
    numbered_positions = Pair{Int, CablePosition}[]
    for cable in cables
        values = _pscad_parameters(cable)
        first_number = _pscad_integer(values, "CABNUM")
        if _pscad_binding(cable) == "master:Cable_CoaxSimpl"
            append!(numbered_positions,
                (first_number + offset - 1) => position
                for (offset, position) in enumerate(
                    _pscad_simplified_positions(values, next_phase)
                ))
        else
            push!(numbered_positions, first_number => _pscad_position(
                values,
                first_number,
                next_phase
            ))
        end
    end
    sort!(numbered_positions; by = first)
    first.(numbered_positions) == collect(eachindex(numbered_positions)) ||
        throw(ArgumentError(
            "PSCAD cable numbers must be consecutive and start at one",
        ))
    positions = last.(numbered_positions)
    earth = _pscad_earth(_pscad_parameters(row.ground))
    system = LineCableSystem(system_id, line_length, positions)
    return earth, system
end
