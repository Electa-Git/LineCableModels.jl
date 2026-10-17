function validate(problem::TransientAnalysisProblem{<:LineParametersProblem},
        ::PSCADTransientFormulation)
    line = problem.line
    validate(line, PSCADFormulation)
    validate(line.earth_props, PSCADFormulation)
    system = line.system
    length(system.designs) <= 12 ||
        throw(ArgumentError("PSCAD cable interfaces support at most twelve cables"))
    active = filter(>(0), system.connection_order)
    allunique(active) || throw(ArgumentError("PSCAD transient export does not support bundled phases"))
    for blueprint in _pscad_blueprints(system)
        length(blueprint.assembly_ranges) == 1 ||
            throw(ArgumentError("PSCAD transient export requires one coaxial assembly per cable"))
        length(only(blueprint.assembly_ranges)) <= 4 ||
            throw(ArgumentError("PSCAD cables support at most four conductor layers"))
    end
    return problem
end

# Parameter table: PowerImpedance.jl `src/PSCAD/emitters/cables.jl:38-49`.
# Native `local/pscad-exchange/m7-master.pslx:14322-14330` defines C1T:
# 1 exposes the core, 2 adds the sheath, 3 adds the armour, and 4 adds the outer conductor.
# Its external pins (`m7-master.pslx:14519-14528`) place the core at (-18, 0),
# sheath at (18, 18), armour at (36, 18), outer at (54, 18), repeated every 54 px.
function _pscad_transient_interface!(canvas, system, name, end_name, labels)
    canvas.next_cell = 6cld(canvas.next_cell, 6)
    x, y = allocate_cell!(canvas)
    counts = [count(t -> t.cable == cable, system.terminal_order)
              for cable in eachindex(system.designs)]
    parameters = ["Name" => name, "NCAB" => string(length(counts))]
    append!(parameters, ["C$(i)T" => string(i <= length(counts) ? counts[i] : 1) for i in 1:12])
    append!(parameters, ["pipe" => "0", "dim_s" => "0", "dim_r" => "0",
        "GroundS" => "0", "GroundA" => "0", "GroundOC" => "0",
        "ResS" => "0.001 [ohm]", "ResA" => "0.001 [ohm]", "ResOC" => "0.001 [ohm]",
        "send_recv" => end_name == "sending" ? "1" : "2",
        "recv_default_type" => "0", "recv_default" => "0"])
    component!(canvas, "cable_interface", (x + 144, y + 18), parameters)
    slots = zeros(Int, length(counts))
    for (terminal, label) in zip(system.terminal_order, labels)
        cable = terminal.cable
        slot = (slots[cable] += 1)
        offset = slot == 1 ? (-18, 54(cable - 1)) : (18(slot - 1), 54(cable - 1) + 18)
        port = (x + 144, y + 18) .+ offset
        direction = slot == 1 ? (-36, 0) : (0, 18(slot - 1))
        _pscad_terminal_stub!(canvas, port, label; offset = direction)
    end
    canvas.next_cell += 6cld(54length(counts) + 72, 216) - 1
    return nothing
end

function _pscad_terminal_stub!(canvas, port, label; offset = (0, 36))
    endpoint = port .+ offset
    if label === nothing
        component!(canvas, "ground", endpoint)
    else
        node_label!(canvas, label, endpoint)
    end
    wire!(canvas, port, endpoint)
    return endpoint
end

# Signal wiring: PowerImpedance.jl `src/PSCAD/blocks.jl:237-248`.
function _pscad_signal_stub!(canvas, port, name; offset = (0, 36))
    endpoint = port .+ offset
    parameters = ["Name" => name]
    component!(canvas, "datalabel", endpoint, parameters)
    wire!(canvas, port, endpoint)
    return endpoint
end

# Meter ports and output table: PowerImpedance.jl `src/PSCAD/emitters/probes.jl:25-35,63-68`.
function _pscad_transient_meters!(canvas, node, title)
    x, y = allocate_cell!(canvas)
    component!(canvas, "ammeter", (x + 72, y + 54), ["Name" => node * "_i"])
    _pscad_terminal_stub!(canvas, (x + 72, y + 54), node; offset = (-36, 0))
    _pscad_terminal_stub!(canvas, (x + 108, y + 54), node * "_branch"; offset = (36, 0))
    component!(canvas, "voltmetergnd", (x + 72, y + 144), ["Name" => node * "_v"])
    _pscad_terminal_stub!(canvas, (x + 72, y + 144), node; offset = (-36, 0))
    for (kind, row, unit) in (("v", y + 54, "V"), ("i", y + 144, "A"))
        component!(canvas, "pgb", (x + 288, row),
            ["Name" => "$title:$kind", "Units" => unit, "Scale" => "1000",
                "UseSignalName" => "0", "enab" => "1"])
        _pscad_signal_stub!(canvas, (x + 288, row), node * "_$kind"; offset = (-36, 0))
    end
    return nothing
end

function _pscad_transient_end!(canvas, assignments::Tuple, problem, end_name, formulation)
    system = problem.line.system
    labels = ["$(end_name)_c$(t.cable)_$(t.terminal)" for t in system.terminal_order]
    _pscad_transient_interface!(canvas, system, system.system_id, end_name, labels)
    for assignment in assignments
        index = only(findall(t -> "cable:$(t.cable):$(t.terminal)" == assignment.terminal,
            system.terminal_order))
        node = labels[index]
        _pscad_transient_meters!(canvas, node, "$end_name:$(assignment.terminal)")
        _pscad_transient_assignment!(canvas, assignment, node * "_branch", problem.window.step)
    end
    return nothing
end

function _pscad_transient_end!(canvas, ::Matched, problem, end_name, formulation)
    system = problem.line.system
    labels = ["$(end_name)_c$(t.cable)_$(t.terminal)" for t in system.terminal_order]
    _pscad_transient_interface!(canvas, system, system.system_id, end_name, labels)
    name = end_name * "_extension"
    parameters = _pscad_main_instance_parameters(system, formulation.frequency)
    c = 1 / sqrt(vacuum_permeability(Float64) * vacuum_permittivity(Float64))
    length_km = 1.1 * c * problem.window.tmax / 2 / 1000
    for (index, (field, value)) in pairs(parameters)
        field == "Name" && (parameters[index] = field => name)
        field == "Length" && (parameters[index] = field => _pscad_value(length_km; sigdigits = 17))
        field == "linc" && (parameters[index] = field => "$length_km [km]")
    end
    x, y = allocate_cell!(canvas)
    binding = "$(system.system_id):CableSystem"
    _pscad_wire!(canvas.schematic, canvas.next_id, canvas.next_id + 1,
        "Cable", binding, binding, binding, parameters;
        x, y, width = 107, height = 128, crc = -1)
    hierarchy = only(findall("/project/hierarchy/call/call", canvas.schematic))
    call = addelement!(hierarchy, "call")
    _pscad_attributes!(call, (link = canvas.next_id + 1, name = binding,
        z = -1, view = "false", instance = 0))
    canvas.next_id += 2
    _pscad_transient_interface!(canvas, system, name, "sending", labels)
    # The extension's far end has no assignments or measurement branches.
    _pscad_transient_interface!(canvas, system, name, "receiving",
        [label * "_open" for label in labels])
    return nothing
end

function _pscad_transient_assignment!(canvas, assignment::TerminalImpedance, node, step)
    x, y = allocate_cell!(canvas)
    if iszero(assignment.R) && iszero(assignment.L)
        node_label!(canvas, node, (x + 72, y + 54))
        component!(canvas, "ground", (x + 108, y + 54))
        wire!(canvas, (x + 72, y + 54), (x + 108, y + 54))
    else
        junction = !iszero(assignment.R) && !iszero(assignment.L) ? node * "_rl" : nothing
        if !iszero(assignment.R)
            component!(canvas, "resistor", (x + 72, y + 54), ["R" => _pscad_value(assignment.R)])
            _pscad_terminal_stub!(canvas, (x + 72, y + 54), node; offset = (-36, 0))
            _pscad_terminal_stub!(canvas, (x + 108, y + 54), junction; offset = (36, 0))
        end
        if !iszero(assignment.L)
            component!(canvas, "inductor", (x + 252, y + 54), ["L" => _pscad_value(assignment.L)])
            _pscad_terminal_stub!(canvas, (x + 252, y + 54), something(junction, node); offset = (-36, 0))
            _pscad_terminal_stub!(canvas, (x + 288, y + 54), nothing; offset = (36, 0))
        end
    end
    return nothing
end

# Source parameters and polarity: PowerImpedance.jl `src/PSCAD/emitters/sources.jl:75-96`.
function _pscad_transient_voltage_source!(canvas, node, amplitude, parameters)
    x, y = allocate_cell!(canvas)
    point = (x + 144, y + 54)
    component!(canvas, "source_1", point,
        [["Name" => node, "Type" => "6", "Grnd" => "0", "Spec" => "0", "Tc" => "0"];
            parameters])
    positive, negative = amplitude < 0 ? (nothing, node) : (node, nothing)
    _pscad_terminal_stub!(canvas, point .+ (-36, 0), positive; offset = (-36, 0))
    _pscad_terminal_stub!(canvas, point .+ (36, 0), negative; offset = (36, 0))
    return point
end

function _pscad_transient_assignment!(canvas, source::Step, node, step)
    # Timed breaker: PowerImpedance.jl `src/PSCAD/emitters/dc_precharge.jl:83-92`.
    # NAME=1 opens and NAME=0 closes. The DC source has no native start-time field.
    _pscad_transient_voltage_source!(canvas, node * "_dc", source.amplitude,
        ["AC" => "0", "Cntrl" => "0", "Vm" => _pscad_value(abs(source.amplitude) / 1000; sigdigits = 17)])
    x, y = allocate_cell!(canvas)
    command = node * "_switch"
    component!(canvas, "breaker1", (x + 144, y + 54),
        ["NAME" => command, "BOpen" => "1", "ENAB" => "0", "OPCUR" => "1",
            "RON" => "0", "ROFF" => "1.0E100", "CLVL" => "0", "ViewB" => "0"])
    _pscad_terminal_stub!(canvas, (x + 108, y + 54), node * "_dc"; offset = (-36, 0))
    _pscad_terminal_stub!(canvas, (x + 180, y + 54), node; offset = (36, 0))
    component!(canvas, "time-sig", (x + 36, y + 144))
    _pscad_signal_stub!(canvas, (x + 72, y + 144), node * "_time"; offset = (0, 36))
    component!(canvas, "compare", (x + 252, y + 144),
        ["X" => _pscad_value(source.t0; sigdigits = 17), "OL" => "1", "OH" => "0",
            "INTR" => "0", "Dim" => "1", "conv" => "0"])
    _pscad_signal_stub!(canvas, (x + 216, y + 144), node * "_time"; offset = (0, 36))
    _pscad_signal_stub!(canvas, (x + 288, y + 144), command; offset = (0, 36))
    return nothing
end

function _pscad_transient_assignment!(canvas, source::Cosine, node, step)
    # `m7-master.pslx:298323-298341`: internal AC magnitude uses L-G RMS kV.
    # Phase uses degrees. The sine convention needs +90 degrees for a cosine.
    phase = mod(rad2deg(source.phase) + 90 + 180, 360) - 180
    _pscad_transient_voltage_source!(canvas, node, source.amplitude,
        ["AC" => "1", "Cntrl" => "0",
            "Vm" => _pscad_value(abs(source.amplitude) / (1000sqrt(2)); sigdigits = 17),
            "f" => _pscad_value(source.frequency; sigdigits = 17),
            "Ph" => _pscad_value(phase; sigdigits = 17)])
    return nothing
end

function _pscad_transient_assignment!(canvas, source::DoubleExponential, node, step)
    # `m7-master.pslx:309472-309631`: surge1 supports custom alpha and beta in 1/us.
    # Its Ipk normalizes the peak. Use the unnormalized waveform's peak in kV
    # as a numeric signal for source_1's external instantaneous DC voltage input.
    point = _pscad_transient_voltage_source!(canvas, node, source.amplitude,
        ["AC" => "0", "Cntrl" => "1"])
    _pscad_signal_stub!(canvas, point .+ (0, 36), node * "_surge")
    x, y = allocate_cell!(canvas)
    peak_time = log(source.beta / source.alpha) / (source.beta - source.alpha)
    peak = exp(-source.alpha * peak_time) - exp(-source.beta * peak_time)
    component!(canvas, "surge1", (x + 144, y + 54),
        ["Name" => node * "_surge", "STD" => "0", "SW" => "-1", "Tstart" => "0",
            "ALPHA" => _pscad_value(source.alpha * 1e-6; sigdigits = 17),
            "BETA" => _pscad_value(source.beta * 1e-6; sigdigits = 17),
            "Ipk" => _pscad_value(abs(source.amplitude) * peak / 1000; sigdigits = 17)])
    _pscad_signal_stub!(canvas, (x + 180, y + 54), node * "_surge")
    return nothing
end
