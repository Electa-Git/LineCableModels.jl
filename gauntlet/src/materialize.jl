const DM = LineCableModels.DataModel
const EP = LineCableModels.EarthProps
const EN = LineCableModels.Engine
const _EPSILON_0 = 8.8541878128e-12

_family(::Val{:coax}) = Coax()
_family(::Val{:overhead}) = Overhead()
_family(::Val{:mixed}) = Mixed()
_family(::Val{:pipe}) = Pipe()
_family(input::PSCADIO.ParsedInput) = _family(Val(input.kind))

function _field(section::PSCADIO.Block, name::AbstractString)
    haskey(section.fields, name) || throw(KeyError(name))
    return section.fields[name]
end

function _field(section::PSCADIO.Block, name::AbstractString, default)
    get(section.fields, name, default)
end

function _material(rho, eps_r, mu_r)
    return LineCableModels.Material(rho, eps_r, mu_r, 20.0, 0.0)
end

function _dielectric(eps_r, mu_r, loss_tangent, frequency)
    omega = 2π * frequency
    conductivity = omega * _EPSILON_0 * eps_r * loss_tangent
    rho = iszero(conductivity) ? 1e30 : inv(conductivity)
    return _material(rho, eps_r, mu_r)
end

_air() = _material(1e30, 1.0, 1.0)
_shell(radius) = max(abs(radius) * 1e-6, 1e-9)

function _component(
        id::AbstractString,
        conductor_inner,
        conductor_outer,
        conductor_material,
        insulation_outer,
        insulation_material,
        reference_frequency
)
    conductor = DM.ConductorGroup(DM.Tubular(
        conductor_inner, conductor_outer, conductor_material
    ))
    insulation = DM.InsulatorGroup(
        DM.Insulator(conductor_outer, insulation_outer, insulation_material);
        reference_frequency
    )
    return DM.CableComponent(id, conductor, insulation)
end

function _coax_design(section::PSCADIO.Block, id::AbstractString)
    reference_frequency = Float64(_field(section, "Frequency for loss tangent", 60.0))
    conductor_inner = Float64(_field(section, "Conductor Inner Radius"))
    conductor_outer = Float64(_field(section, "Conductor Outer Radius"))
    conductor_material = _material(
        _field(section, "Conductor Resistivity"),
        0.0,
        _field(section, "Conductor Permeability")
    )
    layers = Int(_field(section, "Layers"))
    components = DM.CableComponent[]

    if layers == 0
        push!(components,
            _component(
                "conductor",
                conductor_inner,
                conductor_outer,
                conductor_material,
                conductor_outer + _shell(conductor_outer),
                _air(),
                reference_frequency
            ))
        return DM.CableDesign(id, components)
    end

    insulation_one_outer = Float64(_field(section, "Insulator 1 Outer Radius"))
    insulation_one = _dielectric(
        _field(section, "Insulator 1 Relative Permittivity"),
        _field(section, "Insulator 1 Relative Permeability"),
        _field(section, "Insulator 1 Loss Tangent", 0.0),
        reference_frequency
    )
    push!(components,
        _component(
            "core",
            conductor_inner,
            conductor_outer,
            conductor_material,
            insulation_one_outer,
            insulation_one,
            reference_frequency
        ))

    layers >= 2 || return DM.CableDesign(id, components)
    sheath_outer = Float64(_field(section, "Sheath Outer Radius"))
    sheath_material = _material(
        _field(section, "Sheath Resistivity"),
        0.0,
        _field(section, "Sheath Permeability")
    )
    if layers >= 3
        insulation_two_outer = Float64(_field(section, "Insulator 2 Outer Radius"))
        insulation_two = _dielectric(
            _field(section, "Insulator 2 Relative Permittivity"),
            _field(section, "Insulator 2 Relative Permeability"),
            _field(section, "Insulator 2 Loss Tangent", 0.0),
            reference_frequency
        )
    else
        insulation_two_outer = sheath_outer + _shell(sheath_outer)
        insulation_two = _air()
    end
    push!(components,
        _component(
            "sheath",
            insulation_one_outer,
            sheath_outer,
            sheath_material,
            insulation_two_outer,
            insulation_two,
            reference_frequency
        ))
    return DM.CableDesign(id, components)
end

function _position(section::PSCADIO.Block, family::Coax)
    values = Float64.(_field(section, "P1"))
    return values[1], -abs(values[2])
end

function _position(section::PSCADIO.Block, family::Mixed)
    values = Float64.(_field(section, "P1"))
    overhead = _field(section, "Overhead Cable", 0) == 1
    return values[1], overhead ? abs(values[2]) : -abs(values[2])
end

function _case_reduction(name)
    normalized = lowercase(String(name))
    normalized == "retained" && return Retained()
    normalized == "reduced" && return Reduced()
    return NoReduction()
end

function _connections(
        design,
        ::Reduction,
        next_phase::Base.RefValue{Int},
        section::Union{Nothing, PSCADIO.Block} = nothing
)
    result = Dict{String, Int}()
    last_index = length(design.components)
    ground_code = section === nothing ? 0 :
                  Int(_field(section, "Ground Last Layer", 0))
    ground_last = ground_code == 1 || ground_code == 2 && last_index > 1
    for (index, component) in enumerate(design.components)
        eliminated = ground_last && index == last_index
        if eliminated
            result[component.id] = 0
        else
            result[component.id] = next_phase[]
            next_phase[] += 1
        end
    end
    return result
end

function _ports(design, cable::Int, connections)
    return [Port(
                "cable$(cable):$(component.id)",
                cable,
                component.id,
                connections[component.id]
            ) for component in design.components]
end

function _frequency(input::PSCADIO.ParsedInput)
    options = only(PSCADIO.block(input.root, "Frequency Dep. (Phase) Model Options"))
    first_frequency = Float64(_field(options, "Curve Fitting Start Frequency"))
    last_frequency = Float64(_field(options, "Curve Fitting End Frequency"))
    increments = Int(_field(options, "Total Number of Frequency Increments"))
    values = collect(10.0 .^ range(
        log10(first_frequency), log10(last_frequency); length = increments + 1
    ))
    for name in ("Cable Summary", "Line Summary")
        sections = PSCADIO.block(input.root, name)
        isempty(sections) && continue
        steady = _field(only(sections), "Steady State Frequency", nothing)
        steady === nothing || push!(values, Float64(steady))
    end
    return sort!(unique!(values))
end

function _earth(input::PSCADIO.ParsedInput)
    ground = only(PSCADIO.block(input.root, "Line Constants Ground Data"))
    return EP.EarthModel(
        _field(ground, "GroundResistivity"),
        _field(ground, "GroundPermittivity"),
        _field(ground, "GroundPermeability")
    )
end

function _line_length(input::PSCADIO.ParsedInput)
    names = ("Cable Summary", "Line Summary")
    for name in names
        matches = PSCADIO.block(input.root, name)
        isempty(matches) || return 1000 * Float64(_field(
            only(matches), name == "Cable Summary" ? "Cable Length" : "Line Length"
        ))
    end
    throw(KeyError("line length"))
end

function _formulation(input::PSCADIO.ParsedInput, reduction::Reduction)
    return _formulation(_family(input), input, reduction)
end

function _needs_kron(input::PSCADIO.ParsedInput, reduction::Reduction)
    for section in PSCADIO.block(input.root, "Coax Cable")
        layers = Int(_field(section, "Layers"))
        code = Int(_field(section, "Ground Last Layer", 0))
        (code == 1 || code == 2 && layers > 1) && return true
    end
    for tower in PSCADIO.block(input.root, "Line Constants Tower")
        any(wires -> _field(wires, "Eliminate Ground Wires", 1) == 1,
            PSCADIO.block(tower, "GroundWires")) && return true
    end
    for pipe in PSCADIO.block(input.root, "Pipe-Type Cable")
        _field(pipe, "Ground Pipe", 0) != 0 && return true
    end
    return false
end

function _formulation(
        family::Family,
        input::PSCADIO.ParsedInput,
        reduction::Reduction
)
    transposed = false
    towers = PSCADIO.block(input.root, "Line Constants Tower")
    for tower in towers, circuit in PSCADIO.block(tower, "Circuit")

        transposed |= _field(circuit, "Transposed", 0) == 1
    end
    for cable in PSCADIO.block(input.root, "Coax Cable")
        transposed |= _field(cable, "Conductor 1 Transposed", 0) == 1
    end
    earth_impedance = family isa Overhead ?
                      EN.EarthImpedance.Papadopoulos(s = 1, t = 1, Γx = 0) :
                      EN.EarthImpedance.Papadopoulos()
    # The native Papadopoulos admittance kernel is the broadband underground
    # formulation.  Applied to overhead conductors at millihertz frequencies,
    # its improper integral is both the wrong approximation for this campaign
    # and pathologically expensive.  PSCAD's overhead shunt path is represented
    # by the native electrostatic-images formulation instead.
    earth_admittance = family isa Overhead ?
                       EN.EarthAdmittance.Images() :
                       EN.EarthAdmittance.Pollaczek()
    return LineCableModels.Formulation(
        :EMT;
        earth_impedance,
        earth_admittance,
        insulation_admittance = EN.InsulationAdmittance.ParallelRC(),
        options = (
            reduce_bundle = false,
            kron_reduction = _needs_kron(input, reduction),
            ideal_transposition = transposed,
            temperature_correction = false
        )
    )
end

function _assumptions(input::PSCADIO.ParsedInput, family::Family)
    assumptions = Assumption[]
    ground = only(PSCADIO.block(input.root, "Line Constants Ground Data"))
    earth_code = Int(_field(ground, "EarthImpedanceFormula", 0))
    push!(assumptions,
        Assumption(
            :earth,
            "PSCAD earth code $earth_code is evaluated with the closest available " *
            "LineCableModels Papadopoulos impedance/admittance formulation."
        ))
    family isa Overhead && push!(assumptions,
        Assumption(
            :geometry,
            "Tower sag is represented by its two-thirds mean-height correction and " *
            "a vanishing air shell required by the native cable-shaped input."
        ))
    family isa Overhead && push!(assumptions,
        Assumption(
            :admittance,
            "PSCAD's overhead shunt path is represented by the native electrostatic-" *
            "images earth-admittance formulation."
        ))
    !(family isa Overhead) && push!(assumptions,
        Assumption(
            :admittance,
            "PSCAD's underground shunt path is represented by the native Pollaczek " *
            "earth-admittance formulation without displacement currents."
        ))
    family isa Mixed && push!(assumptions,
        Assumption(
            :formulation,
            "The native solver has no cross-interface earth kernel; overhead members " *
            "are mirrored below the interface so the complete mixed case remains computable."
        ))
    family isa Pipe && push!(assumptions,
        Assumption(
            :pipe,
            "The common pipe is represented by the closest independent-coax geometry."
        ))
    any(section -> Int(_field(section, "Layers")) < 3,
        PSCADIO.block(input.root, "Coax Cable")) && push!(assumptions,
        Assumption(
            :insulation,
            "A vanishing air shell closes PSCAD conductors without an explicit outer dielectric."
        ))
    return assumptions
end

function _coax_problem(
        input::PSCADIO.ParsedInput,
        family::Union{Coax, Mixed},
        reduction::Reduction,
        id::AbstractString
)
    sections = PSCADIO.block(input.root, "Coax Cable")
    isempty(sections) && throw(ArgumentError("PSCAD cable input has no coax sections"))
    phase = Ref(1)
    ports = Port[]
    system = nothing
    for (cable_index, section) in enumerate(sections)
        design = _coax_design(section, "$id-cable-$cable_index")
        connections = _connections(design, reduction, phase, section)
        append!(ports, _ports(design, cable_index, connections))
        horizontal, vertical = _position(section, family)
        position = DM.CablePosition(design, horizontal, vertical, connections)
        if system === nothing
            system = DM.LineCableSystem(id, _line_length(input), position)
        else
            DM.add!(system, position)
        end
    end
    problem = EN.LineParametersProblem(
        system;
        temperature = 20.0,
        earth_props = _earth(input),
        frequencies = _frequency(input)
    )
    return problem, ports
end

function _tower_conductors(input::PSCADIO.ParsedInput)
    records = NamedTuple[]
    for (tower_index, tower) in enumerate(
        PSCADIO.block(input.root, "Line Constants Tower")
    )
        for (circuit_index, circuit) in enumerate(PSCADIO.block(tower, "Circuit"))
            count = Int(_field(circuit, "Conductors"))
            radius = Float64(_field(circuit, "Radius"))
            resistance = Float64(_field(circuit, "DCResistance")) / 1000
            permeability = Float64(_field(circuit, "Conductor Relative Permeability"))
            sag = Float64(_field(circuit, "Sag", 0.0))
            phases = Int.(_field(circuit, "Conductor Phase Information"))
            for index in 1:count
                point = Float64.(_field(circuit, "P$index"))
                push!(records,
                    (;
                        radius,
                        resistance,
                        permeability,
                        horizontal = point[1],
                        vertical = point[2] - 2sag / 3,
                        phase = phases[index],
                        grounded = false,
                        group = "circuit:$tower_index:$circuit_index"
                    ))
            end
        end
        for (group_index, wires) in enumerate(PSCADIO.block(tower, "GroundWires"))
            count = Int(something(wires.value))
            radius = Float64(_field(wires, "Radius"))
            resistance = Float64(_field(wires, "DCResistance")) / 1000
            permeability = Float64(_field(wires, "Ground Wire Relative Permeability"))
            sag = Float64(_field(wires, "Sag", 0.0))
            eliminated = _field(wires, "Eliminate Ground Wires", 1) == 1
            phases = Int.(_field(wires, "Ground Wire Phase Information", fill(0, count)))
            for index in 1:count
                point = Float64.(_field(wires, "P$index"))
                push!(records,
                    (;
                        radius,
                        resistance,
                        permeability,
                        horizontal = point[1],
                        vertical = point[2] - 2sag / 3,
                        phase = eliminated ? 0 : phases[index],
                        grounded = eliminated,
                        group = "ground:$tower_index:$group_index:$index"
                    ))
            end
        end
    end
    return records
end

function _overhead_problem(
        input::PSCADIO.ParsedInput,
        reduction::Reduction,
        id::AbstractString
)
    conductors = _tower_conductors(input)
    isempty(conductors) && throw(ArgumentError("PSCAD tower input has no conductors"))
    system = nothing
    ports = Port[]
    phase = 0
    for (index, conductor) in enumerate(conductors)
        rho = conductor.resistance * π * conductor.radius^2
        material = _material(rho, 0.0, conductor.permeability)
        design = DM.CableDesign(
            "$id-conductor-$index",
            _component(
                "conductor",
                0.0,
                conductor.radius,
                material,
                conductor.radius + _shell(conductor.radius),
                _air(),
                60.0
            )
        )
        assigned = conductor.grounded || reduction isa Reduced && conductor.phase == 0 ?
                   0 : (phase += 1)
        connections = Dict("conductor" => assigned)
        push!(ports, Port(
            "$(conductor.group):$index", index, conductor.group, assigned
        ))
        position = DM.CablePosition(
            design, conductor.horizontal, abs(conductor.vertical), connections
        )
        if system === nothing
            system = DM.LineCableSystem(id, _line_length(input), position)
        else
            DM.add!(system, position)
        end
    end
    problem = EN.LineParametersProblem(
        system;
        temperature = 20.0,
        earth_props = _earth(input),
        frequencies = _frequency(input)
    )
    return problem, ports
end

function _mixed_problem(
        input::PSCADIO.ParsedInput,
        reduction::Reduction,
        id::AbstractString
)
    conductors = _tower_conductors(input)
    sections = PSCADIO.block(input.root, "Coax Cable")
    isempty(conductors) && isempty(sections) &&
        throw(ArgumentError(
            "PSCAD mixed input has neither tower nor cable conductors",
        ))
    system = nothing
    ports = Port[]
    next_phase = Ref(1)
    cable_index = 0

    for (index, conductor) in enumerate(conductors)
        cable_index += 1
        rho = conductor.resistance * π * conductor.radius^2
        material = _material(rho, 0.0, conductor.permeability)
        design = DM.CableDesign(
            "$id-overhead-$index",
            _component(
                "conductor",
                0.0,
                conductor.radius,
                material,
                conductor.radius + _shell(conductor.radius),
                _air(),
                60.0
            )
        )
        assigned = conductor.grounded ? 0 : next_phase[]
        conductor.grounded || (next_phase[] += 1)
        connections = Dict("conductor" => assigned)
        push!(ports, Port(
            "$(conductor.group):$index", cable_index, conductor.group, assigned
        ))
        position = DM.CablePosition(
            design, conductor.horizontal, -abs(conductor.vertical), connections
        )
        if system === nothing
            system = DM.LineCableSystem(id, _line_length(input), position)
        else
            DM.add!(system, position)
        end
    end

    for (index, section) in enumerate(sections)
        cable_index += 1
        design = _coax_design(section, "$id-cable-$index")
        connections = _connections(design, reduction, next_phase, section)
        append!(ports, _ports(design, cable_index, connections))
        horizontal, vertical = _position(section, Mixed())
        position = DM.CablePosition(
            design, horizontal, -abs(vertical), connections
        )
        if system === nothing
            system = DM.LineCableSystem(id, _line_length(input), position)
        else
            DM.add!(system, position)
        end
    end

    problem = EN.LineParametersProblem(
        system;
        temperature = 20.0,
        earth_props = _earth(input),
        frequencies = _frequency(input)
    )
    return problem, ports
end

function materialize(
        input::PSCADIO.ParsedInput,
        reduction::Reduction,
        id::AbstractString
)
    return materialize(_family(input), input, reduction, id)
end

materialize(::Coax, input, reduction, id) = _coax_problem(input, Coax(), reduction, id)
materialize(::Mixed, input, reduction, id) = _mixed_problem(input, reduction, id)
materialize(::Overhead, input, reduction, id) = _overhead_problem(input, reduction, id)

function materialize(::Pipe, input, reduction, id)
    pipe = only(PSCADIO.block(input.root, "Pipe-Type Cable"))
    synthetic = PSCADIO.ParsedInput(:coax, PSCADIO.Block("root"))
    for name in ("Cable Summary", "Frequency Dep. (Phase) Model Options",
        "Line Constants Ground Data")
        append!(synthetic.root.children, PSCADIO.block(input.root, name))
    end
    pipe_position = Float64.(_field(pipe, "P1"))
    for child in PSCADIO.block(pipe, "Coax Cable")
        copied = deepcopy(child)
        local_position = Float64.(_field(copied, "P1"))
        copied.fields["P1"] = [
            pipe_position[1] + local_position[1],
            abs(pipe_position[2]) + local_position[2]
        ]
        push!(synthetic.root.children, copied)
    end
    pipe_outer = Float64(_field(pipe, "Outer Insulator Outer Radius"))
    independent_pipe = PSCADIO.Block("Coax Cable")
    merge!(independent_pipe.fields,
        Dict{String, Any}(
            "Cable Number" => Int(_field(pipe, "Number of Inner Cables")) + 1,
            # The native model has no common-enclosure geometry. Keeping the pipe
            # as a separate, nonoverlapping tubular conductor preserves its port
            # and material data while the case remains explicitly approximate.
            "P1" => [pipe_position[1] + 10pipe_outer, abs(pipe_position[2])],
            "Layers" => 1,
            "Ground Last Layer" => Int(_field(pipe, "Ground Pipe", 0)),
            "Conductor Inner Radius" => Float64(_field(
                pipe, "Inner Insulator Outer Radius"
            )),
            "Conductor Outer Radius" => Float64(_field(
                pipe, "Conductor Outer Radius"
            )),
            "Conductor Resistivity" => Float64(_field(pipe, "Conductor Resistivity")),
            "Conductor Permeability" => Float64(_field(pipe, "Conductor Permeability")),
            "Insulator 1 Outer Radius" => pipe_outer,
            "Insulator 1 Relative Permittivity" => Float64(_field(
                pipe, "Outer Ins. Relative Permittivity"
            )),
            "Insulator 1 Relative Permeability" => Float64(_field(
                pipe, "Outer Ins. Relative Permeability"
            )),
            "Insulator 1 Loss Tangent" => Float64(_field(
                pipe, "Outer Ins. Loss Tangent", 0.0
            )),
            "Frequency for loss tangent" => Float64(_field(
                pipe, "Frequency for loss tangent", 60.0
            ))
        ))
    push!(synthetic.root.children, independent_pipe)
    return _coax_problem(synthetic, Coax(), reduction, id)
end

function _fidelity(family::Family, input::PSCADIO.ParsedInput)
    # Promotion to Exact is evidence-driven and frozen only after geometry,
    # reduction, formulation, and held-out numerical equivalence are verified.
    return Approximate()
end
