function _pscad_output_path(system::LineCableSystem, file_name)
    if file_name === nothing
        return joinpath(@__DIR__, "$(system.system_id)_export.pscx")
    end
    requested = isabspath(file_name) ? String(file_name) : abspath(file_name)
    directory, name = splitdir(requested)
    return joinpath(directory, "$(system.system_id)_$name")
end

function _pscad_part_parameters(component, index::Int, angular_frequency)
    names = ("C", "S", "A", "O")
    radii = (("R1", "R2", "R3"), ("R4", "R4", "R5"),
        ("R6", "R6", "R7"), ("R8", "R8", "R9"))
    resistivities = ("RHOC", "RHOS", "RHOA", "RHOO")
    permeabilities = ("PERMC", "PERMS", "PERMA", "PERMO")
    dielectric_permittivities = ("EPS1", "EPS2", "EPS3", "EPS4")
    dielectric_permeabilities = ("PERM1", "PERM2", "PERM3", "PERM4")
    loss_tangents = ("LT1", "LT2", "LT3", "LT4")
    conductor_name = "CONNAM$index"
    inner_name, outer_name, insulation_name = radii[index]
    capacitance = component.insulator_group.shunt_capacitance
    conductance = component.insulator_group.shunt_conductance
    loss = iszero(capacitance) ? zero(capacitance) :
           conductance / (angular_frequency * capacitance)

    parameters = Pair{String, String}[
        conductor_name => uppercasefirst(component.id),
        outer_name => _pscad_value(component.conductor_group.r_ex),
        resistivities[index] => _pscad_value(component.conductor_props.rho),
        permeabilities[index] => _pscad_value(component.conductor_props.mu_r),
        insulation_name => _pscad_value(component.insulator_group.r_ex),
        dielectric_permittivities[index] => _pscad_value(component.insulator_props.eps_r),
        dielectric_permeabilities[index] => _pscad_value(component.insulator_props.mu_r),
        loss_tangents[index] => _pscad_value(loss; maximum = 10)
    ]
    index == 1 && pushfirst!(
        parameters,
        inner_name => _pscad_value(component.conductor_group.r_in)
    )
    index > 1 && push!(parameters, "elim$(index - 1)" => "0")
    push!(parameters, "T$(2index + 1)" => "0.0000")
    index == 1 && append!(parameters, [
        "SemiCL" => "0", "SL2" => "0.0000", "SL1" => "0.0000"
    ])
    return parameters
end

function _empty_pscad_part(index::Int)
    radii = (("R1", "R2", "R3"), ("R4", "R4", "R5"),
        ("R6", "R6", "R7"), ("R8", "R8", "R9"))
    resistivities = ("RHOC", "RHOS", "RHOA", "RHOO")
    permeabilities = ("PERMC", "PERMS", "PERMA", "PERMO")
    epsilons = ("EPS1", "EPS2", "EPS3", "EPS4")
    mus = ("PERM1", "PERM2", "PERM3", "PERM4")
    losses = ("LT1", "LT2", "LT3", "LT4")
    _, outer_name, insulation_name = radii[index]
    parameters = Pair{String, String}[
        "CONNAM$index" => "none",
        outer_name => "0.0",
        resistivities[index] => "0.0",
        permeabilities[index] => "0.0",
        insulation_name => "0.0",
        epsilons[index] => "0.0",
        mus[index] => "0.0",
        losses[index] => "0.0000",
        "T$(2index + 1)" => "0.0000"
    ]
    index > 1 && push!(parameters, "elim$(index - 1)" => "0")
    return parameters
end

function _pscad_cable_parameters(position, index::Int, base_frequency)
    components = position.design_data.components
    length(components) <= 4 || throw(ArgumentError(
        "PSCAD Cable_Coax supports at most four concentric components",
    ))
    length(position.conn) == length(components) || throw(DimensionMismatch(
        "PSCAD phase mapping must match the cable component count",
    ))
    parameters = Pair{String, String}[
        "CABNUM" => string(index),
        "Name" => position.design_data.cable_id,
        "X" => _pscad_value(position.horz),
        "OHC" => position.vert < 0 ? "0" : "1",
        "Y" => position.vert < 0 ? _pscad_value(abs(position.vert)) : "0.0",
        "Y2" => position.vert > 0 ? _pscad_value(position.vert) : "0.0",
        "ShuntA" => "1.0e-11 [mho/m]",
        "FLT" => _pscad_value(base_frequency),
        "RorT" => "0",
        "LL" => string(2length(components) - 1),
        "CROSSBOND" => "0",
        "GROUPNO" => "1",
        "CBC1" => "1",
        "CBC2" => "0",
        "CBC3" => "0",
        "CBC4" => "0",
        "SHRad" => "1",
        "LC" => "3"
    ]
    angular_frequency = 2 * pi * base_frequency
    for (component_index, component) in enumerate(components)
        append!(parameters, _pscad_part_parameters(
            component, component_index, angular_frequency
        ))
        if component_index > 1
            elimination = "elim$(component_index - 1)"
            parameter_index = findlast(pair -> first(pair) == elimination, parameters)
            parameters[parameter_index] = elimination =>
                (position.conn[component_index] == 0 ? "1" : "0")
        end
    end
    for component_index in (length(components) + 1):4
        append!(parameters, _empty_pscad_part(component_index))
    end
    return parameters
end

function _pscad_project(system::LineCableSystem, earth::EarthModel, base_frequency)
    document = XMLDocument()
    project = ElementNode("project")
    setroot!(document, project)
    _pscad_attributes!(project, (
        name = system.system_id,
        version = "5.0.2",
        schema = "",
        Target = "EMTDC"
    ))
    settings = addelement!(project, "paramlist")
    settings["name"] = "Settings"
    for (name, value) in (
        "creator" => "LineCableModels.jl",
        "description" => "Generated by LineCableModels.jl",
        "time_duration" => "0.5",
        "time_step" => "5",
        "sample_step" => "250"
    )
        parameter = addelement!(settings, "param")
        parameter["name"] = name
        parameter["value"] = value
    end

    definitions = addelement!(project, "definitions")
    definition = addelement!(definitions, "Definition")
    _pscad_attributes!(definition,
        (
            id = 100000001,
            classid = "RowDefn",
            name = "CableSystem",
            version = "RowDefn",
            build = "RowDefn",
            crc = -1,
            view = "false"
        ))
    schematic = addelement!(definition, "schematic")
    schematic["classid"] = "RowCanvas"

    identifier = 100000010
    _pscad_user!(schematic, identifier, _PSCAD_FREQUENCY_BINDING,
        [
            "FS" => "0.5", "FE" => "1.0E6", "Numf" => "100",
            "DCCOR" => "1", "CPASS" => "0"
        ];
        x = 576, y = 180)
    identifier += 1
    for (index, position) in enumerate(system.cables)
        _pscad_user!(
            schematic,
            identifier,
            _PSCAD_CABLE_BINDING,
            _pscad_cable_parameters(position, index, base_frequency);
            x = 234 + (index - 1) * 400,
            y = 612
        )
        identifier += 1
    end
    ground = last(earth.layers)
    _pscad_user!(schematic,
        identifier,
        _PSCAD_GROUND_BINDING,
        [
            "EarthForm2" => "0",
            "EarthForm" => "3",
            "EarthForm3" => "2",
            "GrRho" => "0",
            "GRRES" => _pscad_value(ground.rho),
            "GPERM" => _pscad_value(ground.mu_r),
            "K0" => "0.001",
            "K1" => "0.01",
            "alpha" => "0.7",
            "GRP" => _pscad_value(ground.eps_r)
        ];
        x = 504,
        y = 288)
    return document
end
