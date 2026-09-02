const APPLICATION_CASE_ID = "ohl_ugc_transition_v1"
const APPLICATION_DEFAULT_SHARE = 0.5
const APPLICATION_DEFAULT_LENGTH_M = 100.0e3
const APPLICATION_VOLTAGE_KV = 380 / sqrt(3)
const APPLICATION_POWER_MW = 100.0

const APPLICATION_CONNECTIONS = (
    (node=:B3d, element=:tl1, side=2, terminal=1),
    (node=:B3d, element=:c1, side=2, terminal=1),
    (node=:B3q, element=:tl1, side=2, terminal=2),
    (node=:B3q, element=:c1, side=2, terminal=2),
    (node=:B2d, element=:g4, side=1, terminal=1),
    (node=:B2d, element=:tl1, side=1, terminal=1),
    (node=:B2q, element=:g4, side=1, terminal=2),
    (node=:B2q, element=:tl1, side=1, terminal=2),
    (node=:B4, element=:ugc, side=1, terminal=1),
    (node=:B4, element=:c1, side=1, terminal=1),
    (node=:BX, element=:ugc, side=2, terminal=1),
    (node=:BX, element=:ohl, side=1, terminal=1),
    (node=:B5, element=:ohl, side=2, terminal=1),
    (node=:B5, element=:c2, side=1, terminal=1),
    (node=:B6d, element=:tl78, side=1, terminal=1),
    (node=:B6d, element=:c2, side=2, terminal=1),
    (node=:B6q, element=:tl78, side=1, terminal=2),
    (node=:B6q, element=:c2, side=2, terminal=2),
    (node=:B7d, element=:g1, side=1, terminal=1),
    (node=:B7d, element=:tl78, side=2, terminal=1),
    (node=:B7q, element=:g1, side=1, terminal=2),
    (node=:B7q, element=:tl78, side=2, terminal=2),
)

const APPLICATION_BUILDER_OPTIONS = (
    voltageBase=APPLICATION_VOLTAGE_KV,
    power_flow=(is_bounded=(bus_voltage=true,),),
)

struct PreparedPowerFlow
    specification_hash::String
    network::Any
    network_model::Any
    powerflow::Any
end

function terminal_line(length, earth_resistivity)
    return PowerImpedance.overhead_line(
        length=length,
        conductors=PowerImpedance.Conductors(
            organization=:flat,
            nᵇ=3,
            nˢᵇ=1,
            Rᵈᶜ=0.063,
            rᶜ=0.015,
            yᵇᶜ=30,
            Δyᵇᶜ=0,
            Δxᵇᶜ=10,
            Δ̃xᵇᶜ=0,
            dˢᵇ=0,
            dˢᵃᵍ=10
        ),
        groundwires=PowerImpedance.Groundwires(
            nᵍ=2,
            Rᵍᵈᶜ=0.92,
            rᵍ=0.0062,
            Δxᵍ=6.5,
            Δyᵍ=7.5,
            dᵍˢᵃᵍ=10
        ),
        earth_parameters=(1, 1, earth_resistivity),
        transformation=true
    )
end

function application_converter_one()
    return only(PowerImpedance.mmc(
        PowerImpedance.NetworkBuilder.Grid;
        Vᵈᶜ=640,
        vDCbase=640,
        Vₘ=APPLICATION_VOLTAGE_KV,
        P_max=1500,
        P_min=-1500,
        P=-APPLICATION_POWER_MW,
        Q=100,
        Q_max=500,
        Q_min=-500,
        occ=PowerImpedance.PI_control(Kₚ=0.7691, Kᵢ=522.7654),
        ccc=PowerImpedance.PI_control(Kₚ=0.1048, Kᵢ=48.1914),
        pll=PowerImpedance.PI_control(Kₚ=0.28, Kᵢ=12.5664),
        q=PowerImpedance.PI_control(Kₚ=0.1, Kᵢ=31.4159),
        dc=PowerImpedance.PI_control(Kₚ=6, Kᵢ=15),
        timeDelay=200.0e-6,
        padeOrderNum=5,
        padeOrderDen=5
    ))
end

function application_converter_two()
    return only(PowerImpedance.mmc(
        PowerImpedance.NetworkBuilder.Grid;
        Vᵈᶜ=640,
        vDCbase=640,
        Vₘ=APPLICATION_VOLTAGE_KV,
        P_max=1000,
        P_min=-1000,
        P=APPLICATION_POWER_MW,
        Q=100,
        Q_max=1000,
        Q_min=-1000,
        vACbase_LL_RMS=333,
        turnsRatio=333 / 380,
        Lᵣ=0.0461,
        Rᵣ=0.4103,
        Lₐᵣₘ=30.0e-3,
        occ=PowerImpedance.PI_control(Kₚ=0.7691, Kᵢ=522.7654),
        ccc=PowerImpedance.PI_control(Kₚ=0.1048, Kᵢ=48.1914),
        pll=PowerImpedance.PI_control(Kₚ=0.28, Kᵢ=12.5664),
        p=PowerImpedance.PI_control(Kₚ=0.1, Kᵢ=31.4159),
        q=PowerImpedance.PI_control(Kₚ=0.1, Kᵢ=31.4159),
        timeDelay=200.0e-6,
        padeOrderNum=5,
        padeOrderDen=5
    ))
end

function corridor_elements(share, earth_resistivity, length_m)
    regularized = clamp(share, 1e-3, 1 - 1e-3)
    ohl = PowerImpedance.overhead_line(
        length=length_m * (1 - regularized),
        conductors=PowerImpedance.Conductors(
            organization=:flat,
            nᵇ=2,
            nˢᵇ=1,
            Rᵈᶜ=0.0266,
            rᶜ=44.8e-3 / 2,
            yᵇᶜ=18.0,
            Δyᵇᶜ=0.0,
            Δxᵇᶜ=7.3,
            Δ̃xᵇᶜ=0.0,
            dˢᵇ=0.0,
            dˢᵃᵍ=6.0
        ),
        groundwires=PowerImpedance.Groundwires(
            nᵍ=2,
            Rᵍᵈᶜ=0.92,
            rᵍ=0.0062,
            Δxᵍ=7.3,
            Δyᵍ=7.0,
            dᵍˢᵃᵍ=6.0
        ),
        earth_parameters=(1, 1, earth_resistivity),
        transformation=true
    )
    ugc = PowerImpedance.cable(
        length=length_m * regularized,
        positions=[(-0.5, 1), (0.5, 1)],
        C1=PowerImpedance.Conductor(rₒ=0.02622, ρ=2.354e-8, μᵣ=1.035),
        I1=PowerImpedance.Insulator(rᵢ=0.02622, rₒ=0.06006, ϵᵣ=2.67, μᵣ=1.469),
        C2=PowerImpedance.Conductor(rᵢ=0.06006, rₒ=0.06336, ρ=2.14e-7, μᵣ=1.0),
        I2=PowerImpedance.Insulator(rᵢ=0.06336, rₒ=0.06636, ϵᵣ=2.3, μᵣ=1.0),
        C3=PowerImpedance.Conductor(rᵢ=0.06636, rₒ=0.06651, ρ=2.826e-8, μᵣ=1.0),
        I3=PowerImpedance.Insulator(rᵢ=0.06651, rₒ=0.07256, ϵᵣ=2.3, μᵣ=1.0),
        earth_parameters=(1, 1, earth_resistivity),
        transformation=true
    )
    return (; ohl, ugc)
end

function application_network(specification; share=APPLICATION_DEFAULT_SHARE)
    earth = specification["earth_resistivity_ohm_m"]
    corridor = corridor_elements(share, earth, APPLICATION_DEFAULT_LENGTH_M)
    elements = (;
        g1=PowerImpedance.ac_source(
            V=APPLICATION_VOLTAGE_KV,
            P=APPLICATION_POWER_MW,
            P_min=-2000,
            P_max=2000,
            Q_max=1000,
            Q_min=-1000,
            pins=3,
            transformation=true
        ),
        c1=application_converter_one(),
        c2=application_converter_two(),
        corridor.ugc,
        corridor.ohl,
        g4=PowerImpedance.ac_source(
            V=APPLICATION_VOLTAGE_KV,
            P=APPLICATION_POWER_MW,
            P_min=-2000,
            P_max=2000,
            Q_max=1000,
            Q_min=-1000,
            pins=3,
            transformation=true
        ),
        tl1=terminal_line(25.0e3, earth),
        tl78=terminal_line(90.0e3, earth)
    )
    return PowerImpedance.NetworkBuilder.define(
        elements,
        APPLICATION_CONNECTIONS;
        options=APPLICATION_BUILDER_OPTIONS
    )
end

function validate_powerflow_spec(parameters)
    specification = required(parameters, "specification", Dict)
    case_id = string(get(specification, "case_id", APPLICATION_CASE_ID))
    case_id == APPLICATION_CASE_ID || throw(PermanentOperationError(
        "invalid_input",
        "Unsupported power-flow specification: $case_id"
    ))
    normalized = Dict{String,Any}(
        "case_id" => case_id,
        "earth_resistivity_ohm_m" => Float64(get(
            specification,
            "earth_resistivity_ohm_m",
            100.0
        )),
    )
    1 <= normalized["earth_resistivity_ohm_m"] <= 10_000 ||
        throw(PermanentOperationError(
            "invalid_input",
            "Earth resistivity must be between 1 and 10000 ohm m"
        ))
    return Dict{String,Any}("specification" => normalized)
end

powerflow_resource_key(specification) =
    input_hash("powerflow.specification", specification)

function build_prepared_powerflow(context, specification)
    specification_hash = powerflow_resource_key(specification)
    return prepare_resource!(context.prepared_cache, specification_hash) do
        progress!(context, 0.1, "network", message="Constructing reference network")
        network = application_network(specification)
        progress!(context, 0.3, "powerflow", message="Solving reference power flow")
        powerflow = PowerImpedance.compute(
            PowerImpedance.PowerFlowProblem(network),
            PowerImpedance.ACDCPowerFlow()
        )
        progress!(context, 0.65, "linearization", message="Linearizing network model")
        linearization = PowerImpedance.compute(
            PowerImpedance.LinearizationProblem(network, powerflow),
            PowerImpedance.AdmittanceLinearization()
        )
        progress!(context, 1.0, "hot", message="Prepared resource is hot")
        PreparedPowerFlow(
            specification_hash,
            network,
            linearization.network_model,
            powerflow
        )
    end
end

function execute_powerflow_prepare(context, parameters)
    specification = parameters["specification"]
    key = powerflow_resource_key(specification)
    prior_status = prepared_status(context.prepared_cache, key)
    build_prepared_powerflow(context, specification)
    return Dict{String,Any}(
        "prepared_resource_key" => key,
        "specification_hash" => key,
        "cache_status" => prior_status == :hot ? "hit" : "miss",
        "resource_status" => "hot",
    )
end

function validate_impedance_evaluation(parameters)
    prepared = validate_powerflow_spec(parameters)
    prepared["prepared_resource_key"] = string(get(
        parameters,
        "prepared_resource_key",
        ""
    ))
    prepared["ugc_share"] = bounded_real(parameters, "ugc_share", 0.0, 1.0)
    prepared["corridor_length_m"] = bounded_real(
        parameters,
        "corridor_length_m",
        1_000.0,
        1e6
    )
    prepared["length_error_percent"] = bounded_real(
        parameters,
        "length_error_percent",
        0.0,
        50.0
    )
    minimum_frequency = bounded_real(parameters, "minimum_frequency_hz", 1.0, 1e6)
    maximum_frequency = bounded_real(parameters, "maximum_frequency_hz", 1.0, 1e6)
    maximum_frequency > minimum_frequency || throw(PermanentOperationError(
        "invalid_input",
        "maximum_frequency_hz must exceed minimum_frequency_hz"
    ))
    prepared["minimum_frequency_hz"] = minimum_frequency
    prepared["maximum_frequency_hz"] = maximum_frequency
    prepared["frequency_points"] = bounded_integer(parameters, "frequency_points", 2, 2_000)
    return prepared
end

function impedance_cache_inputs(parameters)
    projected = normalize_wire(parameters)
    delete!(projected, "prepared_resource_key")
    return projected
end

function set_corridor_lengths!(network, share, total_length)
    regularized = clamp(share, 1e-3, 1 - 1e-3)
    network.elements.ohl.element_model.length = total_length * (1 - regularized)
    network.elements.ugc.element_model.length = total_length * regularized
    return network
end

function corridor_network_model(prepared, share, total_length)
    if share == APPLICATION_DEFAULT_SHARE &&
            total_length == APPLICATION_DEFAULT_LENGTH_M
        return prepared.network_model
    end

    network = deepcopy(prepared.network)
    set_corridor_lengths!(network, share, total_length)

    # Length changes only the two passive corridor elements. Reuse the
    # prepared active-device linearization and ask PowerImpedance's own element
    # builder for replacement admittance functions for those two elements.
    # This keeps both the power flow and converter steady-state solves out of
    # the evaluation loop without reproducing either line formulation here.
    builder = PowerImpedance.NetworkBuilder
    baseline = prepared.network_model
    admittances = copy(baseline.element_admittances.Y!)
    for name in (:ohl, :ugc)
        index = baseline.indices.elements[name]
        admittances[index] = builder.build(
            getproperty(network.elements, name),
            PowerImpedance.Setpoint()
        )
    end
    lookup = builder.AdmittanceLookup(
        admittances,
        baseline.element_admittances.indices
    )
    return builder.NetworkModel(
        lookup,
        baseline.active_elements,
        baseline.passive_elements,
        baseline.grounded_nodes,
        baseline.retained_nodes,
        baseline.indices
    )
end

function impedance_curve(network_model, parameters)
    frequency_range = (
        parameters["minimum_frequency_hz"],
        parameters["maximum_frequency_hz"],
        parameters["frequency_points"],
    )
    problem = PowerImpedance.PowerImpedanceProblem(
        network_model;
        nodes=[:B5],
        eliminated_elements=[:c2],
        frequency_range
    )
    response = PowerImpedance.compute(problem, PowerImpedance.NodalImpedance())
    nodes = collect(PowerImpedance.response_nodes(response))
    node_index = only(findall(==(:B5), nodes))
    frequencies = real.(PowerImpedance.angular_frequencies(response)) ./ (2π)
    impedance = @view PowerImpedance.response_values(response)[node_index, node_index, :]
    return Dict{String,Any}(
        "frequency_hz" => Float64.(frequencies),
        "magnitude_db_ohm" => Float64.(20 .* log10.(abs.(vec(impedance)))),
    )
end

function execute_impedance_evaluation(context, parameters)
    specification = parameters["specification"]
    authoritative_key = powerflow_resource_key(specification)
    supplied_key = parameters["prepared_resource_key"]
    if !isempty(supplied_key) && supplied_key != authoritative_key
        log!(context, "Ignoring stale prepared-resource hint; reconstructing from specification")
    end
    prepared = build_prepared_powerflow(context, specification)
    base_length = parameters["corridor_length_m"]
    relative_error = parameters["length_error_percent"] / 100
    cases = (
        ("base_minus_error", base_length * (1 - relative_error)),
        ("base", base_length),
        ("base_plus_error", base_length * (1 + relative_error)),
    )
    curves = Dict{String,Any}()
    for (index, (label, length_m)) in enumerate(cases)
        check_canceled(context)
        progress!(context, 0.25 + 0.2 * index, "impedance";
            message="Evaluating $label")
        network_model = corridor_network_model(
            prepared,
            parameters["ugc_share"],
            length_m
        )
        curve = impedance_curve(network_model, parameters)
        curve["corridor_length_m"] = length_m
        curves[label] = curve
    end
    progress!(context, 1.0, "complete")
    return Dict{String,Any}(
        "prepared_resource_key" => authoritative_key,
        "ugc_share" => parameters["ugc_share"],
        "corridor_length_m" => base_length,
        "length_error_percent" => parameters["length_error_percent"],
        "curves" => curves,
    )
end

function register_powerimpedance_operations!(registry)
    register!(registry, OperationSpec(
        "powerflow.prepare",
        validate_powerflow_spec,
        execute_powerflow_prepare;
        timeout_seconds=600,
        cache_policy=:prepared,
        execution_mode=:supervised
    ))
    register!(registry, OperationSpec(
        "impedance.evaluate",
        validate_impedance_evaluation,
        execute_impedance_evaluation;
        cache_key_inputs=impedance_cache_inputs,
        timeout_seconds=600,
        cache_policy=:content,
        execution_mode=:supervised
    ))
    return registry
end
