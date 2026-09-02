const DEFAULT_CABLE_PARAMETERS = Dict{String,Any}(
    "core_radius_m" => 0.01,
    "insulation_thickness_m" => 0.005,
    "sheath_thickness_m" => 0.001,
    "conductor_resistivity_ohm_m" => 1.7241e-8,
    "conductor_relative_permeability" => 1.0,
    "insulation_resistivity_ohm_m" => 1.0e14,
    "insulation_relative_permittivity" => 2.3,
    "temperature_celsius" => 20.0,
)

function cable_inputs(parameters)
    merged = merge(DEFAULT_CABLE_PARAMETERS, parameters)
    core_radius = bounded_real(merged, "core_radius_m", 1e-5, 0.5)
    insulation = bounded_real(merged, "insulation_thickness_m", 1e-5, 0.5)
    sheath = bounded_real(merged, "sheath_thickness_m", 1e-6, 0.1)
    conductor_rho = bounded_real(
        merged,
        "conductor_resistivity_ohm_m",
        1e-10,
        1e-3
    )
    conductor_mu = bounded_real(
        merged,
        "conductor_relative_permeability",
        0.1,
        10_000.0
    )
    insulation_rho = bounded_real(
        merged,
        "insulation_resistivity_ohm_m",
        1.0,
        1e20
    )
    insulation_eps = bounded_real(
        merged,
        "insulation_relative_permittivity",
        1.0,
        100.0
    )
    temperature = bounded_real(merged, "temperature_celsius", -50.0, 150.0)
    return Dict{String,Any}(
        "core_radius_m" => core_radius,
        "insulation_thickness_m" => insulation,
        "sheath_thickness_m" => sheath,
        "conductor_resistivity_ohm_m" => conductor_rho,
        "conductor_relative_permeability" => conductor_mu,
        "insulation_resistivity_ohm_m" => insulation_rho,
        "insulation_relative_permittivity" => insulation_eps,
        "temperature_celsius" => temperature,
    )
end

function build_coaxial_design(parameters)
    core_radius = parameters["core_radius_m"]
    insulation_outer = core_radius + parameters["insulation_thickness_m"]
    sheath_outer = insulation_outer + parameters["sheath_thickness_m"]
    conductor = LineCableModels.Material(
        :conductor,
        parameters["conductor_resistivity_ohm_m"],
        parameters["conductor_relative_permeability"],
        1.0,
        20.0,
        0.00393
    )
    dielectric = LineCableModels.Material(
        :insulator,
        parameters["insulation_resistivity_ohm_m"],
        parameters["insulation_relative_permittivity"]
    )
    radial = LineCableModels.Stack(
        LineCableModels.Group(
            :core,
            LineCableModels.Region(
                :core_metal,
                LineCableModels.Disk(core_radius),
                conductor
            )
        ),
        LineCableModels.Region(
            :insulation,
            LineCableModels.Annulus(core_radius, insulation_outer),
            dielectric
        ),
        LineCableModels.Group(
            :sheath,
            LineCableModels.Region(
                :sheath_metal,
                LineCableModels.Annulus(insulation_outer, sheath_outer),
                conductor
            )
        )
    )
    return LineCableModels.build(
        LineCableModels.CableDesign,
        "playground-coaxial",
        radial
    )
end

function geometry_summary(_, parameters)
    design = build_coaxial_design(parameters)
    regions = Any[
        Dict{String,Any}(
            "name" => string(region.source.tag),
            "primitive" => string(nameof(typeof(region.primitive))),
            "area_m2" => Float64(LineCableModels.area(region.primitive)),
            "terminal" => isnothing(region.terminal) ? nothing : string(region.terminal),
        )
        for region in design.geometry.regions
    ]
    return Dict{String,Any}(
        "cable_id" => design.cable_id,
        "outer_radius_m" => Float64(LineCableModels.outer_radius(design)),
        "terminals" => string.(design.terminal_order),
        "regions" => regions,
    )
end

function constants_row(result, index)
    return Dict{String,Any}(
        "terminal" => string(result.cores[index]),
        "resistance_ohm_per_m" => Float64(result.R[index]),
        "inductance_h_per_m" => Float64(result.L[index]),
        "conductance_s_per_m" => Float64(result.G[index]),
        "capacitance_f_per_m" => Float64(result.C[index]),
    )
end

function validate_cable_constants(parameters)
    values = cable_inputs(parameters)
    values["frequency_hz"] = bounded_real(parameters, "frequency_hz", 1e-3, 1e8)
    return values
end

function compute_cable_constants(context, parameters)
    progress!(context, 0.1, "constructing", message="Constructing CableDesign")
    design = build_coaxial_design(parameters)
    progress!(context, 0.35, "solving", message="Evaluating package formulas")
    problem = LineCableModels.CableConstantsProblem(
        design;
        temperature=parameters["temperature_celsius"],
        frequency=parameters["frequency_hz"]
    )
    result = LineCableModels.compute(
        problem,
        LineCableModels.CableConstantsFormulation()
    )
    progress!(context, 1.0, "complete")
    return Dict{String,Any}(
        "frequency_hz" => parameters["frequency_hz"],
        "constants" => Any[constants_row(result, index)
            for index in eachindex(result.cores)],
    )
end

function validate_frequency_sweep(parameters)
    values = cable_inputs(parameters)
    frequencies = required(parameters, "frequencies_hz", Vector)
    1 <= length(frequencies) <= 5_000 || throw(PermanentOperationError(
        "invalid_input",
        "frequencies_hz must contain between 1 and 5000 values"
    ))
    all(value -> value isa Real && 1e-3 <= value <= 1e8, frequencies) ||
        throw(PermanentOperationError(
            "invalid_input",
            "frequencies_hz must contain finite values between 1e-3 and 1e8 Hz"
        ))
    issorted(frequencies) || throw(PermanentOperationError(
        "invalid_input",
        "frequencies_hz must be sorted"
    ))
    values["frequencies_hz"] = Float64.(frequencies)
    return values
end

function compute_frequency_sweep(context, parameters)
    design = build_coaxial_design(parameters)
    frequencies = parameters["frequencies_hz"]
    rows = Any[]
    for (index, frequency) in enumerate(frequencies)
        check_canceled(context)
        result = LineCableModels.compute(
            LineCableModels.CableConstantsProblem(
                design;
                temperature=parameters["temperature_celsius"],
                frequency
            ),
            LineCableModels.CableConstantsFormulation()
        )
        push!(rows, Dict{String,Any}(
            "frequency_hz" => frequency,
            "constants" => Any[constants_row(result, terminal)
                for terminal in eachindex(result.cores)],
        ))
        progress!(context, index / length(frequencies), "frequency_sweep";
            message="Frequency $index of $(length(frequencies))")
    end
    return Dict{String,Any}("samples" => rows)
end

complex_wire(value) = Dict{String,Any}(
    "real" => Float64(real(value)),
    "imag" => Float64(imag(value)),
)

function tensor_wire(values)
    return Any[
        Any[
            Any[complex_wire(values[row, column, frequency])
                for frequency in axes(values, 3)]
            for column in axes(values, 2)]
        for row in axes(values, 1)
    ]
end

function validate_line_parameters(parameters)
    values = validate_frequency_sweep(parameters)
    values["separation_m"] = bounded_real(parameters, "separation_m", 0.05, 100.0)
    values["depth_m"] = bounded_real(parameters, "depth_m", 0.05, 100.0)
    values["earth_resistivity_ohm_m"] = bounded_real(
        parameters,
        "earth_resistivity_ohm_m",
        0.01,
        1e6
    )
    values["line_length_m"] = bounded_real(parameters, "line_length_m", 0.01, 1e7)
    return values
end

function compute_line_parameters(context, parameters)
    design = build_coaxial_design(parameters)
    separation = parameters["separation_m"]
    depth = parameters["depth_m"]
    connections = [Dict(:core => 1, :sheath => 0), Dict(:core => 2, :sheath => 0)]
    progress!(context, 0.1, "constructing", message="Constructing LineCableSystem")
    system = LineCableModels.build(
        LineCableModels.LineCableSystem,
        [design, design],
        [
            LineCableModels.Pose2(-separation / 2, -depth, 0.0),
            LineCableModels.Pose2(separation / 2, -depth, 0.0),
        ];
        connections,
        system_id="playground-two-cable-system",
        line_length=parameters["line_length_m"]
    )
    problem = LineCableModels.LineParametersProblem(
        system;
        temperature=parameters["temperature_celsius"],
        earth_props=LineCableModels.Earth(
            rho=parameters["earth_resistivity_ohm_m"]
        ),
        frequencies=parameters["frequencies_hz"]
    )
    progress!(context, 0.25, "solving", message="Evaluating line parameters")
    result = LineCableModels.compute(problem, LineCableModels.Formulation())
    progress!(context, 1.0, "complete")
    return Dict{String,Any}(
        "frequencies_hz" => parameters["frequencies_hz"],
        "series_impedance_ohm_per_m" => tensor_wire(LineCableModels.Z(result)),
        "shunt_admittance_s_per_m" => tensor_wire(LineCableModels.Y(result)),
    )
end

function register_linecablemodels_operations!(registry)
    register!(registry, OperationSpec(
        "cable.geometry_summary",
        cable_inputs,
        geometry_summary;
        timeout_seconds=60,
        execution_mode=:supervised
    ))
    register!(registry, OperationSpec(
        "cable.constants",
        validate_cable_constants,
        compute_cable_constants;
        timeout_seconds=60,
        execution_mode=:supervised
    ))
    register!(registry, OperationSpec(
        "cable.frequency_sweep",
        validate_frequency_sweep,
        compute_frequency_sweep;
        timeout_seconds=180,
        execution_mode=:supervised
    ))
    register!(registry, OperationSpec(
        "line.frequency_scan",
        validate_line_parameters,
        compute_line_parameters;
        timeout_seconds=300,
        execution_mode=:supervised
    ))
    return registry
end
