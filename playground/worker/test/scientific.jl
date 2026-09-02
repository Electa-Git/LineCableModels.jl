function wire_complex(value)
    return ComplexF64(value["real"], value["imag"])
end

@testset "LineCableModels scientific parity" begin
    Worker.load_scientific_packages!()
    engine = Worker.LineCableModels
    registry = Worker.default_registry()
    context, _ = test_context()

    constants_parameters = Worker.validate_cable_constants(Dict{String,Any}(
        "frequency_hz" => 50.0,
    ))
    constants = Worker.execute_operation(
        Worker.registered_operation(registry, "cable.constants"),
        context,
        constants_parameters
    )
    design = Worker.build_coaxial_design(constants_parameters)
    direct_constants = engine.compute(
        engine.CableConstantsProblem(
            design;
            temperature=constants_parameters["temperature_celsius"],
            frequency=constants_parameters["frequency_hz"]
        ),
        engine.CableConstantsFormulation()
    )
    @test length(constants["constants"]) == length(direct_constants.cores)
    for index in eachindex(direct_constants.cores)
        row = constants["constants"][index]
        @test row["resistance_ohm_per_m"] ≈ direct_constants.R[index]
        @test row["inductance_h_per_m"] ≈ direct_constants.L[index]
        @test row["conductance_s_per_m"] ≈ direct_constants.G[index]
        @test row["capacitance_f_per_m"] ≈ direct_constants.C[index]
    end

    sweep = Worker.execute_operation(
        Worker.registered_operation(registry, "cable.frequency_sweep"),
        context,
        Dict{String,Any}("frequencies_hz" => [1.0, 50.0, 2_500.0])
    )
    @test getindex.(sweep["samples"], "frequency_hz") == [1.0, 50.0, 2_500.0]
    @test all(sample -> all(
        value -> value isa Real && isfinite(value),
        Iterators.flatten((
            (
                row["resistance_ohm_per_m"],
                row["inductance_h_per_m"],
                row["conductance_s_per_m"],
                row["capacitance_f_per_m"],
            ) for row in sample["constants"]
        ))
    ), sweep["samples"])

    line_parameters = Worker.validate_line_parameters(Dict{String,Any}(
        "frequencies_hz" => [1.0, 50.0, 2_500.0],
        "separation_m" => 1.0,
        "depth_m" => 1.0,
        "earth_resistivity_ohm_m" => 100.0,
        "line_length_m" => 1_000.0,
    ))
    line_result = Worker.execute_operation(
        Worker.registered_operation(registry, "line.frequency_scan"),
        context,
        line_parameters
    )
    line_design = Worker.build_coaxial_design(line_parameters)
    line_system = engine.build(
        engine.LineCableSystem,
        [line_design, line_design],
        [engine.Pose2(-0.5, -1.0, 0.0), engine.Pose2(0.5, -1.0, 0.0)];
        connections=[Dict(:core => 1, :sheath => 0), Dict(:core => 2, :sheath => 0)],
        system_id="scientific-parity-system",
        line_length=1_000.0
    )
    direct_line = engine.compute(
        engine.LineParametersProblem(
            line_system;
            temperature=line_parameters["temperature_celsius"],
            earth_props=engine.Earth(rho=100.0),
            frequencies=line_parameters["frequencies_hz"]
        ),
        engine.Formulation()
    )
    @test line_result["frequencies_hz"] == [1.0, 50.0, 2_500.0]
    @test wire_complex(line_result["series_impedance_ohm_per_m"][1][1][1]) ≈
          engine.Z(direct_line)[1, 1, 1]
    @test wire_complex(line_result["series_impedance_ohm_per_m"][1][2][3]) ≈
          engine.Z(direct_line)[1, 2, 3]
    @test wire_complex(line_result["shunt_admittance_s_per_m"][2][2][2]) ≈
          engine.Y(direct_line)[2, 2, 2]
end

if get(ENV, "LCM_TEST_POWERIMPEDANCE", "0") == "1"
    @testset "PowerImpedance prepared execution" begin
        registry = Worker.default_registry()
        context, _ = test_context()
        supervisor = Worker.ExecutorSupervisor()
        specification = Dict{String,Any}(
            "case_id" => "ohl_ugc_transition_v1",
            "earth_resistivity_ohm_m" => 100.0,
        )
        try
            prepare = Worker.registered_operation(registry, "powerflow.prepare")
            first = Worker.execute_supervised!(
                supervisor,
                prepare,
                context,
                Dict{String,Any}("specification" => specification)
            )
            second = Worker.execute_supervised!(
                supervisor,
                prepare,
                context,
                Dict{String,Any}("specification" => specification)
            )
            @test first["cache_status"] == "miss"
            @test second["cache_status"] == "hit"
            @test first["prepared_resource_key"] == second["prepared_resource_key"]

            impedance = Worker.execute_supervised!(
                supervisor,
                Worker.registered_operation(registry, "impedance.evaluate"),
                context,
                Dict{String,Any}(
                    "specification" => specification,
                    "prepared_resource_key" => first["prepared_resource_key"],
                    "ugc_share" => 0.5,
                    "corridor_length_m" => 100_000.0,
                    "length_error_percent" => 5.0,
                    "minimum_frequency_hz" => 100.0,
                    "maximum_frequency_hz" => 150.0,
                    "frequency_points" => 2,
                )
            )
            @test sort!(collect(keys(impedance["curves"]))) == [
                "base",
                "base_minus_error",
                "base_plus_error",
            ]
            @test impedance["curves"]["base_minus_error"]["corridor_length_m"] ==
                  95_000.0
            @test impedance["curves"]["base"]["corridor_length_m"] == 100_000.0
            @test impedance["curves"]["base_plus_error"]["corridor_length_m"] ==
                  105_000.0
            @test all(
                curve -> length(curve["frequency_hz"]) == 2 &&
                    all(isfinite, curve["magnitude_db_ohm"]),
                values(impedance["curves"])
            )
            @test impedance["curves"]["base_minus_error"]["magnitude_db_ohm"] !=
                  impedance["curves"]["base_plus_error"]["magnitude_db_ohm"]
        finally
            Worker.stop_executor!(supervisor)
        end
    end
end
