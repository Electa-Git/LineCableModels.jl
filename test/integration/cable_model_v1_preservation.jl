@testitem "Cable model v1 / reference payload and numerical preservation" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures
] begin
    using JSON3

    reference=JSON3.read(read(
        joinpath(
            pkgdir(LineCableModels), "test", "fixtures", "reference",
            "cable_model_v1_preservation.json"
        ),
        String
    ))
    @test reference["reference_sha"] ==
          "36dfe57d330eef8dcdb00f35b905954be84401e5"

    design=TestFixtures.mv_cable_design()
    expected_design=reference["design"]
    @test design.cable_id == expected_design["cable_id"]
    @test getproperty.(design.components, :id) ==
          String.(expected_design["terminal_order"])

    material_fields=(:rho, :eps_r, :mu_r, :T0, :alpha)
    conductor_fields=(
        :r_in, :r_ex, :cross_section, :num_wires, :num_turns,
        :resistance, :alpha, :gmr
    )
    insulation_fields=(
        :r_in, :r_ex, :cross_section, :shunt_capacitance,
        :shunt_conductance, :reference_frequency
    )
    insulation_keys=(
        "r_in", "r_ex", "cross_section", "capacitance",
        "conductance", "reference_frequency"
    )
    for (component, expected) in zip(design.components, expected_design["components"])
        @test component.id == expected["id"]
        conductor=expected["conductor"]
        for field in conductor_fields
            @test getproperty(component.conductor_group, field) == conductor[string(field)]
        end
        for field in material_fields
            @test getproperty(component.conductor_props, field) ==
                  conductor["effective_material"][string(field)]
        end

        insulation=expected["insulation"]
        for (field, key) in zip(insulation_fields, insulation_keys)
            @test getproperty(component.insulator_group, field) == insulation[key]
        end
        for field in material_fields
            @test getproperty(component.insulator_props, field) ==
                  insulation["effective_material"][string(field)]
        end
        @test length(component.insulator_group.layers) ==
              length(insulation["raw_layers"])
        for (layer, expected_layer) in
            zip(component.insulator_group.layers, insulation["raw_layers"])
            @test layer.r_in == expected_layer["r_in"]
            @test layer.r_ex == expected_layer["r_ex"]
            for field in material_fields
                @test getproperty(layer.material_props, field) ==
                      expected_layer["material"][string(field)]
            end
        end
    end

    constants=compute(CableConstantsProblem(design), Formulation())
    expected_constants=expected_design["cable_constants"]
    @test constants.R == expected_constants["R"]
    @test constants.L == expected_constants["L"]
    @test constants.C == expected_constants["C"]

    system=TestFixtures.three_phase_system()
    expected_system=reference["system"]
    @test system.system_id == expected_system["system_id"]
    @test system.line_length == expected_system["line_length"]
    @test length(system.cables) == length(expected_system["positions"])
    for (position, expected) in zip(system.cables, expected_system["positions"])
        @test position.design_data.cable_id == expected["cable_id"]
        @test position.horz == expected["horizontal"]
        @test position.vert == expected["vertical"]
        @test position.conn == Int.(expected["connections"])
    end
    actual_primitive_order=[
        (cable=index, terminal=component.id)
        for (index, position) in enumerate(system.cables)
        for component in position.design_data.components
    ]
    expected_primitive_order=[
        (cable=Int(entry["cable"]), terminal=String(entry["terminal"]))
        for entry in expected_system["primitive_order"]
    ]
    @test actual_primitive_order == expected_primitive_order

    problem=TestFixtures.line_parameters_problem(frequencies = [50.0, 500.0])
    input=LineCableModels.Engine.AnalyticalInput(problem)
    expected_input=reference["analytical_input"]
    vector_fields=(
        :freq, :horz, :vert, :r_in, :r_ext, :r_ins_in, :r_ins_ext,
        :rho0_cond, :T0_cond, :alpha_cond, :mu_cond, :eps_cond,
        :rho_ins, :mu_ins, :eps_ins, :tan_ins, :r_ins_layer_in,
        :r_ins_layer_ext, :rho_ins_layer, :eps_ins_layer, :phase_map, :cable_map
    )
    for field in vector_fields
        @test getproperty(input, field) == collect(expected_input[string(field)])
    end
    @test input.phase_map == Int.(expected_system["connection_order"])
    @test input.jω == complex.(
        collect(expected_input["jω"]["real"]),
        collect(expected_input["jω"]["imag"])
    )
    expected_separation=expected_input["horz_sep"]
    @test input.horz_sep == reshape(
        collect(expected_separation["values"]),
        Tuple(Int.(expected_separation["size"]))
    )
    @test input.insulator_layer_ranges == [
        Int(range["first"]):Int(range["last"])
        for range in expected_input["insulator_layer_ranges"]
    ]
    for field in (:line_length, :n_frequencies, :n_phases, :n_cables)
        @test getproperty(input, field) == expected_input[string(field)]
    end
    @test input.earth.vertical_layers == expected_input["earth"]["vertical_layers"]
    @test length(input.earth.layers) == length(expected_input["earth"]["layers"])

    parameters=compute(
        problem,
        Formulation(Val(:analytical); options = (ideal_transposition = false,))
    )
    expected_parameters=reference["line_parameters"]
    function expected_tensor(node)
        return reshape(
            complex.(collect(node["real"]), collect(node["imag"])),
            Tuple(Int.(node["size"]))
        )
    end
    @test parameters.Z.values == expected_tensor(expected_parameters["Z"])
    @test parameters.Y.values == expected_tensor(expected_parameters["Y"])
    @test parameters.f == collect(expected_parameters["frequencies"])
    @test string(basis(parameters)) == expected_parameters["basis"]
end
