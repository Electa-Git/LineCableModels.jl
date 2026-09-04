@testitem "Cable model v1 / DataModel and normalized-input preservation" tags=[:integration] setup=[
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
    @test String.(design.terminal_order) ==
          String.(expected_design["terminal_order"])

    @test design.root isa Stack
    @test all(region -> region isa PlacedRegion, design.geometry.regions)
    @test design.terminal_map == [region.terminal === nothing ? 0 :
           only(findall(==(region.terminal), design.terminal_order))
           for region in design.geometry.regions]

    constants=CableConstants(design)
    constants_problem=CableConstantsProblem(design)
    @test constants_problem.design === design
    @test constants == compute(constants_problem, CableConstantsFormulation())
    @test all(isfinite, Iterators.flatten(
        (constants.R, constants.L, constants.C, constants.G)
    ))

    system=TestFixtures.three_phase_system()
    expected_system=reference["system"]
    @test system.system_id == expected_system["system_id"]
    @test system.line_length == expected_system["line_length"]
    @test length(system.designs) == length(expected_system["positions"])
    for (design_data,
        position,
        connections,
        expected) in zip(
        system.designs,
        system.positions,
        system.connections,
        expected_system["positions"]
    )
        @test design_data.cable_id == expected["cable_id"]
        @test position.x == expected["horizontal"]
        @test position.y == expected["vertical"]
        @test connections == Int.(expected["connections"])
    end
    actual_primitive_order=[(cable = entry.cable, terminal = String(entry.terminal))
                            for entry in system.terminal_order]
    expected_primitive_order=[(cable = Int(entry["cable"]),
                                  terminal = String(entry["terminal"]))
                              for entry in expected_system["primitive_order"]]
    @test actual_primitive_order == expected_primitive_order

    problem=TestFixtures.line_parameters_problem(frequencies = [50.0, 500.0])
    execution=computation_options(LineCableModelsCoaxial, (;))
    blueprints=LineCableModels.Engine.CableBlueprint{eltype(problem)}[LineCableModels.Engine.flatten(
                                                                          LineCableModelsCoaxial(),
                                                                          design,
                                                                          eltype(problem))
                                                                      for design in problem.system.designs]
    workspace=LineParametersWorkspace(
        problem, Formulation(), execution, blueprints)
    input=workspace.input
    expected_input=reference["normalized_input"]
    vector_fields=(:freq, :horz, :vert, :phase_map, :cable_map)
    for field in vector_fields
        @test getproperty(input, field) == collect(expected_input[string(field)])
    end
    cable_fields=(
        :r_in=>:r_in,
        :r_ext=>:r_ext,
        :r_ins_in=>:r_ins_in,
        :r_ins_ext=>:r_ins_ext,
        :rho0_cond=>:rho0_cond,
        :T0_cond=>:T0_cond,
        :alpha_cond=>:alpha_cond,
        :mu_cond=>:mu_cond,
        :mu_ins=>:mu_ins,
        :r_layer_in=>:r_ins_layer_in,
        :r_layer_ext=>:r_ins_layer_ext
    )
    for (field, reference) in cable_fields
        @test isapprox(
            getproperty(input.cable, field),
            collect(expected_input[string(reference)]);
            rtol = 1.0e-12,
            atol = 0.0
        )
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
    @test input.cable.dielectric_ranges ==
          [Int(range["first"]):Int(range["last"])
           for range in expected_input["insulator_layer_ranges"]]
    for field in (:line_length, :n_frequencies, :n_phases, :n_cables)
        @test getproperty(input, field) == expected_input[string(field)]
    end
    @test input.earth.vertical_layers ==
          expected_input["earth"]["vertical_layers"]
    @test length(input.earth.layers) ==
          length(expected_input["earth"]["layers"])

    parameters=compute(
        problem,
        Formulation(options = (ideal_transposition = false,))
    )
    expected_parameters=reference["line_parameters"]
    @test parameters.f == collect(expected_parameters["frequencies"])
    @test string(basis(parameters)) == expected_parameters["basis"]
    @test size(parameters.Z) == (3, 3, 2)
    @test size(parameters.Y) == (3, 3, 2)
    for frequency in eachindex(parameters.f)
        impedance=@view parameters.Z.values[:, :, frequency]
        admittance=@view parameters.Y.values[:, :, frequency]
        @test all(isfinite, impedance)
        @test all(isfinite, admittance)
        @test impedance ≈ transpose(impedance)
        @test admittance ≈ transpose(admittance)
    end
end
