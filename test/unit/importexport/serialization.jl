@testitem "ImportExport / tagged values / lossless scalar and container round trips" tags=[:unit] setup=[
    ImportExportTestSupport,
    UseImportExportSupport
] begin
    using Logging
    import LineCableModels.ImportExport as IE

    values=Any[
        nothing,
        true,
        7,
        1.25,
        Inf,
        -Inf,
        NaN,
        2.0 - 3.0im,
        "cable",
        :phase,
        [1, 2.5, -1.0im],
        (1.0, 2.0),
        Dict(:count=>2, :scale=>0.5)
    ]
    for original in values
        encoded=IE._serialize_value(original)
        decoded=IE._deserialize_value(encoded)
        if original isa AbstractFloat&&isnan(original)
            @test isnan(decoded)
        elseif original isa Symbol
            @test decoded == String(original)
        elseif original isa Tuple
            @test decoded == collect(original)
        elseif original isa Dict
            @test decoded == Dict(:count => 2, :scale => 0.5)
        else
            @test decoded == original
        end
    end

    @test_throws ArgumentError IE._deserialize_value(Dict("__type__" => "Complex"))
    unknown_float=@test_logs (
        :warn,
        r"Unknown SpecialFloat"
    ) IE._deserialize_value(Dict("__type__" => "SpecialFloat", "value" => "huge"))
    @test unknown_float["__type__"] == "SpecialFloat"

    unknown_marker=@test_logs (
        :warn,
        r"Unknown __type__ marker"
    ) IE._deserialize_value(Dict("__type__" => "FutureScalar", "value" => 2))
    @test unknown_marker == Dict(Symbol("__type__") => "FutureScalar", :value => 2)

    @test IE.get_as(Dict{Symbol, Any}(), :absent, 4.0, Float64) == 4.0
    @test ismissing(IE.get_as(Dict(:value => missing), :value, 0.0, Float64))
    @test IE.get_as(Dict("value" => Dict("__type__" => "Int", "value" => 4)), "value", 0, Int) ==
          4
end

@testitem "ImportExport / cable designs / malformed structure equivalence classes" tags=[:unit] setup=[
    ImportExportTestSupport,
    UseImportExportSupport,
    TestFixtures
] begin
    using Logging
    import LineCableModels.ImportExport as IE

    encoded=IE._serialize_value(TestFixtures.mv_cable_design())

    missing_nominal=deepcopy(encoded)
    delete!(missing_nominal, "nominal_data")
    restored=@test_logs (:warn, r"Missing 'nominal_data'") match_mode = :any IE._reconstruct_cabledesign(
        "missing-nominal",
        missing_nominal
    )
    @test isnothing(restored.nominal_data.designation_code)
    @test length(restored.components) == length(encoded["components"])

    malformed_nominal=deepcopy(encoded)
    malformed_nominal["nominal_data"]=Dict("unexpected"=>1)
    @test_throws ErrorException IE._reconstruct_cabledesign(
        "malformed-nominal",
        malformed_nominal
    )

    skipped_component=deepcopy(encoded)
    skipped_component["components"]=Any[3, skipped_component["components"]...]
    restored=@test_logs (:warn, r"Component data at index 1") match_mode = :any IE._reconstruct_cabledesign(
        "skipped-component",
        skipped_component
    )
    @test length(restored.components) == length(encoded["components"])

    no_components=deepcopy(encoded)
    no_components["components"]=Any[]
    @test_throws ErrorException IE._reconstruct_cabledesign("empty", no_components)

    no_valid_components=deepcopy(encoded)
    no_valid_components["components"]=Any[3]
    @test_logs (:warn, r"not a dictionary") match_mode = :any begin
        @test_throws ErrorException IE._reconstruct_cabledesign(
            "no-valid-components",
            no_valid_components
        )
    end

    invalid_first_conductor=deepcopy(encoded)
    invalid_first_conductor["components"][1]["conductor_group"]["layers"]=Any[
        3, invalid_first_conductor["components"][1]["conductor_group"]["layers"][2:end]...]
    @test_throws ErrorException IE._reconstruct_cabledesign(
        "invalid-first-conductor",
        invalid_first_conductor
    )

    skipped_conductor=deepcopy(encoded)
    conductor_layers=Any[
        skipped_conductor["components"][1]["conductor_group"]["layers"]...,
    ]
    conductor_layers[2]=3
    skipped_conductor["components"][1]["conductor_group"]["layers"]=conductor_layers
    @test_logs (:warn, r"Conductor layer data at index 2") match_mode = :any begin
        @test_throws ArgumentError IE._reconstruct_cabledesign(
            "skipped-conductor",
            skipped_conductor
        )
    end

    malformed_conductor=deepcopy(encoded)
    delete!(
        malformed_conductor["components"][1]["conductor_group"]["layers"][2],
        "radius_wire"
    )
    redirect_stderr(devnull) do
        @test_logs (:error, r"Failed to add conductor layer 2") match_mode = :any begin
            @test_throws ErrorException IE._reconstruct_cabledesign(
                "malformed-conductor",
                malformed_conductor
            )
        end
    end

    no_insulators=deepcopy(encoded)
    no_insulators["components"][1]["insulator_group"]["layers"]=Any[]
    @test_throws ErrorException IE._reconstruct_cabledesign("no-insulators", no_insulators)

    invalid_first_insulator=deepcopy(encoded)
    invalid_first_insulator["components"][1]["insulator_group"]["layers"]=Any[
        3, invalid_first_insulator["components"][1]["insulator_group"]["layers"][2:end]...]
    @test_throws ErrorException IE._reconstruct_cabledesign(
        "invalid-first-insulator",
        invalid_first_insulator
    )

    skipped_insulator=deepcopy(encoded)
    insulator_layers=Any[
        skipped_insulator["components"][1]["insulator_group"]["layers"]...,
    ]
    insulator_layers[2]=3
    skipped_insulator["components"][1]["insulator_group"]["layers"]=insulator_layers
    restored=@test_logs (:warn, r"Insulator layer data at index 2") match_mode = :any IE._reconstruct_cabledesign(
        "skipped-insulator",
        skipped_insulator
    )
    @test length(restored.components[1].insulator_group.layers) ==
          length(encoded["components"][1]["insulator_group"]["layers"]) - 1

    malformed_insulator=deepcopy(encoded)
    delete!(
        malformed_insulator["components"][1]["insulator_group"]["layers"][2],
        "r_ex"
    )
    redirect_stderr(devnull) do
        @test_logs (:error, r"Failed to add insulator layer 2") match_mode = :any begin
            @test_throws ErrorException IE._reconstruct_cabledesign(
                "malformed-insulator",
                malformed_insulator
            )
        end
    end
end

@testitem "ImportExport / type reconstruction / resolution and constructor fallbacks" tags=[:unit] setup=[
    ImportExportTestSupport,
    UseImportExportSupport
] begin
    using Logging
    import LineCableModels.ImportExport as IE

    @test IE._resolve_type("Material") === Material
    @test IE._resolve_type("LineCableModels.DataModel.CableConstants") === CableConstants
    @test IE._resolve_type("Vector{Float64}") === Vector{Float64}

    redirect_stderr(devnull) do
        @test_logs (:error, r"Could not resolve type") begin
            @test_throws Exception IE._resolve_type("Definitely.Not.A.Type")
        end
    end

    encoded_constants=Dict(
        "__julia_type__"=>"LineCableModels.DataModel.CableConstants",
        "R"=>IE._serialize_value(1.0),
        "L"=>IE._serialize_value(2.0),
        "C"=>IE._serialize_value(3.0),
        "ignored"=>99
    )
    @test IE._deserialize_value(encoded_constants) == CableConstants(1.0, 2.0, 3.0)

    encoded_material=IE._serialize_value(Material(1.0, 2.0, 3.0, 20.0, 0.01))
    @test IE._deserialize_value(encoded_material) == Material(1.0, 2.0, 3.0, 20.0, 0.01)

    unresolved=redirect_stderr(devnull) do
        @test_logs match_mode = :any (:error, r"Could not resolve type") (
            :error,
            r"Failed to resolve or deserialize type"
        ) IE._deserialize_value(Dict(
            "__julia_type__" => "Definitely.Not.A.Type",
            "field" => 1
        ))
    end
    @test unresolved["field"] == 1

    incomplete=Dict(
        "__julia_type__"=>"LineCableModels.DataModel.CableConstants",
        "R"=>1.0
    )
    returned=redirect_stderr(devnull) do
        @test_logs match_mode = :any (
            :error,
            r"Positional construction failed"
        ) (:error, r"Failed to resolve or deserialize type") IE._deserialize_value(incomplete)
    end
    @test returned === incomplete
end

@testitem "ImportExport / cable layers / semantic reconstruction classes" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport
] begin
    using Logging
    import LineCableModels.ImportExport as IE

    material=IE._serialize_value(Material(1.7241e-8, 2.3, 1.0, 20.0, 0.00393))
    tagged(type_name; fields...) = Dict(
        "__julia_type__"=>"LineCableModels.DataModel.$type_name",
        "r_in"=>IE._serialize_value(0.01),
        "material_props"=>material,
        (string(key)=>IE._serialize_value(value) for (key, value) in fields)...
    )

    cases=(
        ("CircStrands", CircStrands,
            (radius_wire = 0.001, num_wires = 12, lay_ratio = 10.0, lay_direction = -1)),
        ("Tubular", Tubular, (r_ex = 0.02,)),
        ("Strip", Strip, (r_ex = 0.02, width = 0.002, lay_ratio = 12.0, lay_direction = 1)),
        ("Insulator", Insulator, (r_ex = 0.02,)),
        ("Semicon", Semicon, (r_ex = 0.02,))
    )
    for (type_name, expected_type, fields) in cases
        reconstructed=IE._reconstruct_partsgroup(tagged(type_name; fields...))
        @test reconstructed isa expected_type
        @test reconstructed.r_in == 0.01
        @test reconstructed.material_props isa Material
    end

    @test_throws ErrorException IE._reconstruct_partsgroup(Dict("r_in" => 0.0))
    @test_throws ErrorException IE._reconstruct_partsgroup(Dict(
        "__julia_type__" => "LineCableModels.DataModel.Tubular",
        "material_props" => material,
        "r_ex" => 0.02
    ))
    @test_throws ErrorException IE._reconstruct_partsgroup(Dict(
        "__julia_type__" => "LineCableModels.DataModel.Tubular",
        "r_in" => 0.01,
        "r_ex" => 0.02
    ))
    @test_logs (:error, r"Construction failed") begin
        @test_throws ErrorException IE._reconstruct_partsgroup(tagged("CableComponent"))
    end

    malformed_tubular=tagged("Tubular")
    redirect_stderr(devnull) do
        @test_logs (:error, r"Construction failed") begin
            @test_throws ErrorException IE._reconstruct_partsgroup(malformed_tubular)
        end
    end
end
