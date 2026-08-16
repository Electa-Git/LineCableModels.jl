@testitem "ImportExport / VDE parser / ordered tokens and residue" tags=[:unit] setup=[
    DataModelTestSupport,
    UseDataModelSupport
] begin
    parser=LineCableModels.DataModel.vdeparse

    cases=(
        (
            "2XS(F)2Y 3x240/25 12/20kV RM/V",
            Dict(
                :conductor=>"copper conductor",
                :insulation=>"cross-linked PE (XLPE)",
                :waterblocking=>"longitudinally water-proof protection",
                :outer_sheath=>"PE outer sheath",
                :cores=>"3",
                :conductor_cross_section=>"240",
                :metallic_screen_cross_section=>"25",
                :voltage=>"12/20 kV",
                :conductor_type=>"round, stranded, compact"
            )
        ),
        (
            "A2XS(FL)KL2Y 1x630 18/30kV RE",
            Dict(
                :conductor=>"aluminium conductor",
                :sheath=>"aluminium sheath",
                :cores=>"1",
                :conductor_type=>"round, solid"
            )
        ),
        ("N2XSY 1x150", Dict(:designation=>"DIN VDE standard", :cores=>"1"))
    )

    for (code, expected) in cases
        parsed=parser(code)
        @test all(key -> parsed[key] == expected[key], keys(expected))
    end

    @test isempty(parser(" \u00a0 "))
    @test LineCableModels.DataModel.decode_type("RRM") == "round, stranded"
    @test_logs (:warn, r"Unknown conductor type") begin
        @test LineCableModels.DataModel.decode_type("RZ") == "round"
    end
    @test_logs (:warn, r"Unparsed stub residue") begin
        parsed = parser("2XSQ 1x10")
        @test parsed[:unparsed_stub] == "Q"
    end
    @test_logs (:warn, r"Unparsed trailing") begin
        parsed = parser("2XSY 1x10 unsupported")
        @test parsed[:unparsed] == "unsupported"
    end
end

@testitem "ImportExport / MaterialsLibrary / JSON round trip and malformed entries" tags=[:integration] setup=[
    ImportExportTestSupport,
    UseImportExportSupport
] begin
    using JSON3
    using Logging

    library=MaterialsLibrary(add_defaults = false)
    copper=Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
    add!(library, "copper", copper)

    mktempdir() do directory
        requested=joinpath(directory, "materials.dat")
        saved=@test_logs (:warn, r"Forcing extension") (:info, r"saved") save(
            library;
            file_name = requested
        )
        @test saved == joinpath(directory, "materials.json")
        @test isfile(saved)

        restored=MaterialsLibrary(add_defaults = false)
        @test_logs (:info, r"Loading materials") (:info, r"Successfully loaded") load!(
            restored;
            file_name = saved
        )
        @test collect(keys(restored)) == ["copper"]
        @test restored["copper"] == copper

        unsupported=joinpath(directory, "materials.txt")
        write(unsupported, "not json")
        @test_logs (:error, r"only supports") load!(restored; file_name = unsupported)
        @test haskey(restored, "copper")
        @test_throws ErrorException load!(
            restored;
            file_name = joinpath(directory, "missing.json")
        )

        partial=joinpath(directory, "partial.json")
        open(partial, "w") do io
            JSON3.pretty(
                io,
                Dict(
                    "not-a-material"=>3,
                    "unknown-type"=>Dict("__julia_type__"=>"Missing.Type")
                )
            )
        end
        redirect_stderr(devnull) do
            @test_logs (:warn, r"Skipping material") min_level=Logging.Warn match_mode=:any load!(
                restored;
                file_name = partial
            )
        end
        @test isempty(restored)
    end
end

@testitem "ImportExport / XLSX / matrix topology and failure contract" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport
] begin
    using LinearAlgebra
    using XLSX

    frequencies=[50.0, 500.0]
    impedance=Array{ComplexF64}(undef, 2, 2, 2)
    admittance=similar(impedance)
    for k in eachindex(frequencies)
        impedance[:, :, k]=[1.0+2.0im 0.2+0.3im; 0.2+0.3im 1.5+2.5im]
        admittance[:, :, k]=[3.0+4.0im 0.4+0.5im; 0.4+0.5im 3.5+4.5im] .* 1e-6
    end
    parameters=LineParameters(impedance, admittance, frequencies)

    mktempdir() do directory
        output=export_data(:xlsx, parameters; file_name = joinpath(directory, "full.xlsx"))
        @test output == joinpath(directory, "full.xlsx")
        @test isfile(output)
        XLSX.openxlsx(output) do workbook
            @test Set(XLSX.sheetnames(workbook)) ==
                  Set(["Z(1,1)", "Z(1,2)", "Z(2,1)", "Z(2,2)",
                "Y(1,1)", "Y(1,2)", "Y(2,1)", "Y(2,2)"])
        end

        diagonal=LineParameters(
            cat(Diagonal([1.0+2.0im, 2.0+3.0im]); dims = 3),
            cat(Diagonal([3.0+4.0im, 4.0+5.0im]); dims = 3),
            [50.0]
        )
        diagonal_path=@test_logs (:warn, r"Z appears modal") (:warn, r"Y appears modal") export_data(
            :xlsx,
            diagonal;
            file_name = joinpath(directory, "diagonal.xlsx")
        )
        XLSX.openxlsx(diagonal_path) do workbook
            @test Set(XLSX.sheetnames(workbook)) ==
                  Set(["Z(1,1)", "Z(2,2)", "Y(1,1)", "Y(2,2)"])
        end

        @test_logs (:error, r"Failed to export XLSX") begin
            @test export_data(:xlsx, parameters; file_name = directory) === nothing
        end
    end
end

@testitem "ImportExport / TRALIN / export and parser contracts" tags=[:integration] setup=[
    ImportExportTestSupport,
    UseImportExportSupport,
    TestFixtures
] begin
    import LineCableModels.ImportExport: _block_after_anchor,
                                         _extract_tralin_frequencies, _infer_tralin_order,
                                         clean_variable_list, extract_tralin_variable,
                                         take_complex_list

    system=TestFixtures.three_phase_system()
    homogeneous=EarthModel([50.0, 500.0], 100.0, 10.0, 1.0)

    mktempdir() do directory
        exported=@test_logs (:info, r"TRALIN file saved") export_data(
            :tralin,
            system,
            homogeneous;
            freq = [50.0, 500.0],
            file_name = joinpath(directory, "case.f05")
        )
        @test basename(exported) == "tr_case.f05"
        text=read(exported, String)
        @test occursin("RUN-IDENTIFICATION,$(system.system_id)", text)
        @test occursin("UNIFORM,100.0,1,1", text)
        @test length(collect(eachmatch(r"FREQUENCY", text))) == 2
        @test count(line -> startswith(line, "GROUP,PH-"), readlines(exported)) == 3
        @test all(label -> occursin(label, text), ("CORE,", "SHEATH,", "ARMOUR,"))

        layered=EarthModel([50.0], 80.0, 8.0, 1.0; t = 4.0)
        add!(layered, [50.0], 300.0, 12.0, 1.0; t = Inf)
        layered_path=export_data(
            :tralin,
            system,
            layered;
            file_name = joinpath(directory, "layered.f05")
        )
        layered_text=read(layered_path, String)
        @test occursin("HORIZONTAL", layered_text)
        @test occursin("LAYER,TOP,80.0,4.0,1.0,8.0", layered_text)
        @test occursin("LAYER,BOTTOM,300.0,,1.0,12.0", layered_text)

        oversized_design=deepcopy(TestFixtures.mv_cable_design())
        extra=deepcopy(last(oversized_design.components))
        push!(oversized_design.components, extra)
        oversized_position=CablePosition{Float64}(
            oversized_design,
            0.0,
            -1.0,
            [1, 0, 0, 0]
        )
        oversized_system=LineCableSystem("oversized", 1000.0, oversized_position)
        @test_throws ArgumentError export_data(
            :tralin,
            oversized_system,
            homogeneous;
            file_name = joinpath(directory, "oversized.f05")
        )

        short_position=CablePosition{Float64}(
            TestFixtures.mv_cable_design(),
            0.0,
            -1.0,
            [1]
        )
        short_system=LineCableSystem("short-mapping", 1000.0, short_position)
        @test_throws ArgumentError export_data(
            :tralin,
            short_system,
            homogeneous;
            file_name = joinpath(directory, "short.f05")
        )

        report=String[
            "CHARACTERISTICS OF ALL CONDUCTORS",
            " 1 1 1 1 1 core 0.0 0.1",
            "TRALIN package - PAGE 2",
            "FREQUENCY OF HARMONIC CURRENT:",
            " 1 50.0",
            "TRALIN package - PAGE 3",
            "GROUND WIRES ELIMINATED",
            "POTENTIAL COEFFICIENTS (meghoms.kilometer)"
        ]
        append!(report, fill("header", 14))
        push!(report, "1 7.0 + j 8.0")
        push!(report, "SERIES IMPEDANCES - (ohms/kilometer)")
        append!(report, fill("header", 14))
        push!(report, "1 1.0 + j 2.0")
        push!(report, "SHUNT ADMITTANCES (microsiemens/kilometer)")
        append!(report, fill("header", 14))
        push!(report, "1 3.0 + j 4.0")
        push!(report, "SERIES ADMITTANCES (siemens.kilometer)")

        report_path=joinpath(directory, "synthetic.f09")
        write(report_path, join(report, '\n'))
        @test _infer_tralin_order(report) == 1
        @test _infer_tralin_order(report_path) == 1
        @test _extract_tralin_frequencies(report) == [50.0]
        @test _block_after_anchor(report, "FREQUENCY")[1] == " 1 50.0"
        @test_throws ArgumentError _block_after_anchor(report, "missing")
        @test_throws ArgumentError _infer_tralin_order(["CHARACTERISTICS OF ALL CONDUCTORS"])
        @test_throws ArgumentError _extract_tralin_frequencies(["FREQUENCY OF HARMONIC CURRENT:"])
        @test take_complex_list("4 1.5 + j 2.5E+1") == [4.0, 1.5 + 25.0im]
        @test isempty(take_complex_list("header"))
        @test clean_variable_list(Any[[1.0, 2.0 + 3.0im]], 2) ==
              [[2.0 + 3.0im, 0.0im], [0.0im, 0.0im]]
        redirect_stdout(devnull) do
            @test extract_tralin_variable(["no blocks"], 2, "start", "end") ==
                  zeros(ComplexF64, 2, 2)
        end

        parsed_frequencies, parsed_Z,
        parsed_Y, parsed_P=LineCableModels.ImportExport.parse_tralin_file(report_path)
        @test parsed_frequencies == [50.0]
        @test parsed_Z[1, 1, 1] == (1.0 + 2.0im) / 1000
        @test parsed_Y[1, 1, 1] ≈ (3.0 + 4.0im) * 1e-9
        @test parsed_P[1, 1, 1] == (7.0 + 8.0im) * 1e9

        parameters=LineParameters(report_path)
        @test frequencies(parameters) == [50.0]
        symbol_parameters=LineParameters(:tralin, report_path)
        @test Z(symbol_parameters) == Z(parameters)
        @test Y(symbol_parameters) == Y(parameters)
        @test_throws ArgumentError LineParameters(report_path; format = :unknown)
        @test_throws ArgumentError LineParameters(Val(:unsupported), report_path)
    end
end
