@testitem "Extensions / XLSX boundary / unloaded writer is absent" tags = [
    :extension,
    :core_only
] begin
    import LineCableModels
    using RequiredInterfaces: NotImplementedError

    @test Base.get_extension(LineCableModels, :LineCableModelsXLSXExt) === nothing
    @test !isdefined(LineCableModels.ReportBuilder, :_write_xlsx_sheet!)

    impedance = reshape(
        ComplexF64[1 + 2im, 0.2 + 0.3im, 0.2 + 0.3im, 1.5 + 2.5im],
        2,
        2,
        1
    )
    admittance = 1.0e-6 .* reshape(
        ComplexF64[3 + 4im, 0.4 + 0.5im, 0.4 + 0.5im, 3.5 + 4.5im],
        2,
        2,
        1
    )
    parameters = LineCableModels.LineParameters(impedance, admittance, [50.0])
    definition = LineCableModels.XLSXReportDefinition()
    report_builder = LineCableModels.ReportBuilder
    selected = report_builder.select(definition, parameters)
    table = report_builder.tabulate(definition, parameters, selected)
    encoded = report_builder.encode(
        definition,
        parameters,
        selected,
        table,
        nothing
    )

    @test encoded isa report_builder.XLSXWorkbook
    @test encoded.destination == joinpath(pwd(), "ZY_export.xlsx")
    @test getproperty.(encoded.sheets, :name) == [
        "Z(1,1)", "Z(1,2)", "Z(2,1)", "Z(2,2)",
        "Y(1,1)", "Y(1,2)", "Y(2,1)", "Y(2,2)"
    ]
    @test encoded.sheets[1].cells == ["frequency" "Hz" "";
           "R" "Ω/km" "";
           "X" "Ω/km" "";
           "" "" "";
           "frequency" "R" "X";
           "50" "1000" "2000"]
    @test report_builder.encode_cell(definition, missing) == ""
    @test report_builder.encode_cell(definition, 1 / 3) == "0.333333333333"

    @test_throws NotImplementedError LineCableModels.report(
        definition,
        parameters
    )
    for private_name in (
        :_xlsx_string,
        :_xlsx_strings,
        :_xlsx_units,
        :_xlsx_destination,
        :_write_xlsx_sheet!
    )
        @test !isdefined(report_builder, private_name)
    end
end

@testitem "Extensions / XLSX writer / explicit package activation" tags = [
    :extension
] begin
    using XLSX
    import LineCableModels

    extension_module = Base.get_extension(
        LineCableModels,
        :LineCableModelsXLSXExt
    )
    @test extension_module !== nothing
    @test any(
        method -> method.module === extension_module,
        methods(LineCableModels.ReportBuilder.write)
    )
    @test !isdefined(extension_module, :_write_xlsx_sheet!)
end
