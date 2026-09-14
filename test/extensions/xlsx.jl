@testitem "Extensions / XLSX boundary / unloaded writer is absent" tags = [
    :extension,
    :core_only
] begin
    import LineCableModels

    @test Base.get_extension(LineCableModels, :LineCableModelsXLSXExt) === nothing

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

    # Encoding belongs to core; writing exists only when the XLSX extension
    # loads. write is not a mandatory hook for every report definition.
    @test !applicable(report_builder.write,
        definition, parameters, selected, table, nothing, encoded)
    @test_throws MethodError LineCableModels.report(
        definition,
        parameters
    )

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
    # Exercise encoding, first-sheet rename, subsequent sheets and readback
    # through the current XLSX writer with distinguishable channel values.
    parameters=LineCableModels.LineParameters(
        reshape(ComplexF64[1+2im,3+4im],1,1,2),
        reshape(ComplexF64[5e-6+6e-6im,7e-6+8e-6im],1,1,2),[50.,125.])
    mktempdir() do directory
        cd(directory) do
            LineCableModels.report(LineCableModels.XLSXReportDefinition(),parameters)
            workbook=XLSX.readxlsx("ZY_export.xlsx")
            @test XLSX.sheetnames(workbook)==["Z(1,1)","Y(1,1)"]
            @test workbook["Z(1,1)"]["A6:C7"]==["50" "1000" "2000"; "125" "3000" "4000"]
            @test workbook["Y(1,1)"]["A6"]=="50"
            @test workbook["Y(1,1)"]["B6"]!=workbook["Z(1,1)"]["B6"]
        end
    end
end
