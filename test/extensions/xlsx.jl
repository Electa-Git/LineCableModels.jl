@testitem "Extensions / XLSX boundary / unloaded writer is absent" tags = [
    :extension,
    :core_only
] begin
    import LineCableModels

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

    @test_throws MethodError LineCableModels.report(
        LineCableModels.XLSXReport(),
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
end
