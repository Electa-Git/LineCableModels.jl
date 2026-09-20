@testitem "Extensions / XLSX boundary / unloaded writer is absent" tags=[:extension,:core_only] begin
    const RB=LineCableModels.ReportBuilder
    @test Base.get_extension(LineCableModels,:LineCableModelsXLSXExt)===nothing
    observed=ObservedResult(LineParameters(fill(1.0+2im,2,2,1),fill(3.0+4im,2,2,1),[50.]))
    definition=XLSXReportDefinition()
    tables=RB.tabulate(definition,observed)
    encoded=RB.encode(definition,observed,tables,nothing)
    @test length(encoded)==4
    @test all(book -> book isa RB.XLSXWorkbook,encoded)
    @test first(encoded).sheets[1].cells[2,2]==1000
    @test getproperty.(first(encoded).sheets,:name)==["values","std","metadata"]
    @test !applicable(RB.write,definition,encoded)
    @test_throws MethodError report(definition,observed)
end

@testitem "Extensions / XLSX writer / explicit activation" tags=[:extension] begin
    using XLSX
    extension=Base.get_extension(LineCableModels,:LineCableModelsXLSXExt)
    @test extension!==nothing
    @test any(method -> method.module===extension,methods(LineCableModels.ReportBuilder.write))
    parameters=LineParameters(fill(1.0+2im,1,1,1),fill(3e-6+4e-6im,1,1,1),[50.])
    mktempdir() do directory
        files=export_data(:xlsx,parameters;file_name=joinpath(directory,"owned.xlsx"))
        @test length(files)==4
        r=XLSX.readxlsx(only(filter(path -> endswith(path,"_R.xlsx"),files)))
        g=XLSX.readxlsx(only(filter(path -> endswith(path,"_G.xlsx"),files)))
        @test r["values"]["B2"]==1000
        @test g["values"]["B2"]≈0.003
        @test r["values"]["B2"] isa Number
    end
end
