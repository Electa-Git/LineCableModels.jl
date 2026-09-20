@testitem "ReportBuilder / observations and fixed stage order" tags=[:unit] begin
    const RB=LineCableModels.ReportBuilder
    source=LineParameters(fill(1.0+2im,1,1,2),fill(3.0+4im,1,1,2),[1.,2.])
    observed=ObservedResult(source)
    const stages=Symbol[]
    struct StageReport <: RB.AbstractReportDefinition end
    struct Tabulated
        observed::ObservedResult
    end
    struct Illustrated
        tables::Tabulated
    end
    struct Encoded
        illustration::Illustrated
    end
    RB.tabulate(::StageReport,point::ObservedResult;reference=nothing) = begin
        push!(stages,:tabulate); Tabulated(point)
    end
    RB.illustrate(::StageReport,point,tables::Tabulated;reference=nothing) = begin
        @test tables.observed===point
        push!(stages,:illustrate); Illustrated(tables)
    end
    RB.encode(::StageReport,point,tables,illustration::Illustrated;reference=nothing) = begin
        @test illustration.tables===tables && tables.observed===point
        push!(stages,:encode); Encoded(illustration)
    end
    RB.write(::StageReport,encoded::Encoded) = begin
        push!(stages,:write);encoded
    end
    artifact=report(StageReport(),observed)
    @test stages==[:tabulate,:illustrate,:encode,:write]
    @test artifact.observed===observed
    @test artifact.reference===nothing
    @test artifact.tables.observed===observed
    @test artifact.output.illustration===artifact.illustration
    @test fieldnames(typeof(artifact))==(:observed,:reference,:tables,:illustration,:output)
    @test_throws TypeError report(TableReportDefinition(),observed;reference=source)
    @test_throws MethodError report(StageReport(),source)
    illustrated=Ref{Any}(nothing)
    definition=TableReportDefinition((R,);illustration=(point;ydata) -> (illustrated[]=point))
    artifact=report(definition,observed)
    @test illustrated[]===observed
    @test artifact.tables.R[!,2]==[1000.,1000.]
    @test report(TableReportDefinition(),source).observed isa ObservedResult
end

@testitem "ReportBuilder / XLSX numeric quantities and complete matrices" tags=[:integration] begin
    using XLSX, Measurements, LinearAlgebra
    const RB=LineCableModels.ReportBuilder
    impedance=reshape(ComplexF64.(1:8),2,2,2)
    impedance[1,2,2]=0.25+0.05im
    impedance[2,1,2]=0.5+0.1im
    parameters=LineParameters(impedance,fill(3e-6+4e-6im,2,2,2),[50.,500.])
    mktempdir() do directory
        artifact=report(XLSXReportDefinition(file_name=joinpath(directory,"full.xlsx")),parameters)
        @test length(artifact.output)==4
        @test all(isfile,artifact.output)
        resistance=only(filter(path -> endswith(path,"_R.xlsx"),artifact.output))
        XLSX.openxlsx(resistance) do workbook
            @test XLSX.sheetnames(workbook)==["values","std","metadata"]
            sheet=workbook["values"]
            @test sheet["A1"]=="frequency"
            @test sheet["C1"]=="[1,2]"
            @test sheet["D1"]=="[2,1]"
            @test sheet["A2"]==50.0
            @test sheet["C3"]==250.0
            @test sheet["D3"]==500.0
            @test sheet["C3"] isa Number
        end
        @test length(export_data(:xlsx,parameters;file_name=joinpath(directory,"delegated.xlsx")))==4
        diagonal=LineParameters(cat(Diagonal([1.0+2im,2.0+3im]);dims=3),
            cat(Diagonal([3.0+4im,4.0+5im]);dims=3),[50.])
        result=report(XLSXReportDefinition(file_name=joinpath(directory,"diagonal.xlsx")),diagonal)
        workbook=XLSX.readxlsx(only(filter(path -> endswith(path,"_R.xlsx"),result.output)))
        @test workbook["values"]["C2"]==0
        @test workbook["values"]["D2"]==0
        shared=measurement(1e-4,1e-5)
        uncertain=LineParameters(fill(complex(shared,2shared),1,1,2),fill(3e-6+4e-6im,1,1,2),[50.,500.])
        result=report(XLSXReportDefinition(file_name=joinpath(directory,"uncertain.xlsx")),uncertain)
        workbook=XLSX.readxlsx(only(filter(path -> endswith(path,"_R.xlsx"),result.output)))
        @test workbook["values"]["B2"]≈0.1
        @test workbook["std"]["B2"]≈0.01
        @test workbook["std"]["B2"] isa Number
        @test impedance==Z(parameters)
    end
    definition=XLSXReportDefinition()
    @test_throws ArgumentError RB.encode_cell(definition,big"1e400")
    @test_throws ArgumentError RB.encode_cell(definition,big"1e-400")
    @test RB.encode_cell(definition,1/3)==1/3
    @test ismissing(RB.encode_cell(definition,missing))
    @test Base.ispublic(RB,:XLSXSheet)
    @test Base.ispublic(RB,:XLSXWorkbook)
end

@testitem "ReportBuilder / completed result conveniences use observations" tags=[:unit] setup=[TestFixtures] begin
    using DataFrames, Measurements, Tables
    const RB=LineCableModels.ReportBuilder
    constants=CableConstants(1.,2.,3.)
    artifact=report(RB.CableConstantsTableDefinition(),constants)
    @test artifact.observed isa ObservedResult
    @test artifact.tables.constants.R.value==[1000.]
    @test !Tables.istable(typeof(artifact.observed))
    @test parentmodule(which(DataFrame,(ObservedResult,)))===RB
    mc=TestFixtures.cable_monte_carlo_result()
    artifact=report(RB.MonteCarloTableDefinition(),mc)
    @test artifact.observed isa Vector{ObservedResult}
    @test length(artifact.tables)==length(mc)
    @test length(first(artifact.observed).quantities)==8
    @test !Tables.istable(typeof(mc))
    for name in (:TableReportDefinition,:CableConstantsTableDefinition,:LineParametersTableDefinition,
            :BenchmarkTableDefinition,:MonteCarloTableDefinition,:XLSXReportDefinition)
        @test parentmodule(getproperty(RB,name))===RB
    end
end
