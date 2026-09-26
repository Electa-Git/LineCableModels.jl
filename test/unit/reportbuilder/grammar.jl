@testitem "ReportBuilder / observations and fixed stage order" tags=[:unit] begin
    const RB=LineCableModels.ReportBuilder
    source=LineParameters(fill(1.0+2im,1,1,2),fill(3.0+4im,1,1,2),[1.,2.])
    observed=ObservedResult(source)
    const stages=Symbol[]
    struct StageReport <: RB.AbstractReportDefinition end
    struct Selected
        observed::ObservedResult
    end
    struct Tabulated
        observed::ObservedResult
    end
    struct Illustrated
        tables::Tabulated
    end
    struct Encoded
        illustration::Illustrated
    end
    RB.select(::StageReport,point::ObservedResult;reference=nothing) = begin
        push!(stages,:select); Selected(point)
    end
    RB.tabulate(::StageReport,point::ObservedResult,selected::Selected;reference=nothing) = begin
        @test selected.observed===point
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
    @test stages==[:select,:tabulate,:illustrate,:encode,:write]
    empty!(stages)
    @test RB.tabulate(StageReport(),observed).observed===observed
    @test stages==[:select,:tabulate]
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
    @test artifact.tables.Z.R[!,2]==[1000.,1000.]
    @test RB.tabulate(definition,observed)==artifact.tables
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
    @test artifact.tables.constants.R[!,2]==[1000.]
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

@testitem "ReportBuilder / XLSX destinations are preflighted before writes" tags=[:integration] begin
    using XLSX
    RB=LineCableModels.ReportBuilder
    observed=ObservedResult(LineParameters(fill(1.0+2im,1,1,2),fill(3.0+4im,1,1,2),[1.,2.]))
    mktempdir() do directory
        definition=XLSXReportDefinition(file_name=joinpath(directory,"safe.xlsx"))
        books=RB.encode(definition,observed,RB.tabulate(observed),nothing)
        open(last(books).destination,"w") do io
            write(io,"retained destination")
        end
        @test_throws r"already exists" RB.write(definition,books)
        @test !ispath(first(books).destination)
        @test read(last(books).destination,String)=="retained destination"
        @test_throws r"already exists" report(definition,observed)
        @test length(readdir(directory))==1
        overwrite=XLSXReportDefinition(file_name=definition.file_name,overwrite=true)
        written=report(overwrite,observed).output
        @test length(written)==4
        @test XLSX.readxlsx(first(written))["values"]["B2"]==1000.
        @test length(export_data(:xlsx,observed;file_name=definition.file_name,overwrite=true))==4
        bad=RB.XLSXWorkbook(joinpath(directory,"oversized.xlsx"),[RB.XLSXSheet("values",Matrix{Any}(undef,0,16385))])
        @test_throws DimensionMismatch RB.write(definition,[RB.XLSXWorkbook(joinpath(directory,"first.xlsx"),first(books).sheets),bad])
        @test !ispath(joinpath(directory,"first.xlsx"))
        @test !ispath(bad.destination)
        @test_throws ArgumentError RB.write(overwrite,[first(books),first(books)])
    end
end

@testitem "ReportBuilder / source-first quantity tables" tags=[:unit] setup=[TestFixtures] begin
    using DataFrames, Statistics, Measurements
    using LineCableModels.Grammar: observation_product, observation_gridpoint, gridpoint_id
    using LineCableModels.Engine: retain_gridpoint
    RB=LineCableModels.ReportBuilder
    @test LineCableModels.report===RB.report
    @test isempty(Test.detect_ambiguities(RB; recursive = true))
    constants=CableConstants(1.0, 2.0, 3.0)
    line=TestFixtures.two_conductor_results()
    for selection in (nothing, ())
        @test keys(report(constants; values = selection).tables.constants)==(:R, :L, :G, :C)
        tables=report(line; values = selection).tables
        @test keys(tables.Z)==(:R, :X)
        @test keys(tables.Y)==(:G, :B)
    end
    @test report(constants).tables==report(constants; values = nothing).tables
    @test report(line).tables==report(line; values = ()).tables
    units=(length_unit = :kilo,
        quantity_units = (R = :base, L = :milli, G = :micro, C = :micro))
    requested=(R, L, G, C)
    artifact=report(constants; values = requested, units...)
    @test artifact.tables.constants.R[!, 2]==[1000.0]
    @test keys(artifact.tables)==(:constants,)
    @test names(artifact.tables.constants.R)==["frequency", string.(constants.cores)...]
    @test artifact.tables.constants.R.frequency==[constants.frequency]
    @test report(constants; values = @observe(R[1]), units...).tables.constants.R==artifact.tables.constants.R
    @test report(line, requested).tables==report(line; values = requested).tables
    @test keys(report(line; values = R).tables.Z)==(:R,)
    @test keys(report(line; values = R).tables)==(:Z,)
    @test report(line; values = ((Z, real),)).tables==report(line; values = R).tables
    @test keys(report(line; values = Z).tables.Z)==(:R, :X)
    full=report(line; values = requested)
    @test names(full.tables.Z.R)==["frequency", "[1,1]", "[1,2]", "[2,1]", "[2,2]"]
    @test full.tables.Z.R[!, 3]≈1000 .* observe(line, R)[1, 2, :]
    @test full.tables.Z.R[!, 4]≈1000 .* observe(line, R)[2, 1, :]
    @test full.tables.Z.R[!, 3]!=full.tables.Z.R[!, 4]
    selected=report(line; values = @observe(R[[2, 1], [2], [3, 1]]))
    @test names(selected.tables.Z.R)==["frequency", "[2,2]", "[1,2]"]
    @test selected.tables.Z.R.frequency==[1000.0, 10.0]
    @test selected.tables.Z.R[!, 2]≈1000 .* observe(line, R)[2, 2, [3, 1]]
    @test report(line.Z; values = R, frequencies = frequencies(line)).tables.Z.R==full.tables.Z.R
    @test report(line.Y; values = C, frequencies = frequencies(line)).tables.Y.C==full.tables.Y.C
    for bad in (:ydata, :rdata, :requests, :quantities)
        @test_throws r"values" report(line; NamedTuple{(bad,)}((requested,))...)
    end
    for source in (line, ObservedResult(line), [ObservedResult(line)])
        @test_throws r"values" report(source, R; values = R)
        @test_throws r"values" report(source, R; values = nothing)
        @test_throws ArgumentError report(source; frequency_unit = :base, freq_unit = :base)
        @test_throws ArgumentError report(source; backend = :cairo)
        @test_throws ArgumentError report(source; controls = false)
        @test_throws ArgumentError report(source; layout = (1, 1))
    end
    @test_throws MethodError report(42)
    @test_throws ArgumentError report(line; reference = [line])
    @test_throws ArgumentError report(ObservedResult[])
    source_id=gridpoint_id().source_id
    points=[retain_gridpoint(line, gridpoint_id(; source_id, problem_index = i))
            for i in (3, 1)]
    for collection in (points, Tuple(points))
        result=report(collection; values = R)
        @test getproperty.(getproperty.(result.observed, :gridpoint),
            :id)==[observation_gridpoint(point).id for point in points]
        @test [point.gridpoint.id.problem_index for point in result.observed]==[3, 1]
        @test length(result.tables)==2
        @test result.tables[1].Z.R==result.tables[2].Z.R
    end
    study=ParametricResult(
        nothing, points, (problems = [:one, :two], formulations = [:default]),
        ComputationDetails())
    @test report(study; values = R).tables==report(points; values = R).tables
    @test [point.gridpoint.id.problem_index for point in report(study).observed]==[3, 1]
    @test length(report(points[1:1]).observed)==1
    mixed=report((constants, line))
    @test keys(mixed.tables[1])==(:constants,)
    @test keys(mixed.tables[2])==(:Z, :Y)
    mc=TestFixtures.cable_monte_carlo_result()
    @test length(report(mc).tables)==length(mc)
    statistics_report=report(mc; values = @observe((statistics, R, mean)[1]))
    @test length(only(statistics_report.observed).quantities)==1
    @test only(statistics_report.tables).statistics.statistics_R_mean[!, 2]≈1000 .* observe(
        mc, statistics, R, mean, 1)
    @test length(only(report(mc; values = (statistics, R)).observed).quantities)==7
    @test report(mc; values = R).tables[1].constants.R==report(mc).tables[1].constants.R
    @test_throws ArgumentError DataFrame(artifact.observed)
end

@testitem "ReportBuilder / source-first retained units and illustration" tags=[:unit] setup=[TestFixtures] begin
    using DataFrames
    using LineCableModels.Grammar: observation_product
    line=TestFixtures.two_conductor_results()
    quantities=(R, L, G, C)
    unit_options=(length_unit = :base, quantity_units = :base, frequency_unit = :kilo)
    observed=ObservedResult(line, quantities; unit_options...)
    raw=report(line; values = quantities, unit_options...)
    retained=report(observed)
    @test retained.observed===observed
    @test raw.tables==retained.tables
    for selector in quantities
        a=observation_product(raw.observed, selector)
        b=observation_product(retained.observed, selector)
        @test a.values==b.values
        @test a.unit==b.unit
        @test a.coordinates==b.coordinates
    end
    @test report(observed; unit_options...).tables==retained.tables
    converted=report(observed; length_unit = :kilo, freq_unit = :base)
    roundtrip=report(converted.observed; length_unit = :base, frequency_unit = :kilo)
    for selector in quantities
        a=observation_product(roundtrip.observed, selector)
        b=observation_product(observed, selector)
        @test a.values≈b.values
        @test a.unit==b.unit && a.coordinates==b.coordinates
        @test a.available==b.available && a.engineering_zero==b.engineering_zero
    end
    @test report(observed; values = R).tables.Z.R==retained.tables.Z.R
    @test report(observed, R).tables==report(observed; values = R).tables
    @test only(report([observed], R).tables).Z.R==retained.tables.Z.R
    for source in (observed, [observed], (observed,)),
        options in (
            (clip = true,), (clip = nothing,), (atol = nothing,), (frequencies = nothing,))

        @test_throws ArgumentError report(source; options...)
    end
    recorded=Ref{Any}(nothing)
    renderer=(points; ydata, reference = nothing,
        kwargs...)->(recorded[]=(; points, ydata, reference, options = (; kwargs...)))
    @test report(line).illustration===nothing
    @test report(line).output===nothing
    @test recorded[]===nothing
    @test_throws ArgumentError report(line; illustration = renderer, unknown = true)
    @test_throws ArgumentError report(line; values = R, illustration = renderer, plot_options = (ydata = (L,),))
    @test_throws ArgumentError report(line; plot_options = (backend = :cairo,))
    @test recorded[]===nothing
    reference=ObservedResult(line; length_unit = :base)
    artifact=report(line; values = R, reference, illustration = renderer,
        plot_options = (ydata = (R,), backend = :cairo))
    @test recorded[].points===artifact.observed
    @test recorded[].reference===reference
    @test recorded[].ydata==((R, :, :, :),)
    @test recorded[].options==(backend = :cairo,)
    @test artifact.reference===reference
    @test isempty(artifact.observed.errors)
    @test isempty(artifact.reference.errors)
    @test observation_product(reference, R).unit!=observation_product(artifact.observed, R).unit
    converted_reference=report(line; values = R, reference, length_unit = :kilo)
    @test observation_product(converted_reference.reference,
        R).unit==observation_product(converted_reference.observed, R).unit
    @test observation_product(reference, R).unit!=observation_product(converted_reference.reference, R).unit
    @test report(line; values = R, reference = line).reference isa ObservedResult
    report(observed; values = R, illustration = renderer, length_unit = :kilo)
    @test recorded[].points isa ObservedResult
    @test recorded[].ydata==(R,)
    report(observed; illustration = renderer)
    @test recorded[].ydata==quantities
    @test raw.tables==retained.tables
    @test report(TableReportDefinition(quantities), line).tables==report(line; values = quantities).tables
    @test report(LineCableModels.ReportBuilder.LineParametersTableDefinition(quantities),
        line).tables==report(line; values = quantities).tables
end

@testitem "ReportBuilder / direct reported quantity access" tags=[:unit] setup=[TestFixtures] begin
    using DataFrames, Statistics, Measurements
    using LineCableModels.Grammar: gridpoint_id
    const RB=LineCableModels.ReportBuilder
    constants=CableConstants(1e-4,2e-7,3e-10,4e-12;frequency=50)
    r=report(constants;values=(R,L,G,C),length_unit=:kilo,
        quantity_units=(R=:base,L=:milli,G=:micro,C=:micro))
    for (request,expected) in zip((R,L,G,C),(0.1,0.2,0.004,0.3))
        table=r[request]
        @test table isa DataFrame
        @test table===getproperty(r.tables.constants,nameof(request))
        @test r[1,request]===table
        @test r[Int32(1),request]===table
        @test r[big(1),request]===table
        @test size(table)==(1,2)
        @test table.frequency==[50.]
        @test table[1,2]≈expected
        @test all(T -> T<:Number,eltype.(eachcol(table)))
        @test metadata(table,"request")===request
        @test metadata(table,"statistic")===:value
    end
    @test_throws BoundsError r[0,R]
    @test_throws BoundsError r[2,R]
    @test_throws ArgumentError r[true,R]
    @test_throws ArgumentError r[:R]
    @test_throws ArgumentError r["R"]
    @test_throws ArgumentError r[@observe(R[1])]

    selected=report(constants;values=@observe(R[1]))
    @test selected[R]===selected[@observe(R[1])]
    @test_throws ArgumentError selected[L]
    @test_throws ArgumentError selected[@observe(R[:])]
    @test_throws ArgumentError selected[@observe(R[2])]

    z=reshape(complex.(1.:8.,9.:16.),2,2,2)
    source=LineParameters(z,fill(3e-6+4e-6im,2,2,2),[50.,500.])
    full=report(source;values=(R,L,G,C))
    @test propertynames(full[R])==[:frequency,Symbol("[1,1]"),Symbol("[1,2]"),Symbol("[2,1]"),Symbol("[2,2]")]
    @test full[R][!,3]!=full[R][!,4]
    partial=report(source;values=@observe(R[2,1,2:2]))
    @test partial[R]===partial[@observe(R[2,1,2:2])]
    @test partial[R].frequency==[500.]
    @test metadata(partial[R],"coordinates").rows==[2]
    @test_throws ArgumentError partial[L]
    @test_throws ArgumentError partial[@observe(R[1,1,2:2])]

    points=[ObservedResult(constants;gridpoint=(id=gridpoint_id(problem_index=i),)) for i in (9,3)]
    collection=report(points;values=(R,L))
    tables=collection[R]
    @test tables isa Vector{DataFrame}
    @test all(tables[i]===collection[i,R] for i in 1:2)
    @test [metadata(table,"gridpoint_id").problem_index for table in tables]==[9,3]
    @test report(points[1:1];values=R)[R] isa Vector{DataFrame}
    @test report((constants,);values=R)[R] isa Vector{DataFrame}
    @test_throws BoundsError collection[9,R]
    incomplete=ReportArtifact(points,nothing,[collection.tables[1],(;L=collection[2,L])],nothing,nothing)
    @test_throws r"reported result 2.*problem_index = 3.*available reported requests" incomplete[R]
    diagnostic=try incomplete[R]; catch error; sprint(showerror,error); end
    @test occursin(repr(points[2].gridpoint.id),diagnostic)
    @test !occursin("candidate",diagnostic)
    ambiguous=ReportArtifact(r.observed,nothing,(a=r[R],b=copy(r[R])),nothing,nothing)
    @test_throws r"ambiguous.*reported result 1" ambiguous[R]
    undescribed=ReportArtifact(r.observed,nothing,DataFrame(R=[1.]),nothing,nothing)
    @test_throws r"absent quantity descriptors" undescribed[R]
    missing_table=copy(r[R]); allowmissing!(missing_table); missing_table[1,2]=missing
    missing_report=ReportArtifact(r.observed,nothing,missing_table,nothing,nothing)
    @test missing_report[R]===missing_table

    uq=report(TestFixtures.cable_monte_carlo_result();values=((statistics,R,mean),(statistics,R,std)))
    @test uq[1,(statistics,R,mean)]!==uq[1,(statistics,R,std)]
    @test metadata(uq[1,(statistics,R,mean)],"statistic")===:mean
    @test_throws ArgumentError uq[R]
    @test_throws ArgumentError uq[(statistics,L,mean)]

    original=deepcopy(r.observed.quantities)
    r[R][1,2]=42.
    @test isequal(r.observed.quantities,original)
    @test observe(constants,R)==[1e-4]
    @test r[R][1,2]==42.
end
