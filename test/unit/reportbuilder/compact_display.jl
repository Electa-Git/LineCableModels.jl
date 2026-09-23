@testitem "ReportBuilder / compact views preserve values, missing masks and quantities" tags=[:unit] begin
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    using DataFrames
    f=[0.1,10.,50.,1e3,1e6,1e7]
    reference=LineParameters(fill(1.0+1im,1,1,6),zeros(ComplexF64,1,1,6),f;
        details=ComputationDetails(;coordinates=["core"],))
    candidate=LineParameters(fill(1.2+1im,1,1,6),zeros(ComplexF64,1,1,6),f;
        details=ComputationDetails(;coordinates=["core"],))
    artifact=report(BenchmarkTableDefinition((R,G);bands=(:all,:dc,:harmonic,:narrow,:wide)),(;reference,candidate))
    @test artifact[R]===artifact.tables.quantities[1].Z.R
    @test artifact[1,R]===artifact[R]
    @test artifact.reference!==nothing
    original_id=artifact.observed.gridpoint.id
    reference_id=artifact.reference.gridpoint.id
    @test all(row -> row.candidate_id==original_id && row.reference_id==reference_id,
        artifact.observed.errors)
    before=deepcopy(artifact.tables)
    errors=deepcopy(artifact.observed.errors)
    for mime in (MIME"text/plain"(),MIME"text/html"())
        text=sprint(show,mime,artifact)
        @test occursin("R — Relative RMS [%]",text)
        @test occursin("G — Relative RMS [%]",text)
        @test occursin("missing",text)
        @test occursin("eligible terms only",text)
        @test occursin("No controlled performance measurements",text)
        @test !occursin("Candidate 1",text)
        @test !occursin("configuration 1",text)
        @test !occursin("Absolute RMS",text)
        @test !occursin("<svg",text) && !occursin("<img",text)
        for feature in artifact.tables.features
            rendered=sprint((io,frame) -> show(IOContext(io,:limit=>false),MIME"text/html"(),
                frame;summary=false,eltypes=false),feature.relative)
            mime isa MIME"text/html" && @test occursin(rendered,text)
        end
        absolute=sprint((io,value) -> show(io,mime,value;metric=:absolute),artifact)
        @test occursin("Absolute RMS [Ω/m]",absolute)
        @test occursin("Absolute RMS [S/m]",absolute)
        @test_throws ArgumentError show(IOBuffer(),mime,artifact;metric=:bogus)
        @test_throws ArgumentError show(IOBuffer(),mime,artifact;problem=2)
    end
    # These compare actual retained values and availability masks, not just table
    # shapes: formatting must not clip noise, replace missing, or normalize again.
    for (prior,after) in zip(before.features,artifact.tables.features)
        @test isequal(prior.relative,after.relative)
        @test isequal(prior.absolute,after.absolute)
        @test propertynames(after.relative)==[:formula,:all,:dc,:harmonic,:narrow,:wide]
    end
    @test all(ismissing,Matrix(only(filter(f -> f.quantity===:G,artifact.tables.features)).relative[:,2:end]))
    @test all(value -> value≈20,Matrix(only(filter(f -> f.quantity===:R,artifact.tables.features)).relative[:,2:end]))
    @test isequal(before.terms,artifact.tables.terms)
    @test isequal(before.maxima,artifact.tables.maxima)
    @test artifact.observed.gridpoint.id==original_id
    @test artifact.reference.gridpoint.id==reference_id
    @test artifact[R]===artifact.tables.quantities[1].Z.R
    for (a,b) in zip(errors,artifact.observed.errors)
        @test isequal(a.relative,b.relative)
        @test isequal(a.absolute,b.absolute)
        @test isequal(a.settings,b.settings)
    end

    # Descriptions are user-extensible text. HTML must escape them without using
    # rendered text as a scientific key or altering the underlying formula table.
    feature=first(artifact.tables.features)
    feature.relative.formula[1]="test <formula> & option"
    html=sprint(show,MIME"text/html"(),artifact)
    @test occursin("test &lt;formula&gt; &amp; option",html)
    @test !occursin("test <formula>",html)
end

@testitem "ReportBuilder / ordinary scientific display and bounded previews" tags=[:unit] setup=[TestFixtures] begin
    using DataFrames, Measurements, Statistics
    const RB=LineCableModels.ReportBuilder
    const U=LineCableModels.Units
    constants=CableConstants(1e-4,2e-7,3e-10,4e-12;frequency=50)
    r=report(constants;values=(R,L,G,C),length_unit=:kilo,
        quantity_units=(R=:base,L=:milli,G=:micro,C=:micro))
    original=deepcopy(r.observed.quantities)
    originals=map(copy,(r[R],r[L],r[G],r[C]))
    normal=IOContext(IOBuffer(),:limit=>true,:displaysize=>(24,80),:color=>false)
    for mime in (MIME"text/plain"(),MIME"text/html"())
        rendered=sprint(show,mime,r;context=normal)
        for (request,value) in zip((R,L,G,C),("0.1","0.2","0.004","0.3"))
            @test occursin(value,rendered)
            @test occursin(U.label(metadata(r[request],"unit")),rendered)
        end
        @test occursin("50.0",rendered)
        @test occursin("core",rendered)
        @test occursin("Hz",rendered)
        @test !occursin("ReportArtifact(",rendered)
        @test !occursin("illustration",rendered)
        @test !occursin(r"[Cc]andidate|Result 1|Gridpoint 1",rendered)
        @test sprint(show,mime,r;context=:compact=>true)==sprint(show,r)
        mime isa MIME"text/html" && @test length(collect(eachmatch(r"<table",rendered)))==4
    end
    @test occursin("4 tables",sprint(summary,r))
    @test sprint(summary,r)=="Report · 4 tables"
    @test sprint(show,r)=="ReportArtifact(tables=4)"
    @test !occursin('\n',sprint(show,r))
    @test isequal(r.observed.quantities,original)
    @test all(isequal(a,b) for (a,b) in zip(originals,(r[R],r[L],r[G],r[C])))

    design=TestFixtures.coaxial_design()
    completed=CableConstants(design)
    single=report(completed;values=R)
    owner_label=only(LineCableModels.Grammar.observation_labels(single.observed))
    @test sprint(summary,single)=="Report · 1 table"
    @test owner_label==description(CableConstantsFormulation();compact=true)
    other_observed=ObservedResult(CableConstants(design;temperature=40.))
    with_reference=report(single.observed;values=R,reference=other_observed)
    comparison_labels=LineCableModels.Grammar.observation_labels([single.observed,other_observed])
    single_table=single[R]
    table_before=copy(single_table)
    for artifact in (single,with_reference), mime in (MIME"text/plain"(),MIME"text/html"())
        rendered=sprint(show,mime,artifact;context=normal)
        @test occursin(artifact===single ? owner_label : first(comparison_labels),rendered)
        @test !occursin(r"[Cc]andidate|Result 1|Gridpoint 1",rendered)
        @test artifact[1,R]===artifact[R]
        @test isequal(artifact.observed.quantities,single.observed.quantities)
        if artifact===with_reference
            @test occursin("Reference",rendered)
            @test occursin(last(comparison_labels),rendered)
        end
    end
    @test single[R]===single_table
    @test isequal(single_table,table_before)
    @test with_reference.reference===other_observed
    @test with_reference.observed===single.observed

    points=[single.observed,other_observed]
    study=report(points;values=R)
    labels=LineCableModels.Grammar.observation_labels(points)
    @test sprint(summary,study)=="Report · 2 gridpoints · 2 tables"
    for mime in (MIME"text/plain"(),MIME"text/html"())
        rendered=sprint(show,mime,study;context=normal)
        @test all(label -> occursin(label,rendered),labels)
        @test !occursin(r"[Cc]andidate|\[1\]|\[2\]",rendered)
    end

    z=reshape(complex.(1.:800.,801.:1600.),2,2,200)
    source=LineParameters(z,fill(3e-6+4e-6im,2,2,200),collect(1.:200.))
    collection=report([source,source];values=(R,L,G,C),frequency_unit=:kilo)
    for size in ((24,80),(40,120),(10,40),(5,30),(2,8),(1,1))
        context=IOContext(IOBuffer(),:limit=>true,:displaysize=>size,:color=>false)
        shown=sprint(show,MIME"text/plain"(),collection;context)
        @test length(split(shown,'\n'))<=size[1]
        @test all(textwidth(line)<=size[2] for line in split(shown,'\n'))
        @test !endswith(shown,'\n')
        @test !isempty(shown)
    end
    preview=sprint(show,MIME"text/plain"(),collection;context=normal)
    @test occursin("kHz",preview)
    @test occursin("omitted",preview)
    @test length(collect(eachmatch(r"\b0\.001\b",preview)))==4
    for request in (R,L,G,C)
        @test occursin(U.label(metadata(collection[1,request],"quantity")),preview)
    end
    unlimited=IOContext(IOBuffer(),:limit=>false,:displaysize=>(2,8),:color=>false)
    for mime in (MIME"text/plain"(),MIME"text/html"())
        shown=sprint(show,mime,collection;context=unlimited)
        @test occursin("800000",shown) || occursin("800000.0",shown) || occursin("800000.",shown) || occursin("800000e",shown) || occursin("8.0e5",shown)
        @test occursin("[1,2]",shown) && occursin("[2,1]",shown)
        @test occursin("[2]",shown)
        @test !occursin(r"[Cc]andidate",shown)
        @test !occursin("omitted",shown)
        mime isa MIME"text/html" && @test length(collect(eachmatch(r"<table",shown)))==8
    end
    selected=report(source;values=@observe(R[2,1,199:200]),frequency_unit=:kilo)
    selected_html=sprint(show,MIME"text/html"(),selected;context=unlimited)
    @test length(collect(eachmatch(r"<table",selected_html)))==1
    @test occursin("[2,1]",selected_html) && !occursin("[1,2]",selected_html)

    labeled=ObservedResult(constants;gridpoint=(id=nothing,name="Cable <core> & sheath",))
    html=sprint(show,MIME"text/html"(),report(labeled))
    @test occursin("Cable &lt;core&gt; &amp; sheath",html)
    @test !occursin("Cable <core>",html)
    bare=ReportArtifact(r.observed,nothing,DataFrame(value=[1.,2.]),nothing,nothing)
    @test occursin("<table",sprint(show,MIME"text/html"(),bare))
    @test occursin("2.0",sprint(show,MIME"text/plain"(),bare))
    aggregate=ReportArtifact([r.observed,r.observed],nothing,bare.tables,nothing,nothing)
    @test occursin("2 gridpoints",sprint(summary,aggregate))
    @test occursin("1 table",sprint(summary,aggregate))
    @test_throws r"not associated with a gridpoint" aggregate[R]
    @test_throws r"not associated with a gridpoint" aggregate[2,R]
    @test length(collect(eachmatch(r"<table",sprint(show,MIME"text/html"(),aggregate))))==1
    # No terminal width may turn one integer into a different, clipped token.
    large_integer=parse(Int128,"123456789012345678901")
    wide_table=DataFrame(frequency=[50],coefficient=[large_integer],other=[77])
    wide=ReportArtifact(r.observed,nothing,wide_table,nothing,nothing)
    for width in 12:2:72
        context=IOContext(IOBuffer(),:limit=>true,:displaysize=>(12,width),:color=>false)
        text=sprint(show,MIME"text/plain"(),wide;context)
        @test all(textwidth(line)<=width for line in split(text,'\n'))
        @test all(match.match in ("1","3","50",string(large_integer),"77")
            for match in eachmatch(r"[0-9]+",text))
        width==72 && @test occursin(string(large_integer),text)
    end
    for table in (DataFrame(value=Float64[]),DataFrame(value=[missing,missing]))
        artifact=ReportArtifact(r.observed,nothing,table,nothing,nothing)
        for mime in (MIME"text/plain"(),MIME"text/html"())
            @test !isempty(sprint(show,mime,artifact;context=normal))
        end
    end
    uq=report(TestFixtures.cable_monte_carlo_result();values=((statistics,R,mean),(statistics,R,std)))
    for mime in (MIME"text/plain"(),MIME"text/html"())
        shown=sprint(show,mime,uq;context=normal)
        @test occursin("mean",shown) && occursin("std",shown)
    end
    shared=measurement(1e-4,1e-6)
    uncertain=report(CableConstants(shared,2shared,3shared,4shared);values=R,clip=false)
    bigsource=LineParameters(fill(complex(-big"1.23456789e-40",big"0"),1,1,2),
        fill(complex(big"0",big"1e-6"),1,1,2),BigFloat[50,500])
    precise=report(bigsource;values=R,clip=false)
    for artifact in (uncertain,precise)
        table=artifact[R]
        before=deepcopy(table)
        quantities=deepcopy(artifact.observed.quantities)
        for mime in (MIME"text/plain"(),MIME"text/html"())
            shown=sprint(show,mime,artifact;context=unlimited)
            @test occursin(artifact===uncertain ? "±" : "-",shown)
        end
        @test isequal(table,before)
        @test isequal(artifact.observed.quantities,quantities)
    end
    @test Base.nonmissingtype(eltype(precise[R][!,2]))===BigFloat
end
