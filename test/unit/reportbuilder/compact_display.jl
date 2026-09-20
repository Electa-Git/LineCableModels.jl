@testitem "ReportBuilder / compact views preserve values, missing masks and quantities" tags=[:unit] begin
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    using DataFrames
    f=[0.1,10.,50.,1e3,1e6,1e7]
    reference=LineParameters(fill(1.0+1im,1,1,6),zeros(ComplexF64,1,1,6),f;
        details=ComputationDetails(;coordinates=["core"],))
    candidate=LineParameters(fill(1.2+1im,1,1,6),zeros(ComplexF64,1,1,6),f;
        details=ComputationDetails(;coordinates=["core"],))
    artifact=report(BenchmarkTableDefinition((R,G);bands=(:all,:dc,:harmonic,:narrow,:wide)),(;reference,candidate))
    before=deepcopy(artifact.tables)
    errors=deepcopy(artifact.observed.errors)
    for mime in (MIME"text/plain"(),MIME"text/html"())
        text=sprint(show,mime,artifact)
        @test occursin("R — Relative RMS [%]",text)
        @test occursin("G — Relative RMS [%]",text)
        @test occursin("missing",text)
        @test occursin("eligible terms only",text)
        @test occursin("No controlled performance measurements",text)
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
