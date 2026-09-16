@testitem "Descriptions / owner dispatch survives composed reports and saved declarations" tags=[:unit] setup=[FormulaContractModels] begin
    using DataFrames, Statistics
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    import LineCableModels: description, formula_id, formulation_options
    E=LineCableModels.Engine
    IO=LineCableModels.ImportExport

    # Equal text must not become an identity, and a new leaf must reach the real
    # report without being added to a reader or plotting catalogue.
    struct TestEarthLeaf{ID} <: E.EarthImpedanceFormulation end
    formula_id(::TestEarthLeaf{ID}) where {ID} = ID
    description(::Type{<:TestEarthLeaf{ID}};compact::Bool=false) where {ID} =
        compact ? "Display-"*string(ID) : "Test-owned explanation shared by distinct equations"
    description(value::TestEarthLeaf;compact::Bool=false) = description(typeof(value);compact)
    formulation_options(::TestEarthLeaf) = FormulationOptions()
    Base.NamedTuple(value::TestEarthLeaf) = (identifier=formula_id(value),parameters=(;),options=(;))
    normal=Formulation()
    leaves=(TestEarthLeaf{:TestAlpha}(),TestEarthLeaf{:TestBeta}())
    choices=[E.LineParametersFormulation(merge(normal.methods,(earth_impedance=leaf,)),
        normal.options,merge(normal.definitions,(earth_impedance=leaf,))) for leaf in leaves]
    push!(choices,normal)
    @test description(LineCableModelsCoaxial())==description(LineCableModelsCoaxial)=="coaxial"
    @test description(LineCableModelsFEM())==description(LineCableModelsFEM)=="FEM"
    @test description(Formulation(:pscad))==description(PSCAD.PSCADFormulation)=="PSCAD"
    for value in (normal,LineCableModelsFEM(),Formulation(:pscad),MonteCarlo(normal),LinearError(normal))
        @test description(typeof(value))==description(value)
    end
    labels=description([MonteCarlo(normal),LinearError(normal)];roles=[:reference,:candidate])
    @test labels == ["Reference · Monte Carlo", "LEP"]
    @test all(label->!occursin("earth Z",label),labels)
    explicit_default=Formulation(earth_impedance=:unified)
    @test formula_id(normal,R)==formula_id(explicit_default,R)
    @test description([normal,explicit_default];quantity=R)==["coaxial","coaxial"]
    @test description([normal];quantity=R)==description([explicit_default];quantity=R)

    f=[0.1,50.,100.,1e3,1e6,1e7]
    z=reshape(complex.(collect(1.:24.),collect(101.:124.)),2,2,6)
    y=reshape(complex.(collect(201.:224.),collect(301.:324.)),2,2,6)*1e-6
    ports=["a","b"]
    ref=LineParameters(z,y,f;details=ComputationDetails(;coordinates=ports,formulations=NamedTuple(LineCableModelsFEM()),))
    points=[LineParameters(scale*z,y,f;details=ComputationDetails(;coordinates=ports,)) for scale in (1.1,1.2,1.3)]
    candidates=ParametricResult(nothing,points,(problems=[:one],formulations=choices), ComputationDetails((;)))
    source=(reference=ref,candidate=candidates)
    result=report(BenchmarkTableDefinition((R,X,G,B);bands=(:all,:dc,:harmonic,:narrow,:wide)),source)
    @test result.table.formulations.label[1]=="Reference · FEM"
    @test any(contains("TestAlpha"),result.table.formulations.label)
    @test any(contains("TestBeta"),result.table.formulations.label)
    @test count(contains("Test-owned explanation"),result.table.formula_details.selection)==2
    @test length(result.published.candidate.metadata.formulation_sources)==3
    @test result.published.candidate.metadata.formulation_sources[1]===choices[1]
    labels=description(choices[[3,1]];quantity=R)
    @test occursin("earth Z=Unified",labels[1]) && occursin("earth Z=Display-TestAlpha",labels[2])
    @test all(label->!occursin("internal Z",label),labels)
    @test all(feature -> !occursin("earth Y",join(feature.relative.formula)),
        filter(feature -> feature.quantity in (:R,:X),result.table.features))

    for native in (normal,Formulation(earth_impedance=:saad1996),LineCableModelsFEM(),
            Formulation(:pscad),MonteCarlo(normal),LinearError(normal))
        saved=IO.deserialize_value(Val(:formulation),NamedTuple(native))
        @test description([saved];roles=[:reference])==description([native];roles=[:reference])
        @test [(scope,formula_id(value)) for (scope,value) in pairs(saved...)]==
            [(scope,formula_id(value)) for (scope,value) in pairs(native)]
    end
    saved_pair=[IO.deserialize_value(Val(:formulation),NamedTuple(method)) for method in
        (MonteCarlo(normal),LinearError(normal))]
    # Ordinary collection promotion must not erase a retained method's owner.
    @test description(saved_pair;roles=[:reference,:candidate])==
        description([MonteCarlo(normal),LinearError(normal)];roles=[:reference,:candidate])
    for order in (:before,:after)
        native=Formulation(earth_impedance=formula(:carson1926;
            equivalent_earth=formula(:default;order)))
        saved=IO.deserialize_value(Val(:formulation),NamedTuple(native))
        label=only(description([native];quantity=R))
        @test occursin(string(order)*" FrequencyDependent",label)
        @test !occursin("FormulaDefinition{",label)
        @test description([saved];quantity=R)==[label]
        @test formula_id(saved,R)==formula_id(native,R)
    end
    # Historical inspection must not call today's declaration constructor or
    # reject controls that its current live owner would not admit.
    historical=NamedTuple(Formulation(earth_impedance=formula(:carson1926;
        equivalent_earth=formula(:bottommost;order=:before))))
    historical=merge(historical,(requested=merge(historical.requested,
        (earth_impedance=merge(historical.requested.earth_impedance,
            (equivalent_earth=(identifier=:ArchivedRule,order=:archived_order,
                parameters=(saved_parameter=2,),options=(saved_control=true,)),)),)),))
    saved=IO.deserialize_value(Val(:formulation),historical)
    label=only(description([saved];quantity=R))
    @test occursin("ArchivedRule archived_order FrequencyDependent",label)
    @test occursin("saved_control",label)
    # Inspection retains a native selection without executing its equation.
    selected_inner=FormulaContractModels.SurfaceLaw(kinds=(:inner,))
    routed=Formulation(
        internal_impedance=(inner=selected_inner,outer=:default,transfer=:default),
        earth_impedance=(air=:carson1926,earth=:pollaczek1926,mixed=:lucca1994),
        earth_admittance=formula(:default;parameters=(reference=:interface,),options=(integration=(method=:quad,),)))
    for native in (routed,MonteCarlo(routed),LinearError(routed),
            LineCableModelsFEM(options=(physics=:quasi_fw,)),
            Formulation(:pscad;options=(base_frequency=60.0,)))
        saved=IO.deserialize_value(Val(:formulation),NamedTuple(native))
        for quantity in (R,B)
            @test description([native,normal];quantity)==description([saved,normal];quantity)
        end
        @test [description(scope,value;compact=false) for (scope,value) in pairs(native)]==
            [description(scope,value;compact=false) for (scope,value) in pairs(saved...)]
    end
    routed_result=ParametricResult(nothing,points[1:2],
        (problems=[:one],formulations=[routed,normal]), ComputationDetails((;)))
    routed_report=report(BenchmarkTableDefinition((R,B);bands=(:all,)),
        (reference=ref,candidate=routed_result))
    @test any(contains("earth Z(air)=Carson"),first(routed_report.table.features).relative.formula)
    @test any(contains("inner"),routed_report.table.formula_details.selection)
    @test isempty(selected_inner.evaluations) && isempty(selected_inner.preparations)
    single=Formulation(insulation_admittance=FormulaContractModels.InsulationLaw())
    single_data=ParametricResult(nothing,points[1:1],
        (problems=[:one],formulations=[single]), ComputationDetails((;)))
    single_report=report(BenchmarkTableDefinition((R,B);bands=(:all,)),
        (reference=ref,candidate=single_data))
    single_label=only(last(single_report.table.features).relative.formula)
    @test occursin("insulation Y=InsulationLaw",single_label) && occursin("scale",single_label)
    @test !occursin("insulation Y",only(first(single_report.table.features).relative.formula))
    @test description([single];quantity=B)==description([
        IO.deserialize_value(Val(:formulation),NamedTuple(single))];quantity=B)
    @test isempty(selected_inner.evaluations) && isempty(selected_inner.preparations)
    fem_labels=description([LineCableModelsFEM(),LineCableModelsFEM(options=(physics=:quasi_fw,))];
        roles=[:reference,:reference])
    @test occursin("quasi-tem",first(fem_labels)) && occursin("quasi-fw",last(fem_labels))
    @test ismissing(IO.deserialize_value(Val(:formulation),(backend=:unknown,)))
    @test description([missing];roles=[:reference])==["Reference · method unavailable"]
    unknown_inner=IO.deserialize_value(Val(:formulation),
        (kind=:monte_carlo,inner=(backend=:unknown,),options=(;)))
    @test only(description([unknown_inner];roles=[:reference]))=="Reference · Monte Carlo"
    @test any(occursin("method unavailable",description(scope,value;compact=false))
        for (scope,value) in pairs(unknown_inner...))
    # Read retained consumed IDs when present; never advertise a saved inactive
    # route as used, or replace an explicit unknown selection with today's default.
    declared=NamedTuple(normal)
    consumed=merge(declared,(effective=merge(map(_ -> :default,declared.requested),
        (earth_impedance=:saad1996,pipe_impedance=nothing)),))
    decoded=IO.deserialize_value(Val(:formulation),consumed)
    @test occursin("earth Z=Saad",only(description([decoded];quantity=R)))
    @test any(last(scope)==(:pipe_impedance,) && value===nothing for (scope,value) in pairs(decoded...))

    # A programming error is not absent metadata. The real composed consumer
    # must execute the owned description and propagate its failure.
    struct BrokenEarthLeaf <: E.EarthImpedanceFormulation end
    formula_id(::BrokenEarthLeaf)=:BrokenTestLeaf
    description(::BrokenEarthLeaf;compact::Bool=false)=throw(ArgumentError("test-owned description failed"))
    formulation_options(::BrokenEarthLeaf)=FormulationOptions()
    Base.NamedTuple(value::BrokenEarthLeaf)=(identifier=formula_id(value),parameters=(;),options=(;))
    broken=E.LineParametersFormulation(merge(normal.methods,(earth_impedance=BrokenEarthLeaf(),)),
        normal.options,merge(normal.definitions,(earth_impedance=BrokenEarthLeaf(),)))
    bad=ParametricResult(nothing,[points[1]],(problems=[:one],formulations=[broken]), ComputationDetails((;)))
    @test_throws r"test-owned description failed" report(BenchmarkTableDefinition(),(reference=ref,candidate=bad))

    # Compare retained values and coordinates, not merely dimensions or file bytes.
    for row in result.published.comparisons
        selected=filter(term -> term.analysis==findfirst(==(row),result.published.comparisons),eachrow(result.table.terms))
        absolute=observe(row.error,E.absolute_error)
        relative=observe(row.error,E.relative_error)
        for term in selected
            @test term.response==ports[term.row] && term.excitation==ports[term.column]
            @test isequal(term.absolute_rms,absolute[term.row,term.column])
            @test isequal(term.relative_rms_percent,100relative[term.row,term.column])
        end
    end
    @test Z(ref)==z && Y(ref)==y && frequencies(ref)==f
    for frame in (result.table.terms,result.table.maxima,result.table.summary),column in eachcol(frame),value in column
        @test value isa Union{Number,Symbol,AbstractString,Missing}
    end
    for row in eachrow(result.table.maxima)
        evidence=filter(term -> term.analysis==row.analysis,eachrow(result.table.terms))
        largest=only(filter(term -> term.response==row.absolute_term_response &&
            term.excitation==row.absolute_term_excitation,evidence))
        @test row.maximum_absolute_rms==largest.absolute_rms
    end
    rendered=sprint(show,MIME"text/plain"(),result)
    # New leaf descriptions must also reach the composed human-readable report.
    @test all(leaf->occursin(description(leaf;compact=true),rendered),leaves)
end

@testitem "Descriptions / quantity identities ignore unrelated slots and retain composite controls" tags=[:unit] setup=[FormulaContractModels] begin
    IO=LineCableModels.ImportExport
    a=Formulation(earth_impedance=:saad1996)
    b=Formulation(earth_impedance=:xue2018)
    @test formula_id(a,Y)==formula_id(b,Y)
    @test formula_id(a,Z)!=formula_id(b,Z)
    @test formula_id(a,nothing)!=formula_id(b,nothing)
    for source in (a,b,MonteCarlo(a),LinearError(a),LineCableModelsFEM(),Formulation(:pscad))
        saved=IO.deserialize_value(Val(:formulation),NamedTuple(source))
        for quantity in (nothing,R,B)
            @test formula_id(source,quantity)==formula_id(saved,quantity)
        end
    end
    routed=Formulation(earth_impedance=(air=:default,earth=:default,mixed=:xue2018))
    other=Formulation(earth_impedance=(air=:default,earth=:default,mixed=:lucca1994))
    @test formula_id(routed,R)!=formula_id(other,R)
    @test formula_id(routed,B)==formula_id(other,B)
    selected_transfer=FormulaContractModels.SurfaceLaw(kinds=(:transfer,))
    internal=Formulation(internal_impedance=(inner=:default,outer=:default,transfer=selected_transfer))
    @test formula_id(internal,R)!=formula_id(Formulation(),R)
    @test formula_id(internal,Y)==formula_id(Formulation(),Y)
    overridden=Formulation(earth_admittance=formula(:default;parameters=(reference=:interface,)))
    @test formula_id(overridden,Y)!=formula_id(Formulation(),Y)
    @test formula_id(overridden,Z)==formula_id(Formulation(),Z)
    @test formula_id(MonteCarlo(a),Y)!=formula_id(LinearError(a),Y)
    @test ismissing(formula_id(missing,Y))
    incomplete=IO.deserialize_value(Val(:formulation),
        (backend=:coaxial,requested=(earth_admittance=(identifier=:default,),)))
    @test ismissing(formula_id(incomplete,Y))
    @test description([Formulation()];quantity=Y)==[
        "shunt geometry=coaxial; insulation Y=Lossless; semicon Y=Lossless; earth Y=Unified; soil law=Constant; temperature law=Linear"]
end
