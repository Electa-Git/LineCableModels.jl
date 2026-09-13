@testitem "Descriptions / owner dispatch survives composed reports and saved declarations" tags=[:unit] begin
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
    formulation_options(::TestEarthLeaf) = (;)
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
    @test description([MonteCarlo(normal),LinearError(normal)];
        roles=[:reference,:candidate],indices=[0,1])==["Reference · Monte Carlo","F1 · LEP"]

    f=[0.1,50.,100.,1e3,1e6,1e7]
    z=reshape(complex.(collect(1.:24.),collect(101.:124.)),2,2,6)
    y=reshape(complex.(collect(201.:224.),collect(301.:324.)),2,2,6)*1e-6
    ports=["a","b"]
    ref=LineParameters(z,y,f;details=(coordinates=ports,formulations=NamedTuple(LineCableModelsFEM()),))
    points=[LineParameters(scale*z,scale*y,f;details=(coordinates=ports,)) for scale in (1.1,1.2,1.3)]
    candidates=ParametricResult(nothing,points,(problems=[:one],formulations=choices),(;))
    source=(reference=ref,candidate=candidates)
    result=report(BenchmarkTableDefinition((R,X,G,B);bands=(:all,:dc,:harmonic,:narrow,:wide)),source)
    @test result.table.formulations.label[1]=="Reference · FEM"
    @test any(contains("TestAlpha"),result.table.formulations.label)
    @test any(contains("TestBeta"),result.table.formulations.label)
    @test count(contains("Test-owned explanation"),result.table.formula_details.selection)==2
    @test length(result.published.candidate.metadata.formulation_sources)==3
    @test result.published.candidate.metadata.formulation_sources[1]===choices[1]
    @test description(choices[[3,1]];indices=[3,1],quantity=R)==
        ["F3 · earth Z=default","F1 · earth Z=Display-TestAlpha"]
    @test all(feature -> !occursin("earth Y",join(feature.relative.formula)),
        filter(feature -> feature.quantity in (:R,:X),result.table.features))

    for native in (normal,Formulation(earth_impedance=:Saad1996),LineCableModelsFEM(),
            Formulation(:pscad),MonteCarlo(normal),LinearError(normal))
        saved=IO.deserialize_value(Val(:formulation),NamedTuple(native))
        @test description([saved];roles=[:reference])==description([native];roles=[:reference])
        @test [(scope,formula_id(value)) for (scope,value) in pairs(saved...)]==
            [(scope,formula_id(value)) for (scope,value) in pairs(native)]
    end
    saved_pair=[IO.deserialize_value(Val(:formulation),NamedTuple(method)) for method in
        (MonteCarlo(normal),LinearError(normal))]
    # Ordinary collection promotion must not erase a retained method's owner.
    @test description(saved_pair;roles=[:reference,:candidate],indices=[0,1])==
        ["Reference · Monte Carlo","F1 · LEP"]
    for order in (:before,:after)
        native=Formulation(earth_impedance=formula(:Carson1926;
            equivalent_earth=formula(:default;order)))
        saved=IO.deserialize_value(Val(:formulation),NamedTuple(native))
        label=only(description([native];quantity=R))
        @test occursin(string(order)*" FrequencyDependent",label)
        @test !occursin("FormulaDefinition{",label)
        @test description([saved];quantity=R)==[label]
    end
    # Describing a hook must never execute it. This same leaf/route composition
    # is exercised below through real reports, including an outer UQ wrapper.
    called=Ref(false)
    hook=(args...)->(called[]=true; error("description executed a scientific hook"))
    LineCableModels.computation_options(::LineCableModels.FormulaMethod{:default,
        typeof(LineCableModels.Engine.InternalImpedance.internal_impedance),Tuple{Val{:inner}}},
        ::typeof(hook))=(;)
    routed=Formulation(
        internal_impedance=formula(:default;hooks=(inner=hook,)),
        earth_impedance=(air=:Carson1926,earth=:Pollaczek1926,mixed=:Lucca1994),
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
        (problems=[:one],formulations=[routed,normal]),(;))
    routed_report=report(BenchmarkTableDefinition((R,B);bands=(:all,)),
        (reference=ref,candidate=routed_result))
    @test any(contains("earth Z(air)=Carson1926"),first(routed_report.table.features).relative.formula)
    @test any(contains("inner"),routed_report.table.formula_details.selection)
    @test !called[]
    LineCableModels.computation_options(::LineCableModels.FormulaMethod{:default,
        typeof(LineCableModels.Engine.InsulationAdmittance.insulation_material)},
        ::typeof(hook))=(;)
    single=Formulation(insulation_admittance=formula(:default;hooks=(contribution=hook,)))
    single_data=ParametricResult(nothing,points[1:1],
        (problems=[:one],formulations=[single]),(;))
    single_report=report(BenchmarkTableDefinition((R,B);bands=(:all,)),
        (reference=ref,candidate=single_data))
    single_label=only(last(single_report.table.features).relative.formula)
    @test occursin("insulation Y=default",single_label) && occursin("contribution",single_label)
    @test !occursin("insulation Y",only(first(single_report.table.features).relative.formula))
    @test description([single];quantity=B)==description([
        IO.deserialize_value(Val(:formulation),NamedTuple(single))];quantity=B)
    @test !called[]
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
        (earth_impedance=:Saad1996,pipe_impedance=nothing)),))
    decoded=IO.deserialize_value(Val(:formulation),consumed)
    @test occursin("earth Z=Saad1996",only(description([decoded];quantity=R)))
    @test any(last(scope)==(:pipe_impedance,) && value===nothing for (scope,value) in pairs(decoded...))

    # A programming error is not absent metadata. The real composed consumer
    # must execute the owned description and propagate its failure.
    struct BrokenEarthLeaf <: E.EarthImpedanceFormulation end
    formula_id(::BrokenEarthLeaf)=:BrokenTestLeaf
    description(::BrokenEarthLeaf;compact::Bool=false)=throw(ArgumentError("test-owned description failed"))
    formulation_options(::BrokenEarthLeaf)=(;)
    broken=E.LineParametersFormulation(merge(normal.methods,(earth_impedance=BrokenEarthLeaf(),)),
        normal.options,merge(normal.definitions,(earth_impedance=BrokenEarthLeaf(),)))
    bad=ParametricResult(nothing,[points[1]],(problems=[:one],formulations=[broken]),(;))
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
    @test occursin("Per-term RMS maxima",sprint(show,MIME"text/plain"(),result))
end
