@testitem "Descriptions / partial recipes retain only declared ordered children" tags=[:unit] begin
    E = LineCableModels.Engine
    IE = LineCableModels.ImportExport
    for selected in (
            Formulation(earth_impedance=(earth=:unified,)),
            Formulation(earth_admittance=(mixed=:unified, air=:unified)),
            Formulation(internal_impedance=(outer=:default,)),
            Formulation(earth_impedance=(;), internal_impedance=(outer=nothing,)))
        record = NamedTuple(selected)
        saved = IE.deserialize_value(Val(:formulation), record)
        for quantity in (nothing, Z, Y)
            @test formula_id(selected, quantity) == formula_id(saved, quantity)
            @test description(selected, quantity) == description(saved, quantity)
        end
        for slot in (:earth_impedance, :earth_admittance, :internal_impedance)
            children = getproperty(selected.methods, slot)
            children isa NamedTuple || continue
            routes = [last(scope)[2] for (scope, _) in pairs(selected)
                if length(last(scope)) == 2 && first(last(scope)) == slot]
            @test Tuple(routes) == keys(children)
            @test keys(getproperty(record.methods, slot)) == keys(children)
        end
    end
    record = NamedTuple(Formulation(earth_impedance=(earth=:unified,)))
    malformed = merge(record, (methods=merge(record.methods,
        (earth_impedance=(ocean=record.methods.earth_impedance.earth,),)),))
    @test_throws ArgumentError pairs(LineParametersFormulation, malformed)
end

@testitem "Descriptions / owner dispatch survives composed reports and saved declarations" tags=[:unit] setup=[FormulaFixtures] begin
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
    pscad=Formulation(:pscad;options=(base_frequency=60.0,))
    pscad_capture=E.completed_formulation(pscad)
    pscad_saved=IO.deserialize_value(Val(:formulation),NamedTuple(pscad))
    @test E.completed_formulation(pscad_saved,NamedTuple(pscad)).formulation_fields==pscad_capture.formulation_fields
    @test first(pscad_capture.formulation_fields.Z).value==description(PSCAD.PSCADFormulation;compact=true)
    @test any(control -> control.text==description(PSCAD.PSCADFormulation,Val(:base_frequency),60.0;compact=true),
        first(pscad_capture.formulation_fields.Z).control_fields)
    @test description(PSCAD.PSCADFormulation,Val(:reduce_bundle),false;compact=true)==
        description(LineParametersFormulation,Val(:reduce_bundle),false;compact=true)
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
    ref=LineParameters(z,y,f;details=ComputationDetails(;coordinates=ports,E.completed_formulation(LineCableModelsFEM())...,))
    source_id=LineCableModels.Grammar.gridpoint_id().source_id
    completed(value,selection,index)=E.retain_gridpoint(value,
        LineCableModels.Grammar.gridpoint_id(;source_id,formulation_index=index);
        fields=E.completed_formulation(selection))
    points=[completed(LineParameters(scale*z,y,f;details=ComputationDetails(;coordinates=ports,)),choices[index],index)
        for (index,scale) in enumerate((1.1,1.2,1.3))]
    candidates=ParametricResult(nothing,points,(problems=[:one],formulations=choices), ComputationDetails((;)))
    source=(reference=ref,candidate=candidates)
    result=report(BenchmarkTableDefinition((R,X,G,B);bands=(:all,:dc,:harmonic,:narrow,:wide)),source)
    @test any(contains("TestAlpha"),result.tables.formulations.label)
    @test any(contains("TestBeta"),result.tables.formulations.label)
    observed_labels=LineCableModels.Grammar.observation_labels(result.observed;request=R)
    @test occursin(description(leaves[1];compact=true),observed_labels[1])
    @test all(label -> !occursin("internal Z",label),observed_labels)
    @test length(observed_labels)==length(result.observed)
    @test length(result.observed)==3
    @test result.observed[1].gridpoint.formulations==NamedTuple(choices[1])
    @test all(leaf -> occursin("Test-owned explanation",description(leaf;compact=false)),leaves)
    labels=description(choices[[3,1]];quantity=R)
    @test occursin("earth Z=Unified",labels[1]) && occursin("earth Z=Display-TestAlpha",labels[2])
    @test all(label->!occursin("internal Z",label),labels)
    @test all(feature -> !occursin("earth Y",join(feature.relative.formula)),
        filter(feature -> feature.quantity in (:R,:X),result.tables.features))

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
    selected_inner=FormulaFixtures.SurfaceLaw(kinds=(:inner,))
    routed=Formulation(
        internal_impedance=(inner=selected_inner,outer=:default,transfer=:default),
        earth_impedance=(air=:carson1926,earth=:pollaczek1926,mixed=:lucca1994),
        earth_admittance=formula(:default;options=(integration=(method=:quad,options=(rtol=1e-9,)),)))
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
    routed_result=ParametricResult(nothing,[completed(points[1],routed,1),completed(points[2],normal,2)],
        (problems=[:one],formulations=[routed,normal]), ComputationDetails((;)))
    routed_report=report(BenchmarkTableDefinition((R,B);bands=(:all,)),
        (reference=ref,candidate=routed_result))
    @test any(contains("earth Z(air)=Carson"),first(routed_report.tables.features).relative.formula)
    @test haskey(first(routed_report.observed).gridpoint.formulations.methods.internal_impedance,:inner)
    @test isempty(selected_inner.evaluations) && isempty(selected_inner.preparations)
    single=Formulation(insulation_admittance=FormulaFixtures.InsulationLaw())
    single_data=ParametricResult(nothing,[completed(points[1],single,1)],
        (problems=[:one],formulations=[single]), ComputationDetails((;)))
    single_report=report(BenchmarkTableDefinition((R,B);bands=(:all,)),
        (reference=ref,candidate=single_data))
    single_label=only(last(single_report.tables.features).relative.formula)
    @test occursin(description(single.methods.insulation_admittance;compact=true),single_label)
    @test !occursin("scale",single_label) # a constant control is not a legend difference
    insulation=only(filter(field -> field.meaning==(:insulation_admittance,),only(single_report.observed).gridpoint.formulation_fields.Y))
    @test any(control -> last(control.scope)===:scale,insulation.control_fields)
    @test !occursin("insulation Y",only(first(single_report.tables.features).relative.formula))
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
    consumed=merge(declared,(methods=merge(declared.methods,
        (earth_impedance=NamedTuple(E.EarthImpedance.Formula(:saad1996)),
            pipe_impedance=nothing)),))
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
    @test_throws r"test-owned description failed" E.completed_formulation(broken)
    @test only(report(BenchmarkTableDefinition(),(reference=ref,candidate=bad)).observed).gridpoint.formulations==NamedTuple(choices[1])

    # Compare retained values and coordinates, not merely dimensions or file bytes.
    for point in result.observed,row in point.errors
        selected=filter(term -> term.candidate_id==row.candidate_id && term.request==row.request &&
            term.band==row.band && term.normalization==row.normalization,eachrow(result.tables.terms))
        for term in selected
            @test term.response==ports[term.row] && term.excitation==ports[term.column]
            @test isequal(term.absolute_rms,row.absolute[term.row,term.column])
            @test isequal(term.relative_rms_percent,100row.relative[term.row,term.column])
        end
    end
    @test Z(ref)==z && Y(ref)==y && frequencies(ref)==f
    for row in eachrow(result.tables.maxima)
        evidence=filter(term -> term.candidate_id==row.candidate_id && term.request==row.request && term.band==row.band && term.normalization==row.normalization,eachrow(result.tables.terms))
        largest=only(filter(term -> term.response==row.absolute_term_response &&
            term.excitation==row.absolute_term_excitation,evidence))
        @test row.maximum_absolute_rms==largest.absolute_rms
    end
    rendered=sprint(show,MIME"text/plain"(),result)
    # New leaf descriptions must also reach the composed human-readable report.
    @test all(leaf->occursin(description(leaf;compact=true),rendered),leaves)
end

@testitem "Descriptions / quantity identities ignore unrelated slots and retain composite controls" tags=[:unit] setup=[FormulaFixtures] begin
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
    selected_transfer=FormulaFixtures.SurfaceLaw(kinds=(:transfer,))
    internal=Formulation(internal_impedance=(inner=:default,outer=:default,transfer=selected_transfer))
    @test formula_id(internal,R)!=formula_id(Formulation(),R)
    @test formula_id(internal,Y)==formula_id(Formulation(),Y)
    overridden=Formulation(earth_admittance=formula(:default;
        options=(integration=(method=:quad,options=(rtol=1e-9,)),)))
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

@testitem "Descriptions / detached differences use compact owner dispatch" tags=[:unit] begin
    import LineCableModels: description,formula_id,formulation_options
    using LineCableModels.Engine: completed_formulation,retain_gridpoint
    using LineCableModels.Grammar: gridpoint_id,observation_labels,observation_groups
    E=LineCableModels.Engine
    IO=LineCableModels.ImportExport
    raw=LineParameters(fill(1.0+2im,1,1,3),fill(3e-6+4e-6im,1,1,3),[1.,10.,100.])
    source_id=gridpoint_id().source_id
    function point(selection,index;rho=100.,problem=1)
        fields=merge(completed_formulation(selection),(inputs=(radius=.01,rho=rho,
            field_descriptions=(radius=(name="radius",unit="m"),rho=(name="electrical resistivity",unit="Ω·m"))),))
        ObservedResult(retain_gridpoint(raw,gridpoint_id(;source_id,problem_index=problem,formulation_index=index);fields))
    end
    a=Formulation(earth_impedance=:unified,earth_admittance=:unified)
    fem=Formulation(:LineCableModelsFEM)
    points=[point(a,1),point(fem,2)]
    for request in (R,X,G,B)
        method=request in (R,X) ? a.methods.earth_impedance : a.methods.earth_admittance
        @test observation_labels(points;request)==[description(method;compact=true),description(fem;compact=true)]
    end
    @test all(field -> field.meaning ∉ ((:earth_impedance,),(:earth_admittance,)),points[2].gridpoint.formulation_fields.all)
    physical=[point(a,1),point(a,1;rho=500.,problem=2)]
    @test observation_labels(physical;request=R)==["electrical resistivity=100.0 Ω·m","electrical resistivity=500.0 Ω·m"]
    @test length(observation_groups(physical;request=R))==2
    other=Formulation(earth_impedance=:pollaczek1926)
    mixed=[physical[1],point(other,2;rho=500.,problem=2)]
    labels=observation_labels(mixed;request=R)
    @test occursin(description(other.methods.earth_impedance;compact=true),labels[2])
    @test occursin("500.0 Ω·m",labels[2])
    reductions=[point(Formulation(options=(kron_reduction=choice,)),i) for (i,choice) in enumerate((false,true))]
    @test observation_labels(reductions;request=R)==[
        description(LineParametersFormulation,Val(:kron_reduction),v;compact=true) for v in (false,true)]
    physics=[Formulation(:LineCableModelsFEM;options=(physics=choice,)) for choice in (:quasi_tem,:quasi_fw)]
    @test observation_labels([point(f,i) for (i,f) in enumerate(physics)];request=B)==[
        description(LineCableModelsFEM,Val(:physics),f.options.data.physics;compact=true) for f in physics]
    reordered=Formulation(options=(ideal_transposition=true,kron_reduction=true,reduce_bundle=true))
    @test observation_labels([point(a,1),point(reordered,2)];request=R)==fill(description(a.methods.earth_impedance;compact=true),2)

    layer(rho)=(rho=rho,field_descriptions=(rho=(name="electrical resistivity",unit="Ω·m"),))
    layered=[ObservedResult(merge(p.gridpoint,(inputs=(layers=[layer(100.),layer(rho)],),)),p.quantities,p.errors,p.timings)
        for (p,rho) in zip(physical,(200.,300.))]
    @test observation_labels(layered;request=R)==["electrical resistivity[2]=200.0 Ω·m","electrical resistivity[2]=300.0 Ω·m"]

    nested_radius=[ObservedResult(merge(p.gridpoint,(inputs=(designs=[(regions=[(
            r=radius,field_descriptions=(r=(name="radius",unit="m"),))],)],),)),
            p.quantities,p.errors,p.timings)
        for (p,radius) in zip(physical,(0.005,0.01))]
    @test observation_labels(nested_radius;request=R)==
        ["radius[1][1]=0.005 m","radius[1][1]=0.01 m"]

    gamma=Formulation(earth_impedance=formula(:unified;options=(Γ=1.0,)))
    gamma_labels=observation_labels([point(a,1),point(gamma,2)];request=R)
    @test all(occursin("Γ=",label) for label in gamma_labels)
    @test occursin(description(typeof(gamma.methods.earth_impedance),Val(:Γ),1.0;compact=true),gamma_labels[2])
    gamma_samples=[0.1,0.2,0.3]
    gamma_vector=Formulation(earth_impedance=formula(:unified;options=(Γ=gamma_samples,)),
        earth_admittance=formula(:unified;options=(Γ=gamma_samples,)))
    for request in (R,B)
        gamma_vector_labels=observation_labels([point(a,1),point(gamma_vector,2)];request)
        selected=request===R ? gamma_vector.methods.earth_impedance : gamma_vector.methods.earth_admittance
        @test occursin(description(typeof(selected),Val(:Γ),gamma_samples;compact=true),last(gamma_vector_labels))
    end
    @test description(E.EarthImpedance.Formula{:default};compact=true)==description(a.methods.earth_impedance;compact=true)
    @test description(E.EarthAdmittance.Formula{:default};compact=true)==description(a.methods.earth_admittance;compact=true)

    # New scientific owner, deliberately unrelated identifier and compact text.
    const compact_calls=Bool[]
    const poison=Ref(false)
    struct DetachedNamingLeaf <: E.EarthImpedanceFormulation end
    formula_id(::DetachedNamingLeaf)=:fixture_unrelated_identifier
    function description(::DetachedNamingLeaf;compact::Bool=false)
        poison[] && error("live description reopened")
        push!(compact_calls,compact)
        compact ? "Independent field solution" : "Verbose explanation for this fixture"
    end
    formulation_options(::DetachedNamingLeaf)=FormulationOptions()
    Base.NamedTuple(::DetachedNamingLeaf)=(identifier=:fixture_unrelated_identifier,parameters=(;),options=(;))
    f=E.LineParametersFormulation(merge(a.methods,(earth_impedance=DetachedNamingLeaf(),)),a.options,
        merge(a.definitions,(earth_impedance=DetachedNamingLeaf(),)))
    extended=point(f,3)
    @test !isempty(compact_calls) && all(compact_calls)
    poison[]=true
    Z(raw) .= NaN
    expected=[description(a.methods.earth_impedance;compact=true),"Independent field solution"]
    @test observation_labels([points[1],extended];request=R)==expected
    restored=IO.deserialize_value(IO.serialize_value(extended))
    @test observation_labels([points[1],restored];request=R)==expected
    @test restored.gridpoint==extended.gridpoint
end
