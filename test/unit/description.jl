@testitem "Descriptions / partial recipes retain only declared ordered children" tags=[:unit, :engine] begin
    E = LineCableModels.Engine
    for selected in (
            Formulation(earth_impedance=(earth=:unified,)),
            Formulation(earth_admittance=(mixed=:unified, air=:unified)),
            Formulation(earth_impedance=(;)))
        record = NamedTuple(selected)
        for slot in (:earth_impedance, :earth_admittance)
            children = getproperty(selected.formulas, slot)
            children isa NamedTuple || continue
            routes = [last(scope)[2] for (scope, _) in pairs(selected)
                if length(last(scope)) == 2 && first(last(scope)) == slot]
            @test Tuple(routes) == keys(children)
            @test keys(getproperty(record.formulas, slot)) == keys(children)
        end
    end
    record = NamedTuple(Formulation(earth_impedance=(earth=:unified,)))
    malformed = merge(record, (formulas=merge(record.formulas,
        (earth_impedance=(ocean=record.formulas.earth_impedance.earth,),)),))
    @test_throws ArgumentError pairs(LineParametersFormulation, malformed)
end

@testitem "Descriptions / owner dispatch survives composed reports and saved declarations" tags=[:unit, :pscad, :slow] setup=[FormulaFixtures] begin
    using DataFrames, Statistics, PythonCall
    const PSCAD = Base.get_extension(LineCableModels, :LineCableModelsPSCADExt)
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    import LineCableModels: description, formula_id, formulation_options
    E=LineCableModels.Engine

    # Equal text must not become an identity, and a new leaf must appear in the real
    # report without being added to a reader or plotting lookup table.
    struct TestEarthLeaf{ID} <: E.EarthImpedanceFormulation end
    formula_id(::TestEarthLeaf{ID}) where {ID} = ID
    description(::Type{<:TestEarthLeaf{ID}};compact::Bool=false) where {ID} =
        compact ? "Display-"*string(ID) : "Test-owned explanation shared by distinct equations"
    description(value::TestEarthLeaf;compact::Bool=false) = description(typeof(value);compact)
    formulation_options(::TestEarthLeaf) = FormulationOptions()
    Base.NamedTuple(value::TestEarthLeaf) = (identifier=formula_id(value),parameters=(;),options=(;))
    normal=Formulation()
    leaves=(TestEarthLeaf{:TestAlpha}(),TestEarthLeaf{:TestBeta}())
    choices=[E.LineParametersFormulation(merge(normal.formulas,(earth_impedance=leaf,)),
        normal.options) for leaf in leaves]
    push!(choices,normal)
    @test description(LineCableModelsCoaxial())==description(LineCableModelsCoaxial)=="coaxial"
    @test description(Formulation(:pscad))==description(PSCAD.PSCADFormulation)=="PSCAD"
    pscad=Formulation(:pscad;options=(base_frequency=60.0,))
    pscad_capture=E.completed_formulation(pscad)
    @test first(pscad_capture.formulation_fields.Z).value==description(PSCAD.PSCADFormulation;compact=true)
    @test any(control -> control.text==description(PSCAD.PSCADFormulation,Val(:base_frequency),60.0;compact=true),
        first(pscad_capture.formulation_fields.Z).control_fields)
    @test description(PSCAD.PSCADFormulation,Val(:reduce_bundle),false;compact=true)==
        description(LineParametersFormulation,Val(:reduce_bundle),false;compact=true)
    for value in (normal,Formulation(:pscad),MonteCarlo(normal),LinearError(normal))
        @test description(typeof(value))==description(value)
    end
    labels=description([MonteCarlo(normal),LinearError(normal)];roles=[:reference,:result])
    @test labels == ["Reference · Monte Carlo; internal Z=Schelkunoff (options=(inner = NamedTuple(), outer = NamedTuple(), transfer = NamedTuple())); earth Z=Unified (options=(Γ = 0,)); earth Y=Unified (options=(Γ = 0,))",
        "LEP; internal Z=Schelkunoff (options=(inner = NamedTuple(), outer = NamedTuple(), transfer = NamedTuple())); earth Z=Unified (options=(Γ = 0,)); earth Y=Unified (options=(Γ = 0,))"]
    @test all(label->occursin("earth Z=Unified (options=(Γ = 0,))",label),labels)
    explicit_default=Formulation(earth_impedance=:unified)
    @test formula_id(normal,R)==formula_id(explicit_default,R)
    @test description([normal,explicit_default];quantity=R)==
        fill("internal Z=Schelkunoff (options=(inner = NamedTuple(), outer = NamedTuple(), transfer = NamedTuple())); earth Z=Unified (options=(Γ = 0,))",2)
    @test description([normal];quantity=R)==description([explicit_default];quantity=R)

    f=[0.1,50.,100.,1e3,1e6,1e7]
    z=reshape(complex.(collect(1.:24.),collect(101.:124.)),2,2,6)
    y=reshape(complex.(collect(201.:224.),collect(301.:324.)),2,2,6)*1e-6
    ports=["a","b"]
    ref=LineParameters(z,y,f;terminals=ports,line_length=1.0,
        details=ComputationDetails(;formulation=Formulation(),choices=()))
    source_id=LineCableModels.Commons.gridpoint_id().source_id
    completed(value,selection,index)=E.retain_gridpoint(value,
        LineCableModels.Commons.gridpoint_id(;source_id,formulation_index=index);
        fields=(formulation=selection,choices=()))
    points=[completed(LineParameters(scale*z,y,f; line_length=1.0,terminals=ports),choices[index],index)
        for (index,scale) in enumerate((1.1,1.2,1.3))]
    results=ParametricResult(points,
        fieldtype(ParametricResult, :failures)(),
        CartesianIndices((1,length(choices))), ComputationDetails((;)))
    source=(reference=ref,result=results)
    result=report(BenchmarkTableDefinition((R,X,G,B);bands=(:all,:dc,:harmonic,:narrow,:wide)),source)
    @test any(contains("TestAlpha"),result.tables.formulations.label)
    @test any(contains("TestBeta"),result.tables.formulations.label)
    observed_labels=LineCableModels.Commons.observation_labels(result.observed;request=R)
    @test occursin(description(leaves[1];compact=true),observed_labels[1])
    @test all(label -> !occursin("internal Z",label),observed_labels)
    @test length(observed_labels)==length(result.observed)
    @test length(result.observed)==3
    @test result.observed[1].gridpoint.formulations==NamedTuple(choices[1])
    labels=description(choices[[3,1]];quantity=R)
    @test occursin("earth Z=Unified",labels[1]) && occursin("earth Z=Display-TestAlpha",labels[2])
    @test all(label->occursin("internal Z=Schelkunoff (options=(inner = NamedTuple(), outer = NamedTuple(), transfer = NamedTuple()))",label),labels)
    @test all(feature -> !occursin("earth Y",join(feature.relative.formula)),
        filter(feature -> feature.quantity in (:R,:X),result.tables.features))

    for order in (:before,:after)
        native=Formulation(earth_impedance=formula(:carson1926;
            equivalent_earth=formula(:default;order)))
        label=only(description([native];quantity=R))
        @test occursin(string(order)*" FrequencyDependent",label)
        @test !occursin("FormulaDefinition{",label)
    end
    routed=Formulation(
        earth_impedance=(air=:carson1926,earth=:pollaczek1926,mixed=:lucca1994),
        earth_admittance=formula(:default;options=(integration=(method=:quad,options=(rtol=1e-9,)),)))
    routed_result=ParametricResult([completed(points[1],routed,1),completed(points[2],normal,2)],
        fieldtype(ParametricResult, :failures)(),
        CartesianIndices((1,2)), ComputationDetails((;)))
    routed_report=report(BenchmarkTableDefinition((R,B);bands=(:all,)),
        (reference=ref,result=routed_result))
    @test any(contains("earth Z(air)=Carson"),first(routed_report.tables.features).relative.formula)
    # Inspection retains a native selection without executing its equation.
    selected_inner=FormulaFixtures.SurfaceLaw()
    single=Formulation(insulation_admittance=FormulaFixtures.InsulationLaw(),
        internal_impedance=selected_inner)
    single_data=ParametricResult([completed(points[1],single,1)],
        fieldtype(ParametricResult, :failures)(),
        CartesianIndices((1,1)), ComputationDetails((;)))
    single_report=report(BenchmarkTableDefinition((R,B);bands=(:all,)),
        (reference=ref,result=single_data))
    single_label=only(last(single_report.tables.features).relative.formula)
    @test occursin(description(single.formulas.insulation_admittance;compact=true),single_label)
    @test !occursin("scale",single_label) # a constant control is not a legend difference
    insulation=only(filter(field -> field.meaning==(:insulation_admittance,),only(single_report.observed).gridpoint.formulation_fields.Y))
    @test any(control -> last(control.scope)===:scale,insulation.control_fields)
    @test !occursin("insulation Y",only(first(single_report.tables.features).relative.formula))
    @test isempty(selected_inner.evaluations) && isempty(selected_inner.state_inputs)
    pscad_labels=description([Formulation(:pscad),Formulation(:pscad;options=(base_frequency=60.0,))];
        roles=[:reference,:reference])
    @test occursin("base_frequency = 50.0",first(pscad_labels)) &&
        occursin("base_frequency = 60.0",last(pscad_labels))
    @test description([missing];roles=[:reference])==["Reference · method unavailable"]

    # A programming error is not absent metadata. The real composed consumer
    # must execute the owned description and propagate its failure.
    struct BrokenEarthLeaf <: E.EarthImpedanceFormulation end
    formula_id(::BrokenEarthLeaf)=:BrokenTestLeaf
    description(::BrokenEarthLeaf;compact::Bool=false)=throw(ArgumentError("test-owned description failed"))
    formulation_options(::BrokenEarthLeaf)=FormulationOptions()
    Base.NamedTuple(value::BrokenEarthLeaf)=(identifier=formula_id(value),parameters=(;),options=(;))
    broken=E.LineParametersFormulation(merge(normal.formulas,(earth_impedance=BrokenEarthLeaf(),)),
        normal.options)
    bad=ParametricResult([points[1]],
        fieldtype(ParametricResult, :failures)(),
        CartesianIndices((1,1)), ComputationDetails((;)))
    @test_throws r"test-owned description failed" E.completed_formulation(broken)
    @test only(report(BenchmarkTableDefinition(),(reference=ref,result=bad)).observed).gridpoint.formulations==NamedTuple(choices[1])

    # Compare retained values and coordinates, not merely dimensions or file bytes.
    for point in result.observed,row in point.errors
        selected=filter(term -> term.result_id==row.result_id && term.request==row.request &&
            term.band==row.band && term.normalization==row.normalization,eachrow(result.tables.terms))
        for term in selected
            @test term.response==ports[term.row] && term.excitation==ports[term.column]
            @test isequal(term.absolute_rms,row.absolute[term.row,term.column])
            @test isequal(term.relative_rms_percent,100row.relative[term.row,term.column])
        end
    end
    @test Z(ref)==z && Y(ref)==y && frequencies(ref)==f
    for row in eachrow(result.tables.maxima)
        evidence=filter(term -> term.result_id==row.result_id && term.request==row.request && term.band==row.band && term.normalization==row.normalization,eachrow(result.tables.terms))
        largest=only(filter(term -> term.response==row.absolute_term_response &&
            term.excitation==row.absolute_term_excitation,evidence))
        @test row.maximum_absolute_rms==largest.absolute_rms
    end
    rendered=sprint(show,MIME"text/plain"(),result)
    # New leaf descriptions must also reach the composed human-readable report.
    @test all(leaf->occursin(description(leaf;compact=true),rendered),leaves)
end

@testitem "Descriptions / quantity identities ignore unrelated slots and retain composite controls" tags=[:unit, :uq] setup=[FormulaFixtures] begin
    a=Formulation(earth_impedance=:saad1996)
    b=Formulation(earth_impedance=:xue2018)
    @test formula_id(a,Y)==formula_id(b,Y)
    @test formula_id(a,Z)!=formula_id(b,Z)
    @test formula_id(a,nothing)!=formula_id(b,nothing)
    routed=Formulation(earth_impedance=(air=:default,earth=:default,mixed=:xue2018))
    other=Formulation(earth_impedance=(air=:default,earth=:default,mixed=:lucca1994))
    @test formula_id(routed,R)!=formula_id(other,R)
    @test formula_id(routed,B)==formula_id(other,B)
    internal=Formulation(internal_impedance=FormulaFixtures.SurfaceLaw())
    @test formula_id(internal,R)!=formula_id(Formulation(),R)
    @test formula_id(internal,Y)==formula_id(Formulation(),Y)
    overridden=Formulation(earth_admittance=formula(:default;
        options=(integration=(method=:quad,options=(rtol=1e-9,)),)))
    @test formula_id(overridden,Y)!=formula_id(Formulation(),Y)
    @test formula_id(overridden,Z)==formula_id(Formulation(),Z)
    @test formula_id(MonteCarlo(a),Y)!=formula_id(LinearError(a),Y)
    @test ismissing(formula_id(missing,Y))
    @test description([Formulation()];quantity=Y)==[
        "shunt geometry=equivalent annular layer; insulation Y=Lossless; semicon Y=Lossless; earth Y=Unified (options=(Γ = 0,)); soil law=Constant; temperature law=Linear"]
end

@testitem "Descriptions / detached differences use compact owner dispatch" tags=[:unit, :pscad, :slow] begin
    using PythonCall
    const PSCAD = Base.get_extension(LineCableModels, :LineCableModelsPSCADExt)
    import LineCableModels: description,formula_id,formulation_options
    using LineCableModels.Engine: retain_gridpoint
    using LineCableModels.Commons: gridpoint_id,observation_labels,observation_groups
    E=LineCableModels.Engine
    IO=LineCableModels.ImportExport
    raw=LineParameters(fill(1.0+2im,1,1,3),fill(3e-6+4e-6im,1,1,3),[1.,10.,100.]; line_length=1.0)
    source_id=gridpoint_id().source_id
    function point(selection,index;rho=100.,problem=1)
        fields=(formulation=selection,choices=(),inputs=(radius=.01,rho=rho,
            field_descriptions=(radius=(name="radius",unit="m"),rho=(name="electrical resistivity",unit="Ω·m"))))
        ObservedResult(retain_gridpoint(raw,gridpoint_id(;source_id,problem_index=problem,formulation_index=index);fields))
    end
    a=Formulation(earth_impedance=:unified,earth_admittance=:unified)
    constants=CableConstantsFormulation()
    points=[point(a,1),point(constants,2)]
    for request in (R,X,G,B)
        method=request in (R,X) ? a.formulas.earth_impedance : a.formulas.earth_admittance
        @test observation_labels(points;request)==[description(method;compact=true),description(constants;compact=true)]
    end
    @test all(field -> field.meaning ∉ ((:earth_impedance,),(:earth_admittance,)),points[2].gridpoint.formulation_fields.all)
    physical=[point(a,1),point(a,1;rho=500.,problem=2)]
    @test observation_labels(physical;request=R)==["electrical resistivity=100.0 Ω·m","electrical resistivity=500.0 Ω·m"]
    @test length(observation_groups(physical;request=R))==2
    other=Formulation(earth_impedance=:pollaczek1926)
    mixed=[physical[1],point(other,2;rho=500.,problem=2)]
    labels=observation_labels(mixed;request=R)
    @test occursin(description(other.formulas.earth_impedance;compact=true),labels[2])
    @test occursin("500.0 Ω·m",labels[2])
    reductions=[point(Formulation(options=(kron_reduction=choice,)),i) for (i,choice) in enumerate((false,true))]
    @test observation_labels(reductions;request=R)==[
        description(LineParametersFormulation,Val(:kron_reduction),v;compact=true) for v in (false,true)]
    frequencies=[Formulation(:pscad;options=(base_frequency=choice,)) for choice in (50.0,60.0)]
    @test observation_labels([point(f,i) for (i,f) in enumerate(frequencies)];request=B)==[
        description(PSCAD.PSCADFormulation,Val(:base_frequency),f.options.data.base_frequency;compact=true)
        for f in frequencies]
    reordered=Formulation(options=(ideal_transposition=false,kron_reduction=true,reduce_bundle=true))
    @test observation_labels([point(a,1),point(reordered,2)];request=R)==fill(description(a.formulas.earth_impedance;compact=true),2)

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
    @test occursin(description(typeof(gamma.formulas.earth_impedance),Val(:Γ),1.0;compact=true),gamma_labels[2])
    gamma_samples=[0.1,0.2,0.3]
    gamma_vector=Formulation(earth_impedance=formula(:unified;options=(Γ=gamma_samples,)),
        earth_admittance=formula(:unified;options=(Γ=gamma_samples,)))
    for request in (R,B)
        gamma_vector_labels=observation_labels([point(a,1),point(gamma_vector,2)];request)
        selected=request===R ? gamma_vector.formulas.earth_impedance : gamma_vector.formulas.earth_admittance
        @test occursin(description(typeof(selected),Val(:Γ),gamma_samples;compact=true),last(gamma_vector_labels))
    end
    @test description(E.EarthImpedance.Formula{:default};compact=true)==description(a.formulas.earth_impedance;compact=true)
    @test description(E.EarthAdmittance.Formula{:default};compact=true)==description(a.formulas.earth_admittance;compact=true)

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
    f=E.LineParametersFormulation(merge(a.formulas,(earth_impedance=DetachedNamingLeaf(),)),a.options)
    extended=point(f,3)
    @test !isempty(compact_calls) && all(compact_calls)
    poison[]=true
    Z(raw) .= NaN
    expected=[description(a.formulas.earth_impedance;compact=true),"Independent field solution"]
    @test observation_labels([points[1],extended];request=R)==expected
    restored=IO.deserialize_value(IO.serialize_value(extended))
    @test observation_labels([points[1],restored];request=R)==expected
    @test restored.gridpoint==extended.gridpoint
end
