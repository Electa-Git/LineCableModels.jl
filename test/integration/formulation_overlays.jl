@testitem "Makie addons / complete formulation overlays preserve labels, cells and exports" tags=[:visual] begin
    using CairoMakie
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    f=[1.,10.,100.]
    z=reshape(complex.(collect(1.:12.),collect(21.:32.)),2,2,3)
    y=1e-6im .* z
    options=(backend=:cairo,display_plot=false,controls=false,length_unit=:base,
        quantity_units=:base,freq_unit=:base,open_export=false)
    common=(backend=:coaxial,options=(reduce_bundle=false,),)
    selections=[merge(common,(requested=(earth_impedance=(air=NamedTuple(formula(:Carson1926)),
        earth=NamedTuple(formula(id)),mixed=NamedTuple(formula(:Lucca1994))),),)) for id in (:default,:Pollaczek1926)]
    reference=LineParameters(PhaseDomain,copy(z),copy(y),copy(f);details=(coordinates=["a","b"],))
    points=[LineParameters(PhaseDomain,factor*z,factor*y,copy(f);details=(coordinates=["a","b"],)) for factor in (2,3)]
    candidates=ParametricResult(nothing,points,(problems=[:one],formulations=selections),(;))
    baseline=(result=reference,metadata=(port_order=["a","b"],formulation=selections[1],axes=nothing))
    publication=report(BenchmarkTableDefinition(),(reference=baseline,candidate=candidates))
    series_attributes=((marker=:circle, markersize=8), (marker=nothing,),
        (linestyle=:dash, marker=nothing))
    plots=LineCableModels.plot(publication,(Z,);options...,series_attributes)
    @test length(plots)==2
    @test sum(length(page.axes) for page in plots)==8
    expected_labels=publication.table.formulations.label
    for (q,page) in zip((real,imag),plots)
        @test Set(values(page.addon_state.labels))==Set(expected_labels)
        for ((i,j),panel) in page.addon_state.panel_data
            curves=filter(item -> item isa Makie.Lines,panel.axis.scene.plots)
            @test length(curves)==3
            markers=only(filter(item -> item isa Makie.Scatter,panel.axis.scene.plots))
            @test markers[1][] == first(curves)[1][]
            @test last(curves).linestyle[] == Makie.to_linestyle(:dash)
            for (factor,curve) in zip((1,2,3),curves)
                @test first.(curve[1][]) ≈ f
                @test last.(curve[1][]) ≈ factor .* q.(z[i,j,:])
            end
        end
        @test length(page.addon_state.formulations.records)==length(selections)
    end
    filtered=LineCableModels.plot(publication,(R,);formulations=2,options...)
    for (all_axis,filtered_axis) in zip(first(plots).axes,filtered.axes)
        all_curves=filter(item -> item isa Makie.Lines,all_axis.scene.plots)
        selected=filter(item -> item isa Makie.Lines,filtered_axis.scene.plots)
        @test length(selected)==2
        @test selected[2].color[] == all_curves[3].color[]
    end
    @test Set(values(filtered.addon_state.labels))==Set(expected_labels[[1,3]])
    signed=LineCableModels.plot(publication;ydata=(G,),options...,controls=true)
    @test haskey(signed.controls,:ylog)
    signed.controls[:ylog].active[]=true
    @test all(axis -> axis.yscale[] === Makie.pseudolog10,signed.axes)
    keyword_ydata=LineCableModels.plot(reference,points[1];
        ydata=(R,),series_labels=("reference","candidate"),options...)
    @test length(keyword_ydata.axes)==4
    foreign=LineParameters(PhaseDomain,z,y,f;details=(coordinates=["a","b"],
        formulations=(schema_version=3,selections=(constitutive=(identifier=:default,),),
            assumptions=(equations=repeat("field equations ",100),))))
    foreign_plot=LineCableModels.plot(candidates,(R,);reference=foreign,options...)
    @test "Reference · method unavailable" in values(foreign_plot.addon_state.labels)
    @test all(label -> !occursin("field equations",label),values(foreign_plot.addon_state.labels))
    @test foreign_plot.addon_state.formulations.reference === details(foreign)

    @test_throws r"explicitly saved" LineCableModels.plot(publication,(R,);pair=(2,1),options...)
    @test_throws r"not retained" LineCableModels.plot(publication,(R,);band=(100.,200.),options...)
    @test_throws r"no retained samples" LineCableModels.plot(publication,(R,);band=:wide,options...)
    duplicated=ParametricResult(nothing,[points[1],points[1]],candidates.axes,(;))
    overlay=LineCableModels.plot(duplicated,(R,);reference,options...)
    @test all(length(filter(item -> item isa Makie.Lines,axis.scene.plots))==3 for axis in overlay.axes)
    multi=ParametricResult(nothing,[points[1],points[1],points[2],points[2]],
        (problems=[:one,:two],formulations=selections),(;))
    @test_throws r"select problem" LineCableModels.plot(multi,(R,);options...)
    chosen=LineCableModels.plot(multi,(R,);problem=2,options...)
    @test chosen.addon_state.formulations.problem==2
    references=ParametricResult(nothing,[reference,reference],(problems=[:one,:two],formulations=selections[1:1]),(;))
    paired=report(BenchmarkTableDefinition(pairing=[(1,1),(2,2),(1,3),(2,4)]),
        (reference=references,candidate=multi))
    a=LineCableModels.plot(paired,(R,);problem=1,options...)
    b=LineCableModels.plot(paired,(R,);problem=2,options...)
    @test first(first(a.axes).scene.plots).color[] == first(first(b.axes).scene.plots).color[]
    mktempdir() do root
        for controls in (false,true)
            page=LineCableModels.plot(publication,(R,);options...,controls,series_attributes)
            xlims!(first(page.axes),2,70)
            ylims!(first(page.axes),0,25)
            limits=[axis.finallimits[] for axis in page.axes]
            exported=export_svg(page;path=joinpath(root,"grid-$controls.svg"),open_file=false)
            @test isfile(exported)
            @test occursin("<svg",read(exported,String))
            @test [axis.finallimits[] for axis in page.axes]==limits
            @test Z(reference)==z
        end
    end
end

@testitem "Makie addons / owner-defined composite slots survive real reports and legends" tags=[:visual] begin
    using CairoMakie
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    import LineCableModels: description, formula_id, formulation_options
    const E=LineCableModels.Engine
    const IO=LineCableModels.ImportExport
    struct TestFormula{ID} <: AbstractFormulation end
    TestFormula(id::Symbol)=TestFormula{id}()
    formula_id(::Type{<:TestFormula{ID}}) where {ID}=ID
    formula_id(value::TestFormula)=formula_id(typeof(value))
    description(::Type{<:TestFormula};compact::Bool=false)=compact ? "Shared display" : "Test-owned scientific explanation"
    description(value::TestFormula;compact::Bool=false)=description(typeof(value);compact)
    formulation_options(::TestFormula)=(;)
    formulation_options(::Type{<:TestFormula},record::NamedTuple)=(;)
    Base.NamedTuple(value::TestFormula)=(identifier=formula_id(value),)
    Base.pairs(::Type{<:TestFormula};quantity=nothing)=pairs((east=TestFormula,west=TestFormula))
    struct TestOwner{T} <: AbstractFormulation
        selection::T
    end
    description(::Type{<:TestOwner};compact::Bool=false)="Test owner"
    description(value::TestOwner;compact::Bool=false)=description(typeof(value);compact)
    formula_id(::Type{<:TestOwner})=:test_owner
    formula_id(::TestOwner)=:test_owner
    formulation_options(::TestOwner)=(;)
    formulation_options(::Type{TestOwner},record::NamedTuple)=record.options
    description(::Type{TestOwner},::Val{:channel})="channel"
    Base.pairs(::Type{TestOwner};quantity=nothing)=pairs((channel=TestFormula,))
    Base.pairs(::Type{TestOwner},record::NamedTuple;quantity=nothing)=
        pairs(E.LineParametersFormulation,record;quantity,owner=TestOwner)
    Base.pairs(value::TestOwner;quantity=nothing)=pairs(TestOwner,
        (methods=(channel=value.selection,),requested=(channel=map(formulation_options,value.selection),),options=(;));quantity)

    a=Formulation(TestFormula,(west=:right,east=:left))
    b=Formulation(TestFormula,(east=:right,west=:left))
    @test keys(a)==(:east,:west)
    @test formula_id(a.east)!=formula_id(b.east)
    @test description(a.east;compact=true)==description(b.east;compact=true)
    selected,controls=IO.deserialize_value(Val(:formulation),TestFormula,map(NamedTuple,a),map(NamedTuple,a))
    saved=TestOwner=>(methods=(channel=selected,),requested=(channel=controls,),options=(;))
    for compact in (false,true)
        @test description([TestOwner(a)];compact)==description([saved];compact)
    end
    f=[1.,10.,100.]
    z=reshape(complex.(1.:12.,21.:32.),2,2,3)
    y=reshape(complex.(101.:112.,201.:212.),2,2,3)*1e-6
    reference=LineParameters(z,y,f;details=(coordinates=["a","b"],formulations=NamedTuple(LineCableModelsFEM()),))
    points=[LineParameters(k*z,k*y,f;details=(coordinates=["a","b"],)) for k in (1.,2.,3.)]
    choices=[TestOwner(a),TestOwner(b),saved]
    candidates=ParametricResult(nothing,points,(problems=[:one],formulations=choices),(;))
    artifact=report(BenchmarkTableDefinition((R,B);bands=(:all,)),(reference=reference,candidate=candidates))
    options=(backend=:cairo,display_plot=false,controls=false,open_export=false,length_unit=:base)
    for source in (candidates,artifact)
        page=LineCableModels.plot(source;ydata=(R,),formulations=[3,1],options...)
        names=[page.addon_state.labels[group] for group in page.addon_state.order]
        candidate_names=source===artifact ? names[2:end] : names
        @test candidate_names==[
            "F3 · channel(east)=Shared display; channel(west)=Shared display",
            "F1 · channel(east)=Shared display; channel(west)=Shared display"]
        curves=filter(item->item isa Makie.Lines,first(page.axes).scene.plots)
        @test length(curves)==(source===artifact ? 3 : 2)
        @test last.(curves[end-1][1][])≈3real.(z[1,1,:])
        @test last.(curves[end][1][])≈real.(z[1,1,:])
    end
    # This checks literal route completeness, not table==plot: both consumers
    # used to agree while silently dropping unchanged/default branches.
    internal=(inner=:default,outer=:default,transfer=:default)
    earth=(air=:default,earth=:Pollaczek1926,mixed=:default)
    physical=Formulation(internal_impedance=internal,earth_impedance=earth,
        earth_admittance=(air=:default,earth=:default,mixed=:default))
    for retained in (physical,IO.deserialize_value(Val(:formulation),NamedTuple(physical)))
        data=ParametricResult(nothing,points[1:1],(problems=[:one],formulations=[retained]),(;))
        for (request,names) in ((R,("internal Z(inner)=default","internal Z(outer)=default",
                "internal Z(transfer)=default","earth Z(air)=default",
                "earth Z(earth)=Pollaczek1926","earth Z(mixed)=default")),
                (B,("earth Y(air)=default","earth Y(earth)=default","earth Y(mixed)=default")))
            page=LineCableModels.plot(data;ydata=(request,),options...)
            label=only(values(page.addon_state.labels))
            @test all(occursin(name,label) for name in names)
            @test !occursin(request===R ? "earth Y" : "internal Z",label)
        end
    end
    @test Z(reference)==z && Y(reference)==y
end

@testitem "Makie addons / formulation legends follow the plotted family without collapsing entries" tags=[:visual] begin
    using CairoMakie
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    f = [1.0, 10.0, 100.0]
    tensor = fill(1.0+2im, 1, 1, 3)
    reference = LineParameters(PhaseDomain, tensor, 1e-6tensor, f;
        details=(coordinates=["a"],))
    ids = ((:default,:default), (:Pollaczek1926,:Pollaczek1926),
        (:Saad1996,:default), (:WedepohlWilcox1973,:default), (:Xue2018,:Xue2018))
    records = [(backend=:coaxial, requested=(
        earth_impedance=(identifier=z,), earth_admittance=(identifier=y,)),
        options=(reduce_bundle=false,)) for (z,y) in ids]
    original = deepcopy(records)
    candidates = ParametricResult(nothing,fill(reference,5),
        (problems=[:one],formulations=records),(;))
    baseline = (result=reference,metadata=(port_order=["a"],formulation=records[1],axes=nothing))
    artifact = report(BenchmarkTableDefinition(),(reference=baseline,candidate=candidates))
    options = (;backend=:cairo,display_plot=false,controls=false,open_export=false,
        length_unit=:base)
    for source in (candidates,artifact)
        extra = source === candidates ? (;reference) : (;)
        for quantity in (R,L,X,Z,G,C,B,Y)
            built = LineCableModels.plot(source; ydata=(quantity,),options...,extra...)
            pages = built isa UIPlot ? (built,) : built
            family = quantity in (R,L,X,Z) ? "Z" : "Y"
            for page in pages
                other = family == "Z" ? "Y" : "Z"
                names = [page.addon_state.labels[group] for group in page.addon_state.order]
                @test length(names) == 6
                @test all(!occursin("earth $other",label) for label in names)
                @test all(occursin("earth $family",label) for label in names[2:end])
                @test count(label -> occursin("=default",label),names[2:end]) == (family == "Z" ? 1 : 3)
                curves = filter(plot -> plot isa Makie.Lines,first(page.axes).scene.plots)
                @test length(curves) == 6
                @test all(curve -> curve[1][] == first(curves)[1][],curves)
                @test allunique([curve.color[] for curve in curves])
                @test length(last(only(page.legend.entrygroups[]))) == 6
            end
        end
        filtered = LineCableModels.plot(source; ydata=(G,),formulations=[4,2,1],options...,extra...)
        names = [filtered.addon_state.labels[group] for group in filtered.addon_state.order]
        expected = ["F4","F2","F1"]
        @test all(startswith.(names[2:end],expected))
        @test count(label -> occursin("=default",label),names[2:end]) == 2
        reordered=LineCableModels.plot(source;ydata=(R,),formulations=[3,1],options...,extra...)
        reordered_names=[reordered.addon_state.labels[group] for group in reordered.addon_state.order]
        @test startswith.(reordered_names[2:end],["F3","F1"])==[true,true]
        override = ("ref","a","b","c","d","e")
        custom = LineCableModels.plot(source; ydata=(R,G),series_labels=override,options...,extra...)
        @test all(page -> Set(values(page.addon_state.labels)) == Set(override),custom)
    end
    @test records == original
end
