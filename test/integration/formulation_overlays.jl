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
    series_attributes=((marker=:circle, markersize=8), (;), (linestyle=:dash,))
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
        @test page.addon_state.formulations.records==selections
    end
    filtered=LineCableModels.plot(publication,(R,);formulations=2,options...)
    for (all_axis,filtered_axis) in zip(first(plots).axes,filtered.axes)
        all_curves=filter(item -> item isa Makie.Lines,all_axis.scene.plots)
        selected=filter(item -> item isa Makie.Lines,filtered_axis.scene.plots)
        @test length(selected)==2
        @test selected[2].color[] == all_curves[3].color[]
    end
    @test Set(values(filtered.addon_state.labels))==Set(expected_labels[[1,3]])
    foreign=LineParameters(PhaseDomain,z,y,f;details=(coordinates=["a","b"],
        formulations=(schema_version=3,selections=(constitutive=(identifier=:default,),),
            assumptions=(equations=repeat("field equations ",100),))))
    foreign_plot=LineCableModels.plot(candidates,(R,);reference=foreign,options...)
    @test "Reference F1" in values(foreign_plot.addon_state.labels)
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
