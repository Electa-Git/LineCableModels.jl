@testitem "Makie addons / construction-time panel legends preserve source ownership" tags=[:visual] begin
    using CairoMakie

    frequency = [1.0, 10.0, 100.0]
    impedance = [complex(i + 2j + k, i + j + 2k) for i in 1:2, j in 1:2, k in 1:3]
    admittance = impedance .* 1e-6
    reference = LineParameters(copy(impedance), copy(admittance), frequency)
    candidate = LineParameters(2impedance, 2admittance, frequency)
    sources = (; reference, candidate)
    requests = ((R, 1, 1:2, :),)
    options = (; backend=:cairo, display_plot=false, controls=false,
        layout=(1, 2), panel_titles=("Self resistance", "Mutual resistance"),
        series_labels=("baseline", "alternative"), legend_position=nothing,
        legend_overflow=:show_all, length_unit=:base, quantity_units=:base,
        freq_unit=:base, clip=false)

    plot = Makie.plot(sources, requests; options...,
        panel_legends=(1, 1) => :right)
    @test Set(keys(plot.panel_legends)) == Set(((1, 1),))
    @test plot.panel_legends[(1, 1)] isa Makie.Legend
    @test plot.addon_state.guides[(:legend,(1,1))].position[] == :right
    @test plot.panel_legends[(1, 1)].orientation[] == :vertical

    local_labels = ("local baseline", "local alternative")
    inside = Makie.plot(sources, requests; options...,
        panel_legends=(
            (1, 1) => (position=:inside, halign=:left, valign=:bottom,
                title="Local sources", legend_labels=local_labels,
                overflow=:show_all, labelsize=13),
            (1, 2) => false,
        ))
    legend = inside.panel_legends[(1, 1)]
    @test Set(keys(inside.panel_legends)) == Set(((1, 1),))
    @test legend.halign[] == :left
    @test legend.valign[] == :bottom
    @test legend.labelsize[] == 13
    @test first(only(legend.entrygroups[])) == "Local sources"
    @test Set(values(inside.addon_state.panel_data[(1, 1)].labels)) == Set(local_labels)
    @test [entry.label[] for entry in last(only(legend.entrygroups[]))] ==
        collect(local_labels)
    @test Set(values(inside.addon_state.labels)) == Set(("baseline", "alternative"))
    @test Set(values(inside.addon_state.panel_data[(1, 2)].labels)) ==
        Set(("baseline", "alternative"))

    moved = panellegend!(inside, (1, 1); position=:bottom, overflow=:show_all)
    @test moved === inside.panel_legends[(1, 1)]
    @test moved === legend
    @test moved.orientation[] == :vertical # native row-major grid
    @test moved.nbanks[] >= 1
    @test first(only(moved.entrygroups[])) == "Local sources"
    @test inside.addon_state.guides[(:legend,(1,1))].position[] == :bottom
    @test panellegend!(inside, (1, 1); position=nothing) === nothing
    @test !haskey(inside.panel_legends, (1, 1))
    @test inside.addon_state.guides[(:legend,(1,1))].position[]===nothing

    configured = Makie.plot(sources, requests; options...,
        panel_legends=Dict(
            (1, 1) => nothing,
            (1, 2) => (position=:top, overflow=:show_all,
                legend_labels=Dict("baseline"=>"renamed"), title="Reactance sources"),
        ))
    @test Set(keys(configured.panel_legends)) == Set(((1, 2),))
    @test configured.panel_legends[(1, 2)].orientation[] == :vertical
    @test configured.panel_legends[(1, 2)].nbanks[] >= 1
    @test first(only(configured.panel_legends[(1, 2)].entrygroups[])) ==
        "Reactance sources"
    @test Set(values(configured.addon_state.panel_data[(1, 2)].labels)) ==
        Set(("renamed", "alternative"))
    @test [entry.label[] for entry in last(only(
        configured.panel_legends[(1, 2)].entrygroups[]))] == ["renamed", "alternative"]
    @test Set(values(configured.addon_state.labels)) == Set(("baseline", "alternative"))
    panellegend!(configured, (1, 2); position=:left, overflow=:show_all)
    @test configured.addon_state.guides[(:legend,(1,2))].position[] == :left
    @test configured.panel_legends[(1, 2)].orientation[] == :vertical

    for rendered in (plot, inside, configured)
        Makie.colorbuffer(rendered.figure)
        @test [axis.title[] for axis in rendered.axes] == ["Self resistance", "Mutual resistance"]
        for (index, axis) in enumerate(rendered.axes)
            curves = filter(item -> item isa Makie.Lines, axis.scene.plots)
            @test length(curves) == 2
            expected = real.(impedance[1,index,:])
            @test last.(curves[1][1][]) ≈ expected
            @test last.(curves[2][1][]) ≈ 2expected
            @test first.(curves[1][1][]) ≈ frequency
        end
    end
    @test Z(reference) == impedance
    @test Y(reference) == admittance
    @test Z(candidate) == 2impedance
    @test Y(candidate) == 2admittance

    for requested in (42, ((0, 1) => :right), ((true, 1) => :right),
            ((1, 1) => 42), ((1, 1) => (position=:inside, anchor=(:left,))))
        @test_throws ArgumentError Makie.plot(sources, requests; options...,
            panel_legends=requested)
    end
    @test_throws ArgumentError Makie.plot(sources, requests; options...,
        panel_legends=(2, 1) => :right)
end

@testitem "Makie addons / guide restoration retains native edits and common placement" tags=[:visual] setup=[TestFixtures] begin
    using CairoMakie
    design=TestFixtures.coaxial_design()
    p=preview(design;backend=:cairo,display_plot=false,controls=false,
        legend_position=:right,colorbar_position=:right)
    legend=p.legend
    bars=copy(p.colorbars)
    legend.labelsize[]=17
    bars[1].labelsize[]=19
    bars[1].label[]="Edited native scale"
    axis=only(p.axes)
    Makie.limits!(axis,-0.04,0.03,-0.02,0.04)
    view=axis.targetlimits[]
    bars[2].blockscene.visible[]=false
    for _ in 1:3
        @test figurelegend!(p;position=nothing)===nothing
        @test isempty(figurecolorbars!(p;position=nothing))
        @test figurelegend!(p;position=:left)===legend
        restored=figurecolorbars!(p;position=:left,
            group_attributes=(valign=0.25,margin=(2,3,4,5)))
        @test all(a===b for (a,b) in zip(restored,bars))
        @test !bars[2].blockscene.visible[]
        @test bars[1].blockscene.visible[]
        @test axis.targetlimits[]==view
        @test legend.labelsize[]==17
        @test bars[1].labelsize[]==19
        @test bars[1].label[]=="Edited native scale"
        @test length(p.addon_state.guide_docks)==1
    end
    companion=filter(p.figure.content) do block
        block isa Makie.Label && block.text[]=="Edited native scale"
    end
    @test length(companion)==1
    @test only(companion).fontsize[]==19
    @test p.addon_state.shell.body.colsizes[3]==Makie.Fixed(0)
    @test_throws ArgumentError figurecolorbars!(p;position=(2,2))
    @test_throws ArgumentError paneltitle!(p,(9,9),"absent")
    @test axis.targetlimits[]==view

    raw=TestFixtures.two_conductor_results(;frequencies=[1.,10.,100.])
    q=Makie.plot(raw,raw,raw,raw; ydata=((R,1,1,:),),
        series_labels=("first complete description", "second complete description",
            "third complete description", "fourth complete description"),
        backend=:cairo,display_plot=false,controls=false,legend_position=:right,fig_size=(1200,700))
    original=q.legend
    figurelegend!(q;position=:bottom)
    @test q.legend===original
    @test q.legend.orientation[]===:vertical
    @test q.legend.nbanks[]>1
    @test length(last(only(q.legend.entrygroups[])))==4
    q.legend.nbanks[]=1
    figurelegend!(q;position=:top)
    @test q.legend.nbanks[]==1 # an explicit native edit releases automatic wrapping
    label=Makie.rich("native ",Makie.rich("label";color=:red))
    figuretitle!(q,label)
    @test q.title.text[]===label
    paneltitle!(q,(1,1),label)
    @test only(q.axes).title[]===label
end

@testitem "Makie addons / guide creation order shares complete scale geometry" tags=[:visual] setup=[TestFixtures] begin
    using CairoMakie
    design=TestFixtures.coaxial_design()
    options=(backend=:cairo,display_plot=false,controls=false,size=(1000,650),
        legend_position=nothing,colorbar_position=nothing)
    first_legend=preview(design;options...)
    first_scale=preview(design;options...)
    for (p,legend_first) in ((first_legend,true),(first_scale,false))
        if legend_first
            figurelegend!(p;position=:right,valign=0.5)
            figurecolorbars!(p;position=:right,group_attributes=(valign=0.5,))
        else
            figurecolorbars!(p;position=:right,group_attributes=(valign=0.5,))
            figurelegend!(p;position=:right,valign=0.5)
        end
        @test length(p.colorbars)==3
        @test p.legend isa Makie.Legend
        @test length(p.addon_state.guide_docks)==1
        @test all(isfinite,only(p.axes).targetlimits[].widths)
    end
    @test all(isapprox.(first_legend.figure.scene.viewport[].widths,first_scale.figure.scene.viewport[].widths;atol=1))
    @test all(isapprox.(only(first_legend.axes).scene.viewport[].widths,only(first_scale.axes).scene.viewport[].widths;atol=1))
    @test all(isapprox.(first_legend.legend.layoutobservables.computedbbox[].origin,
        first_scale.legend.layoutobservables.computedbbox[].origin;atol=1))
    hidden=preview(design;backend=:cairo,display_plot=false,controls=false,display_colorbars=false)
    @test isempty(hidden.colorbars)
    @test length(figurecolorbars!(hidden;position=:bottom))==3
    figurelegend!(hidden;position=nothing)
    figurecolorbars!(hidden;position=:right,group_attributes=(valign=0.25,))
    group=hidden.addon_state.guides[(:colorbars,nothing)].layout[]
    bounds=group.layoutobservables.computedbbox[]
    frame=hidden.addon_state.inside_bbox[]
    @test isapprox(bounds.origin[2],frame.origin[2]+0.25*(frame.widths[2]-bounds.widths[2]);atol=2)
end
