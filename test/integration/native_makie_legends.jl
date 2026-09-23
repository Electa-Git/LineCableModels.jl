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

@testitem "Makie addons / independent scale arrangement and sibling spacing" tags=[:visual] setup=[TestFixtures] begin
    using CairoMakie
    design=TestFixtures.coaxial_design()
    options=(backend=:cairo,display_plot=false,controls=false,size=(1200,850),
        colorbar_attributes=(vertical=false,width=160,height=14))
    @testset "placement chooses default cells independently of bar orientation" begin
        p=preview(design;options...,legend_position=nothing,colorbar_position=:bottom)
        Makie.colorbuffer(p.figure)
        bounds=[bar.layoutobservables.computedbbox[] for bar in p.colorbars]
        @test maximum(b.origin[2] for b in bounds)-minimum(b.origin[2] for b in bounds)<=1
        @test all(bounds[i+1].origin[1]>bounds[i].origin[1]+bounds[i].widths[1] for i in 1:2)
    end
    @testset "moving and direct construction converge" begin
        p=preview(design;options...,legend_position=nothing,colorbar_position=:right)
        bars=copy(p.colorbars)
        figurecolorbars!(p;position=:bottom)
        Makie.colorbuffer(p.figure)
        @test all(a===b for (a,b) in zip(bars,p.colorbars))
        @test all(bar -> !bar.labelvisible[],bars)
        bounds=[bar.layoutobservables.computedbbox[] for bar in bars]
        @test maximum(b.origin[2] for b in bounds)-minimum(b.origin[2] for b in bounds)<=1
    end
    @testset "complete sibling bounds have a minimum skip" begin
        p=preview(design;options...,legend_position=:right,
            legend_attributes=(valign=:center,margin=(0,0,0,0)),colorbar_position=:right,
            colorbar_group_attributes=(valign=:center,margin=(0,0,0,0)))
        Makie.colorbuffer(p.figure)
        legend=p.legend.layoutobservables.computedbbox[]
        scales=p.addon_state.guides[(:colorbars,nothing)].layout[].layoutobservables.computedbbox[]
        skip=legend.origin[2]-(scales.origin[2]+scales.widths[2])
        @info "complete guide baseline" skip legend scales
        @test skip>=11
    end
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

@testitem "Makie addons / guide arrangement packing and option validation" tags=[:visual] begin
    using CairoMakie
    ext=Base.get_extension(LineCableModels,:LineCableModelsMakieExt)
    for (extents,gap,span,fractions) in (([20.,30.],12.,62.,[.5,.5]),
            ([20.,30.],12.,120.,[.25,.25]),([20.,30.],0.,100.,[0.,1.]),
            ([20.],25.,100.,[.5]),([20.,30.,10.],7.5,120.,[1.,0.,.5]))
        result=ext._addon_pack_guides(extents,span,gap,fractions)
        @test result.extent==sum(extents)+(length(extents)-1)*gap
        @test first(result.starts)>=0
        @test last(result.starts)+last(extents)<=result.span
        @test all(result.starts[i+1]>=result.starts[i]+extents[i]+gap for i in 1:length(extents)-1)
        if all(==(first(fractions)),fractions)
            @test first(result.starts)==first(fractions)*(span-result.extent)
        end
    end
    @test ext._addon_pack_guides(Float64[],100.,12.,Float64[]).extent==0
    @test ext._addon_pack_guides([20.,30.],100.,12.,[0.,1.]).starts==[0.,70.]
    @test ext._addon_pack_guides([20.,30.],30.,12.,[0.,1.]).span==62
    @test ext._addon_pack_guides([20.,30.,10.],200.,15.,fill(.5,3)).extent-
        ext._addon_pack_guides([20.,30.,10.],200.,5.,fill(.5,3)).extent==20
    @test ext._addon_guide_spacing((rowgap=2.5,))==(rowgap=2.5,colgap=12.)
    @test ext._addon_guide_spacing((colgap=3.5,),(rowgap=2.5,colgap=8.))==(rowgap=2.5,colgap=3.5)
    for value in (true,-1,Inf,NaN,(rowgap=true,),(rowgap=-1,),(colgap=Inf,),(extra=1,))
        @test_throws ArgumentError ext._addon_guide_spacing(value)
    end
    for value in ((layout=(true,3),),(layout=(1,2),),(layout=(0,3),),(rowgap=NaN,),
            (colgap=-1,),(rowgap=true,),(layout=:horizontal,),(horizontal_strip=true,))
        @test_throws ArgumentError ext._addon_colorbar_group_attributes(value,3)
    end
end

@testitem "Makie addons / guide arrangement reflows retained complete items" tags=[:visual] setup=[TestFixtures] begin
    using CairoMakie
    ext=Base.get_extension(LineCableModels,:LineCableModelsMakieExt)
    design=TestFixtures.coaxial_design()
    options=(backend=:cairo,display_plot=false,controls=false,size=(1400,950),
        legend_position=nothing,colorbar_attributes=(vertical=false,width=160,height=14),
        colorbar_group_attributes=(layout=(1,3),rowgap=12.5,colgap=16.5,margin=(0,0,0,0)))
    p=preview(design;options...,colorbar_position=:bottom,guide_spacing=(rowgap=14,colgap=18))
    direct=preview(design;options...,colorbar_position=:bottom)
    guide=p.addon_state.guides[(:colorbars,nothing)]
    bars=copy(p.colorbars)
    companions=[item.companion for item in guide.items]
    plots=copy(only(p.axes).scene.plots)
    Makie.limits!(only(p.axes),-.03,.025,-.025,.03)
    view=only(p.axes).targetlimits[]
    frame=only(p.axes).layoutobservables.computedbbox[].widths
    listeners=length(guide.subscriptions)
    coordinates()=[Tuple(Makie.GridLayoutBase.gridcontent(item.layout).span.rows)[1]=>
        Tuple(Makie.GridLayoutBase.gridcontent(item.layout).span.cols)[1] for item in guide.items if item.visible[]]
    for position in (:right,:bottom,:top,:left,:bottom)
        figurecolorbars!(p;position)
        Makie.colorbuffer(p.figure)
        @test size(guide.layout[])==(1,3)
        @test coordinates()==[1=>1,1=>2,1=>3]
        @test all(a===b for (a,b) in zip(bars,p.colorbars))
        @test all(a===b for (a,b) in zip(companions,[item.companion for item in guide.items]))
        @test only(p.axes).targetlimits[]==view
        @test all(isapprox.(only(p.axes).layoutobservables.computedbbox[].widths,frame;atol=1))
        @test only(p.axes).scene.plots==plots
        @test all(bar -> bar.width[]==160 && bar.height[]==14 && !bar.vertical[],bars)
        @test length(guide.subscriptions)==listeners
    end
    # Compare local geometry; independent windows retain their own frame sizes.
    Makie.colorbuffer(direct.figure)
    function local_boxes(q)
        g=q.addon_state.guides[(:colorbars,nothing)]
        origin=g.layout[].layoutobservables.computedbbox[].origin
        [(Tuple(item.layout.layoutobservables.computedbbox[].origin-origin),
            Tuple(item.layout.layoutobservables.computedbbox[].widths)) for item in g.items]
    end
    @test all(all(isapprox.(a[1],b[1];atol=1)) && all(isapprox.(a[2],b[2];atol=1))
        for (a,b) in zip(local_boxes(p),local_boxes(direct)))
    @test all(bar -> !bar.labelvisible[],bars)
    @test all(item -> item.companion.blockscene.visible[],guide.items)
    figurecolorbars!(p;position=:right,group_attributes=(layout=nothing,))
    @test size(guide.layout[])==(3,1)
    figurecolorbars!(p;position=:bottom)
    @test size(guide.layout[])==(1,3)
    figurecolorbars!(p;group_attributes=(layout=(2,5),))
    @test size(guide.layout[])==(1,3) # spare guide capacity is not a matrix domain
    figurecolorbars!(p;group_attributes=(layout=(3,1),),vertical=true)
    @test size(guide.layout[])==(3,1)
    @test all(bar -> bar.vertical[] && bar.width[]==160 && bar.height[]==14,bars)
    figurecolorbars!(p;group_attributes=(layout=(1,3),),vertical=false)
    bars[1].label[]="A deliberately long retained material-property label"
    bars[1].labelsize[]=23
    bars[1].labelcolor[]=:red
    bars[2].blockscene.visible[]=false
    Makie.colorbuffer(p.figure)
    @test size(guide.layout[])==(1,2)
    @test !guide.items[2].companion.blockscene.visible[]
    @test Makie.GridLayoutBase.gridcontent(guide.items[2].layout).parent===nothing
    @test length(p.colorbars)==3
    for position in (nothing,:right,:bottom)
        figurecolorbars!(p;position)
        @test !bars[2].blockscene.visible[]
        @test !guide.items[2].companion.blockscene.visible[]
    end
    bars[2].blockscene.visible[]=true
    @test coordinates()==[1=>1,1=>2,1=>3]
    @test bars[1].label[]=="A deliberately long retained material-property label"
    @test bars[1].labelsize[]==23
    @test Makie.to_color(bars[1].labelcolor[])==Makie.to_color(:red)
    @test length(guide.subscriptions)==listeners
    @test only(p.axes).targetlimits[]==view

    before=(position=guide.position[],group=p.addon_state.colorbar_group_attributes[],
        spacing=p.addon_state.guide_spacing[],size=p.figure.scene.viewport[],
        attributes=guide.attributes[],boxes=local_boxes(p))
    for invalid in ((group_attributes=(layout=(1,2),),),(guide_spacing=-1,),
            (group_attributes=(bogus=1,),),(group_attributes=(colgap=true,),),
            (vertical=:horizontal,),(width=NaN,),(position=(2,2),),(unknown_attribute=1,))
        @test_throws ArgumentError figurecolorbars!(p;position=:left,invalid...)
        @test guide.position[]==before.position
        @test p.addon_state.colorbar_group_attributes[]==before.group
        @test p.addon_state.guide_spacing[]==before.spacing
        @test p.figure.scene.viewport[]==before.size
        @test guide.attributes[]==before.attributes
        @test local_boxes(p)==before.boxes
    end
    figurecolorbars!(p;group_attributes=(rowgap=nothing,colgap=nothing),guide_spacing=(colgap=19.5,))
    @test p.addon_state.guide_spacing[]==(rowgap=14.,colgap=19.5)
    @test guide.layout[].default_colgap==Makie.Fixed(19.5)
    figurelegend!(p;guide_spacing=(rowgap=15.5,))
    @test p.addon_state.guide_spacing[]==(rowgap=15.5,colgap=19.5)
    @test guide.layout[].default_rowgap==Makie.Fixed(15.5)
    # Structural reflow must not replay constructor appearance over native edits.
    guide.layout[].halign[]=:left
    guide.layout[].alignmode[]=Makie.Outside(3,4,5,6)
    appearance=guide.layout[].alignmode[]
    figurecolorbars!(p;group_attributes=(layout=(1,3),))
    figurelegend!(p;guide_spacing=(rowgap=15.5,))
    @test guide.layout[].halign[].x==0
    @test guide.layout[].alignmode[]==appearance
    geometry=local_boxes(p)
    previous_group=p.addon_state.colorbar_group_attributes[]
    previous_size=p.figure.scene.viewport[]
    @test_throws ArgumentError figurecolorbars!(p;group_attributes=(width=1,))
    @test p.addon_state.colorbar_group_attributes[]==previous_group
    @test p.figure.scene.viewport[]==previous_size
    @test local_boxes(p)==geometry
    tickformat=bars[1].tickformat[]
    ticks=[bar.ticks[] for bar in bars]
    @test_throws ErrorException figurecolorbars!(p;position=:right,
        ticks=[1.,2.],
        tickformat=values -> error("injected native tick formatter failure"))
    @test all(bar -> bar.tickformat[]===tickformat,bars)
    @test [bar.ticks[] for bar in bars]==ticks
    @test guide.position[]===:bottom
    @test only(p.axes).targetlimits[]==view
    @test !isempty(Makie.colorbuffer(p.figure))
end

@testitem "Makie addons / guide arrangement measures visible scales and independent gaps" tags=[:visual] setup=[TestFixtures] begin
    using CairoMakie
    design=TestFixtures.coaxial_design()
    options=(backend=:cairo,display_plot=false,controls=false,size=(1500,1100),
        colorbar_attributes=(vertical=false,width=160,height=14))
    p=preview(design;options...,legend_position=:bottom,legend_attributes=(halign=.25,margin=(0,0,0,0)),
        colorbar_position=:bottom,colorbar_group_attributes=(layout=(1,3),colgap=16.5,halign=.25,margin=(0,0,0,0)),
        guide_spacing=(rowgap=14.5,colgap=18.5),guide_gap=8)
    guide=p.addon_state.guides[(:colorbars,nothing)]
    # Native block scenes use the figure's logical pixels as their data coordinates.
    # Measure rendered glyphs, including native endpoint text, not only allocations.
    function visible_bounds(item)
        boxes=[Makie.boundingbox(item.bar.blockscene,plot -> !plot.visible[],:data)]
        item.companion.blockscene.visible[] && push!(boxes,
            Makie.boundingbox(item.companion.blockscene,plot -> !plot.visible[],:data))
        low=ntuple(i -> minimum(box -> box.origin[i],boxes),2)
        high=ntuple(i -> maximum(box -> box.origin[i]+box.widths[i],boxes),2)
        (;low,high)
    end
    function check_items(p,dimension,minimum_gap)
        Makie.colorbuffer(p.figure)
        g=p.addon_state.guides[(:colorbars,nothing)]
        active=filter(item -> item.visible[],g.items)
        bounds=visible_bounds.(active)
        for (item,box) in zip(active,bounds)
            reserved=item.layout.layoutobservables.computedbbox[]
            @test all(box.low[i]>=reserved.origin[i]-1 for i in 1:2)
            @test all(box.high[i]<=reserved.origin[i]+reserved.widths[i]+1 for i in 1:2)
        end
        for i in 1:length(bounds)-1
            gap=dimension==1 ? bounds[i+1].low[1]-bounds[i].high[1] : bounds[i].low[2]-bounds[i+1].high[2]
            @test gap>=minimum_gap-1
        end
        bounds
    end
    boxes=check_items(p,1,16.5)
    # With identical font sizes these horizontal items need only their visible
    # height, not a second copy of the native tick/property-label protrusion.
    @test first(guide.items).layout.layoutobservables.computedbbox[].widths[2]<=
        first(boxes).high[2]-first(boxes).low[2]+3
    legend=p.legend.layoutobservables.computedbbox[]
    @test minimum(b -> b.low[1],boxes)-(legend.origin[1]+legend.widths[1])>=17.5
    @test p.addon_state.guide_gap==(8.,8.,8.,8.)
    dimensions=[(bar.width[],bar.height[]) for bar in p.colorbars]
    inner=guide.layout[].addedcolgaps |> copy
    figurelegend!(p;guide_spacing=(colgap=30.5,))
    @test guide.layout[].addedcolgaps==inner
    @test [(bar.width[],bar.height[]) for bar in p.colorbars]==dimensions
    boxes=check_items(p,1,16.5)
    legend=p.legend.layoutobservables.computedbbox[]
    @test minimum(b -> b.low[1],boxes)-(legend.origin[1]+legend.widths[1])>=29.5
    p.colorbars[2].label[]="A much longer physical scale caption"
    p.colorbars[2].labelsize[]=25
    check_items(p,1,16.5)
    p.colorbars[1].spinewidth[]=4
    check_items(p,1,16.5)
    figurecolorbars!(p;position=:right,group_attributes=(layout=(3,1),rowgap=12.5,valign=.25))
    figurelegend!(p;position=:right,valign=.25,guide_spacing=(rowgap=14.5,))
    boxes=check_items(p,2,12.5)
    legend=p.legend.layoutobservables.computedbbox[]
    @test legend.origin[2]-maximum(b -> b.high[2],boxes)>=13.5
    for position in (nothing,:right)
        figurelegend!(p;position)
        check_items(p,2,12.5)
    end
    # Zero or one visible outer sibling has no spacing contribution.
    figurecolorbars!(p;position=nothing)
    Makie.colorbuffer(p.figure)
    size_before=p.figure.scene.viewport[].widths
    legend_before=p.legend.layoutobservables.computedbbox[]
    figurelegend!(p;guide_spacing=0)
    figurelegend!(p;guide_spacing=70)
    Makie.colorbuffer(p.figure)
    @test all(isapprox.(size_before,p.figure.scene.viewport[].widths;atol=1))
    @test all(isapprox.(legend_before.widths,p.legend.layoutobservables.computedbbox[].widths;atol=1))
    @test p.legend.margin[]==(0,0,0,0)
    @test p.addon_state.guide_gap==(8.,8.,8.,8.)
    # Orientation changes the native bar, never its sibling cells. Measure
    # vertical bars in both arrangements as well as the horizontal cases above.
    figurelegend!(p;position=nothing)
    guide.items[2].bar.label[]="Vertical scale"
    for (layout,dimension,gap) in (((1,3),1,16.5),((3,1),2,12.5))
        figurecolorbars!(p;position=:right,group_attributes=(;layout),
            vertical=true,width=14,height=160)
        @test size(guide.layout[])==layout
        check_items(p,dimension,gap)
    end
end

@testitem "Makie addons / scale bars share edges and baselines independently of labels" tags=[:visual] setup=[TestFixtures] begin
    using CairoMakie
    design=TestFixtures.coaxial_design()
    options=(backend=:cairo,display_plot=false,controls=false,
        colorbar_attributes=(vertical=false,width=160,height=14))
    p=preview(design;options...,colorbar_position=:right)
    bars=copy(p.colorbars)
    guide=p.addon_state.guides[(:colorbars,nothing)]
    views=[axis.targetlimits[] for axis in p.axes]
    listeners=length(guide.subscriptions)
    function aligned(p,dimension,indices=eachindex(p.colorbars))
        Makie.colorbuffer(p.figure)
        bounds=[p.colorbars[i].layoutobservables.computedbbox[] for i in indices if p.colorbars[i].blockscene.visible[]]
        for edge in (box -> box.origin[dimension],box -> box.origin[dimension]+box.widths[dimension])
            values=edge.(bounds)
            @test maximum(values)-minimum(values)<=1
        end
    end
    aligned(p,1)
    function local_origins(p)
        Makie.colorbuffer(p.figure)
        group=p.addon_state.guides[(:colorbars,nothing)].layout[]
        [bar.layoutobservables.computedbbox[].origin-group.layoutobservables.computedbbox[].origin for bar in p.colorbars]
    end
    right=local_origins(p)
    direct=preview(design;options...,colorbar_position=:bottom,
        colorbar_group_attributes=(layout=(1,3),))
    figurecolorbars!(p;position=:bottom,group_attributes=(layout=(1,3),))
    @test all(all(isapprox.(a,b;atol=1)) for (a,b) in zip(local_origins(p),local_origins(direct)))
    figurecolorbars!(p;position=:right,group_attributes=(layout=(3,1),))
    @test all(all(isapprox.(a,b;atol=1)) for (a,b) in zip(right,local_origins(p)))
    # Native caption and endpoint font changes must reserve space without
    # independently recentering the unequal decorated bar footprints.
    bars[1].labelsize[]=26
    bars[1].labelfont[]=:bold
    bars[2].ticklabelsize[]=23
    bars[2].ticklabelfont[]=:bold
    positions=first(bars[3].ticks[])
    bars[3].ticks[]=(positions,["long endpoint $(value)" for value in positions])
    aligned(p,1)
    for halign in (:left,:center,:right)
        figurecolorbars!(p;group_attributes=(;halign))
        aligned(p,1)
    end
    bars[2].blockscene.visible[]=false
    aligned(p,1)
    bars[2].blockscene.visible[]=true
    aligned(p,1)
    for position in (:bottom,:right,:bottom)
        figurecolorbars!(p;position,group_attributes=(layout=position===:bottom ? (1,3) : (3,1),))
        aligned(p,position===:bottom ? 2 : 1)
    end
    figurecolorbars!(p;group_attributes=(layout=(2,2),))
    aligned(p,1,(1,3))
    aligned(p,2,(1,2))
    @test all(a===b for (a,b) in zip(bars,p.colorbars))
    @test [axis.targetlimits[] for axis in p.axes]==views
    @test length(guide.subscriptions)==listeners
end

@testitem "Makie addons / guide arrangement membership native edits and export" tags=[:visual] setup=[TestFixtures] begin
    using CairoMakie
    ext=Base.get_extension(LineCableModels,:LineCableModelsMakieExt)
    design=TestFixtures.coaxial_design()
    options=(backend=:cairo,display_plot=false,controls=false,size=(1400,1000),open_export=false)
    hidden=preview(design;options...,colorbar_position=nothing)
    original_content=copy(hidden.figure.content)
    original_scenes=copy(hidden.figure.scene.children)
    original_size=hidden.figure.scene.viewport[]
    @test_throws ErrorException figurecolorbars!(hidden;position=:bottom,ticks=[1.,2.],
        tickformat=values -> error("injected native constructor failure"))
    @test hidden.figure.content==original_content
    @test hidden.figure.scene.children==original_scenes
    @test hidden.figure.scene.viewport[]==original_size
    @test hidden.addon_state.guides[(:colorbars,nothing)].position[]===nothing
    @test isempty(hidden.colorbars)
    @test length(figurecolorbars!(hidden;position=:bottom))==3
    for count in 0:3
        p=preview(design;options...,legend_position=nothing,colorbar_position=nothing)
        p.addon_state=merge(p.addon_state,(color_scales=p.addon_state.color_scales[1:count],))
        figurecolorbars!(p;position=:bottom,group_attributes=(layout=(2,4),),guide_spacing=12.5)
        Makie.colorbuffer(p.figure)
        @test length(p.colorbars)==count
        guide=p.addon_state.guides[(:colorbars,nothing)]
        if count==0
            @test isempty(p.addon_state.guide_docks)
        else
            @test size(guide.layout[])==(1,count)
            for bar in p.colorbars
                bar.blockscene.visible[]=false
            end
            @test isempty(p.addon_state.guide_docks)
            @test length(p.colorbars)==count
            @test all(item -> !item.companion.blockscene.visible[],guide.items)
            figurecolorbars!(p;position=nothing)
            figurecolorbars!(p;position=:right)
            @test length(p.colorbars)==count
            @test isempty(p.addon_state.guide_docks)
            foreach(bar -> bar.blockscene.visible[]=true,p.colorbars)
            @test size(guide.layout[])==(1,count)
        end
    end
    p=preview(design;options...,legend_position=:right,colorbar_position=:bottom,
        colorbar_group_attributes=(layout=(1,3),colgap=16.5),guide_spacing=14.5)
    guide=p.addon_state.guides[(:colorbars,nothing)]
    bars=copy(p.colorbars)
    bars[1].width[]=213.5
    bars[1].vertical[]=true
    @test size(guide.layout[])==(1,3)
    @test bars[1].width[]==213.5
    @test bars[1].height[]==ext._ADDON_COLORBAR_DOCK_LENGTH
    bars[1].vertical[]=false
    @test bars[1].width[]==213.5
    @test bars[1].height[]==Makie.Auto()
    @test isapprox(bars[1].layoutobservables.computedbbox[].widths[2],bars[1].size[];atol=1)
    bars[2].blockscene.visible[]=false
    bars[3].label[]="Exported scale caption"
    bars[3].labelsize[]=21
    Makie.limits!(only(p.axes),-.024,.036,-.025,.035)
    view=only(p.axes).targetlimits[]
    native=(size=p.figure.scene.viewport[],labels=[bar.label[] for bar in bars],
        fonts=[bar.labelfont[] for bar in bars],labelsize=[bar.labelsize[] for bar in bars],
        visibility=[item.visible[] for item in guide.items],
        scenes=[bar.blockscene.visible[] for bar in bars],subscriptions=length(guide.subscriptions))
    mktempdir() do directory
        for theme in (:default,:publication)
            path=export_svg(p;path=joinpath(directory,string(theme)*".svg"),theme,open_file=false)
            @test isfile(path)
            @test occursin("<svg",read(path,String))
            @test_throws ErrorException ext._addon_export_presentation!(p,theme) do
                @test guide.items[2].visible[]==false
                @test !guide.items[2].companion.blockscene.visible[]
                Makie.colorbuffer(p.figure)
                boxes=[Makie.boundingbox(bar.blockscene,plot -> !plot.visible[],:data)
                    for bar in bars if bar.blockscene.visible[]]
                @test length(boxes)==2
                @test boxes[2].origin[1]-(boxes[1].origin[1]+boxes[1].widths[1])>=15.5
                error("injected SVG writer failure")
            end
            @test p.figure.scene.viewport[]==native.size
            @test only(p.axes).targetlimits[]==view
            @test [bar.label[] for bar in bars]==native.labels
            @test [bar.labelfont[] for bar in bars]==native.fonts
            @test [bar.labelsize[] for bar in bars]==native.labelsize
            @test [item.visible[] for item in guide.items]==native.visibility
            @test [bar.blockscene.visible[] for bar in bars]==native.scenes
            @test length(guide.subscriptions)==native.subscriptions
            @test all(a===b for (a,b) in zip(bars,p.colorbars))
        end
    end
    # The shared setting reaches result panels and a caller-owned native canvas.
    raw=TestFixtures.two_conductor_results(;frequencies=[1.,10.,100.])
    result=LineCableModels.plot(raw,raw;ydata=((R,1,1:2,:),),layout=(1,2),
        backend=:cairo,display_plot=false,controls=false,guide_spacing=(rowgap=7.5,),
        panel_legends=(1,1)=>:right)
    @test result.addon_state.guide_spacing[]==(rowgap=7.5,colgap=12.)
    panel=result.panel_legends[(1,1)]
    figurelegend!(result;guide_spacing=(colgap=17.5,))
    @test result.panel_legends[(1,1)]===panel
    @test result.addon_state.guide_spacing[]==(rowgap=7.5,colgap=17.5)
    @test_throws ArgumentError panellegend!(result,(1,1);guide_spacing=5)
    canvas=LineCableModels.plotwindow(;title="Native guide spacing",backend=:cairo,display_plot=false,controls=false,guide_spacing=14.5) do grid
        axis=Axis(grid[1,1]);lines!(axis,[1.,2.,3.],[2.,3.,4.];label="native curve")
        axis
    end
    @test canvas.addon_state.guide_spacing[]==(rowgap=14.5,colgap=14.5)
    figurelegend!(canvas;position=:right,valign=.5)
    panellegend!(canvas,1;position=:right,valign=.5)
    Makie.colorbuffer(canvas.figure)
    a=canvas.legend.layoutobservables.computedbbox[]
    b=canvas.panel_legends[1].layoutobservables.computedbbox[]
    @test a.origin[2]-(b.origin[2]+b.widths[2])>=13.5
    original_view=only(canvas.axes).targetlimits[]
    canvas.legend.blockscene.visible[]=false
    @test !canvas.legend.blockscene.visible[]
    canvas.legend.blockscene.visible[]=true
    @test only(canvas.axes).targetlimits[]==original_view
    scales=LineCableModels.show_material_scale(;backend=:cairo,display_plot=false,controls=false,guide_spacing=9.5)
    @test size(scales.addon_state.guides[(:colorbars,nothing)].layout[])==(3,1)
    @test scales.addon_state.guide_spacing[]==(rowgap=9.5,colgap=9.5)
    collection=preview(fill(design,2);options...,layout=(1,2),colorbar_position=:bottom,
        colorbar_group_attributes=(layout=(1,3),),guide_spacing=10.5)
    @test length(collection.axes)==2
    @test size(collection.addon_state.guides[(:colorbars,nothing)].layout[])==(1,3)
end
