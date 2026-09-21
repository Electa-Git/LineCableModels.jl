@testitem "Makie addons / completed overlays preserve matrices, identities and exports" tags=[:visual] begin
    using CairoMakie
    using LineCableModels.Engine: retain_gridpoint, completed_formulation
    using LineCableModels.Grammar: gridpoint_id, observation_labels
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    f=[1.,10.,100.]
    z=reshape(complex.(collect(1.:12.),collect(21.:32.)),2,2,3)
    y=1e-6im .* z
    source_id=gridpoint_id().source_id
    choices=[Formulation(earth_impedance=:saad1996,earth_admittance=:unified),
        Formulation(earth_impedance=:xue2018,earth_admittance=:xue2018)]
    points=[retain_gridpoint(LineParameters(k*z,k*y,f),gridpoint_id(;source_id,formulation_index=index);
        fields=merge(completed_formulation(choices[index]),(inputs=(radius=.01,resistivity=100.),coordinates=["a","b"])))
        for (index,k) in enumerate((2,3))]
    reference=LineParameters(z,y,f)
    candidates=ParametricResult(nothing,points,(problems=[:unavailable],formulations=choices),ComputationDetails())
    completed=LineCableModels.Engine.compare(reference,candidates,[Z,Y,R,L,G,C];bands=(:all,:wide))
    timings=[(candidate_id=point.details.data.gridpoint,seconds=Float64(index),scope=:compute_call_wall,
        context=(id=:benchmark_title_probe,)) for (index,point) in enumerate(points)]
    observed=observables(candidates;comparisons=completed,timings,length_unit=:base,clip=false)
    observed_reference=ObservedResult(reference;length_unit=:base,clip=false)
    artifact=report(BenchmarkTableDefinition(bands=(:all,:wide)),observed;reference=observed_reference)
    @test artifact.tables.execution.seconds==[1.,2.]
    @test artifact.tables.execution.candidate_formulation==[1,2]
    @test length.(getproperty.(observed,:errors))==[12,12]
    @test all(error -> error.reference_id==observed_reference.gridpoint.id,Iterators.flatten(getproperty.(observed,:errors)))
    options=(backend=:cairo,display_plot=false,controls=false,open_export=false)
    pages=LineCableModels.plot(artifact,(Z,);options...)
    @test length(pages)==2
    @test getproperty.(pages,:export_name)==["benchmark_title_probe — Series resistance","benchmark_title_probe — Series reactance"]
    @test sum(length(page.axes) for page in pages)==8
    for (transform,page) in zip((real,imag),pages)
        @test length(page.addon_state.observed)==3
        @test any(startswith("Reference"),values(page.addon_state.labels))
        for ((i,j),panel) in page.addon_state.panel_data
            curves=filter(item -> item isa Makie.Lines,panel.axis.scene.plots)
            @test length(curves)==3
            @test Makie.to_color(last(curves).color[])==Makie.to_color(:black)
            for (factor,curve) in zip((2,3,1),curves)
                @test first.(curve[1][])≈f
                @test last.(curve[1][])≈factor.*transform.(z[i,j,:])
            end
        end
    end
    filtered=LineCableModels.plot(artifact,(R,);formulations=2,options...)
    @test first(filtered.addon_state.observed).gridpoint.id.formulation_index==2
    @test first(filtered.addon_state.observed).timings.seconds==2.
    @test all(error -> error.candidate_id.formulation_index==2,first(filtered.addon_state.observed).errors)
    for (all_axis,selected_axis) in zip(first(pages).axes,filtered.axes)
        all_curves=filter(item -> item isa Makie.Lines,all_axis.scene.plots)
        selected=filter(item -> item isa Makie.Lines,selected_axis.scene.plots)
        @test length(selected)==2
        @test selected[1].color[]==all_curves[2].color[]
    end
    blocked=LineCableModels.plot(artifact,(R,);blocks=(1,2),options...)
    @test length(blocked)==2
    @test all(length(page.axes)==2 for page in blocked)
    @test LineCableModels.plot(artifact,(R,);title="My comparison",options...).export_name=="My comparison"
    signed=LineCableModels.plot(artifact;ydata=(G,),options...,controls=true)
    signed.controls[:ylog].active[]=true
    @test all(axis -> axis.yscale[](1e-18)>0 && axis.yscale[](-1e-18)<0,signed.axes)
    @test_throws ArgumentError LineCableModels.plot(artifact,(R,);band=(100.,200.),options...)
    @test_throws ArgumentError LineCableModels.plot(artifact,(R,);band=:wide,options...)
    raw_page=LineCableModels.plot(candidates;ydata=(R,),reference,length_unit=:base,options...)
    @test length(raw_page.axes)==4
    @test length(first(raw_page.addon_state.observed).quantities)==4
    # Rendering and exporting must not consult either poisoned numerical source.
    Z(reference).=NaN;Z(points[1]).=NaN
    rebuilt=report(BenchmarkTableDefinition(),artifact.observed;reference=artifact.reference)
    @test isequal(rebuilt.tables.terms,artifact.tables.terms)
    @test rebuilt.tables.execution.seconds==[1.,2.]
    mktempdir() do directory
        for controls in (false,true)
            page=LineCableModels.plot(artifact,(R,);options...,controls)
            xlims!(first(page.axes),2,70);ylims!(first(page.axes),0,25)
            limits=[axis.finallimits[] for axis in page.axes]
            path=export_svg(page;path=joinpath(directory,"grid-$controls.svg"),open_file=false)
            @test occursin("<svg",read(path,String))
            @test [axis.finallimits[] for axis in page.axes]==limits
        end
    end
end

@testitem "Makie addons / grouping retains independent uncertainty interpretations" tags=[:visual] begin
    using CairoMakie,Measurements
    using LineCableModels.Engine: retain_gridpoint,completed_formulation
    using LineCableModels.Grammar: gridpoint_id,observation_groups
    source_id=gridpoint_id().source_id
    shared=measurement(1.,.1)
    values=(shared,shared,measurement(1.,.1),measurement(1.,.4))
    points=[retain_gridpoint(LineParameters(fill(complex(value,2value),1,1,3),fill(1e-6im,1,1,3),[1.,10.,100.]),
        gridpoint_id(;source_id,formulation_index=index);fields=completed_formulation(Formulation()))
        for (index,value) in enumerate(values)]
    observed=observables(points)
    groups=observation_groups(observed;request=R)
    @test getproperty.(groups,:members)==[[1,2],[3],[4]]
    page=LineCableModels.plot(observed;ydata=(R,),backend=:cairo,display_plot=false,controls=false,open_export=false)
    @test length(page.addon_state.observed)==4
    @test page.addon_state.displayed_indices==[1,3,4]
    @test count(item -> item isa Makie.Lines,only(page.axes).scene.plots)==3
    @test count(item -> item isa Makie.Errorbars,only(page.axes).scene.plots)==3
end

@testitem "Makie addons / quantity assumptions share report and plot grouping" tags=[:visual] begin
    using CairoMakie
    using LineCableModels.Engine: retain_gridpoint,completed_formulation
    using LineCableModels.Grammar: gridpoint_id,observation_groups
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    source_id=gridpoint_id().source_id
    choices=[Formulation(earth_impedance=z,earth_admittance=y) for (z,y) in
        ((:unified,:unified),(:pollaczek1926,:pollaczek1926),(:saad1996,:unified),(:wedepohl1973,:unified),(:xue2018,:xue2018))]
    f=[1.,50.,1e3,1e7]
    base=LineParameters(fill(1.0+2im,1,1,4),fill(3e-6+4e-6im,1,1,4),f)
    points=[retain_gridpoint(base,gridpoint_id(;source_id,formulation_index=index);
        fields=merge(completed_formulation(choice),(inputs=(radius=.01,resistivity=100.),)))
        for (index,choice) in enumerate(choices)]
    artifact=report(BenchmarkTableDefinition(),(reference=LineParameters(Z(base),Y(base),f),candidate=points))
    for (request,count) in ((R,5),(X,5),(G,3),(B,3))
        @test length(observation_groups(artifact.observed;request))==count
        page=LineCableModels.plot(artifact;ydata=(request,),backend=:cairo,display_plot=false,controls=false,open_export=false)
        @test length(page.addon_state.displayed_indices)==count+1
        labels=collect(values(page.addon_state.labels))
        other=request in (R,X) ? "earth Y" : "earth Z"
        @test all(!occursin(other,label) for label in labels)
    end
    for feature in artifact.tables.features
        @test size(feature.relative,1)==(feature.quantity in (:Z,:R,:L,:X) ? 5 : 3)
    end
    @test length(artifact.tables.quantities)==5
    # Composite branches survive captured owner descriptions, including defaults.
    physical=Formulation(internal_impedance=(inner=:default,outer=:default,transfer=:default),
        earth_impedance=(air=:default,earth=:pollaczek1926,mixed=:default),
        earth_admittance=(air=:default,earth=:default,mixed=:default))
    observed=ObservedResult(retain_gridpoint(base,gridpoint_id();fields=completed_formulation(physical)))
    for (request,labels) in ((R,("internal Z(inner)=Schelkunoff","internal Z(outer)=Schelkunoff",
        "internal Z(transfer)=Schelkunoff","earth Z(air)=Unified","earth Z(earth)=Pollaczek","earth Z(mixed)=Unified")),
        (B,("earth Y(air)=Unified","earth Y(earth)=Unified","earth Y(mixed)=Unified")))
        page=LineCableModels.plot(observed;ydata=(request,),backend=:cairo,display_plot=false,controls=false,open_export=false)
        @test all(occursin(label,only(values(page.addon_state.labels))) for label in labels)
    end
end
