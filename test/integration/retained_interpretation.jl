@testitem "Makie / retained units coordinates and ragged frequencies" tags=[:visual] begin
    using CairoMakie
    z=reshape(complex.(1.:12.,101.:112.),2,2,3)
    a=ObservedResult(LineParameters(z,2z,[1.,10.,100.]);length_unit=:base)
    b=ObservedResult(LineParameters(z,2z,[1.,10.,100.]),
        ((R,[2,1],[2,1],[1,3]),(X,[2,1],[2,1],[1,3]),
         (G,[2,1],[2,1],[2,3]),(B,[2,1],[2,1],[2,3]));
        length_unit=:kilo,frequency_unit=:kilo)
    options=(backend=:cairo,display_plot=false,controls=false,open_export=false)
    pages=LineCableModels.plot([a,b];ydata=(R,G),options...)
    @test length(pages)==2
    for (page,selector,scale,samples) in zip(pages,(R,G),(1.,2.),([1,3],[2,3]))
        @test length(page.axes)==4
        for (axis,(i,j)) in zip(page.axes,((1,1),(1,2),(2,1),(2,2)))
            curves=filter(p -> p isa Makie.Lines,axis.scene.plots)
            @test length(curves)==2
            expected1=Float32.(real.(z[i,j,:]).*scale)
            expected2=expected1[samples]
            @test getindex.(curves[1][1][],2)≈expected1
            @test getindex.(curves[2][1][],2)≈expected2
            @test getindex.(curves[1][1][],1)≈[1.,10.,100.]
            @test getindex.(curves[2][1][],1)≈[1.,10.,100.][samples]
        end
    end
    @test b.quantities[1].coordinates.rows==[2,1]
    @test b.quantities[1].coordinates.frequency_unit==LineCableModels.Units.units(:kilo,:hertz)
end

@testitem "Makie / quantity capacity and original modal panel identities" tags=[:visual] begin
    using CairoMakie,LinearAlgebra
    options=(backend=:cairo,display_plot=false,controls=false,open_export=false)
    z=fill(0.0+0im,3,3,2)
    for i in 1:3
        z[i,i,:].=i+2im
    end
    raw=LineParameters(z,z.*1e-6,[1.,10.])
    retained=ObservedResult(raw,(R,L);clip=false,length_unit=:base)
    # The source owner declares modal coordinates; the renderer must not infer diag.
    modal=ObservedResult(retained.gridpoint,
        [merge(q,(coordinates=merge(q.coordinates,(domain=:ModalDomain,)),)) for q in retained.quantities],
        retained.errors,retained.timings)
    pages=LineCableModels.plot(modal;ydata=(R,L),layout=(2,2),options...)
    @test length(pages)==8
    @test [p.addon_state.panel_page.dimensions for p in pages]==repeat([(2,2),(2,1),(1,2),(1,1)],2)
    @test sum(length(p.axes) for p in pages)==18
    @test Set(k for p in pages[1:4] for k in keys(p.addon_state.panel_data))==Set(Iterators.product(1:3,1:3))
    @test length(LineCableModels.plot(modal;ydata=(R,L),layout=(1,1),options...))==18
    @test length(LineCableModels.plot(modal;ydata=(R,L),layout=(1,2),options...))==12
    @test_throws ArgumentError LineCableModels.plot(modal;blocks=(2,2),options...)
    residual=ObservedResult(modal.gridpoint,[merge(q,(values=q.values .+ 1e-12,)) for q in modal.quantities],[],(;))
    full=LineCableModels.plot(residual;ydata=(R,),options...)
    @test length(full.axes)==9
    @test last(only(filter(p -> p isa Makie.Lines,full.addon_state.panel_data[(1,2)].axis.scene.plots))[1][][1])≈1e-12
    diagonal=ObservedResult(raw,((R,diag,:,:),(L,diag,:,:));clip=false)
    flow=LineCableModels.plot(diagonal;ydata=((R,diag),),layout=(1,2),options...)
    @test [length(p.axes) for p in flow]==[2,1]
    @test collect(keys(last(flow).addon_state.panel_data))==[(3,3)]
    # An owner can retain a rectangular matrix even though line solvers are square.
    q=first(modal.quantities)
    rect=merge(q,(values=q.values[:,1:2,:],available=true,engineering_zero=false,missing_reason=nothing,
        coordinates=merge(q.coordinates,(extent=(3,2,2),columns=[1,2])),))
    rectangular=ObservedResult((id=nothing,),[rect],[],(;))
    @test length(LineCableModels.plot(rectangular;options...).axes)==6
    one=LineParameters(reshape([1+2im],1,1,1),reshape([1+2im],1,1,1),[1.])
    single=LineCableModels.plot(one;ydata=(R,),options...)
    @test length(single.axes)==1
    @test any(p -> p isa Makie.Scatter && length(p[1][])==1,only(single.axes).scene.plots)
end

@testitem "Makie / retained assembly unions and element-index points" tags=[:visual] begin
    using CairoMakie,Measurements
    options=(backend=:cairo,display_plot=false,controls=true,open_export=false)
    a=CableConstants([:core,:sheath],[1.,2.],[.1,.2],[.01,.02],[.001,.002],50.)
    b=CableConstants([:core,:screen],[3.,4.],[.3,.4],[.03,.04],[.003,.004],60.)
    observed=[ObservedResult(a;length_unit=:base),ObservedResult(b;length_unit=:base)]
    page=LineCableModels.plot(observed;ydata=(R,),options...)
    axis=only(page.axes)
    @test axis.dim1_conversion[] isa Makie.CategoricalConversion
    @test !haskey(page.controls,:xlog)
    @test isempty(filter(p -> p isa Makie.Lines,axis.scene.plots))
    points=filter(p -> p isa Makie.Scatter,axis.scene.plots)
    @test first.(points[1][1][])==[1.,2.]
    @test first.(points[2][1][])==[1.,3.]
    @test last.(points[2][1][])==[3.,4.]
    @test last.(axis.dim1_conversion[].int_to_category)==["core","sheath","screen"]
    raw=LineCableModels.plot((a,b);ydata=(R,),length_unit=:base,options...)
    @test length(raw.axes)==1
    template=first(first(observed).quantities)
    function indexed(values,indices=())
        q=merge(template,(request=R,values,coordinates=(kind=:array,indices,extent=size(values)),
            available=true,engineering_zero=false,missing_reason=nothing,thresholds=nothing,clipped=false))
        ObservedResult((id=nothing,),[q],[],(;))
    end
    vector=LineCableModels.plot(indexed([-1.,-2.,-3.],([4,7,9],));options...)
    @test first.(only(filter(p -> p isa Makie.Scatter,only(vector.axes).scene.plots))[1][])==[4.,7.,9.]
    @test only(vector.axes).xlabel[]=="Original element index"
    scalar=LineCableModels.plot(indexed(2. ± .5);options...)
    @test only(scalar.axes).xlabel[]=="Retained element position"
    @test length(only(filter(p -> p isa Makie.Scatter,only(scalar.axes).scene.plots))[1][])==1
    @test_throws ArgumentError LineCableModels.plot(indexed(ones(2,2));options...)
end

@testitem "ObservedResult / retained band associations and sample intersections" begin
    using LineCableModels.Engine: retain_gridpoint,compare
    using LineCableModels.Grammar: gridpoint_id,observation_product
    base=LineParameters(reshape(complex.([1.,2.,3.],[2.,3.,4.]),1,1,3),fill(1+2im,1,1,3),[1.,10.,100.])
    reference=retain_gridpoint(base,gridpoint_id())
    candidates=[retain_gridpoint(base,gridpoint_id()) for _ in 1:2]
    errors=compare(reference,candidates,[R];bands=((10.,100.),))
    a=ObservedResult(candidates[1],((R,:,:,[1,3]),(X,:,:,[1,3]));comparisons=errors[1:1])
    b=ObservedResult(candidates[2],((R,:,:,[1,2]),(X,:,:,[1,2]));comparisons=errors[2:2])
    ref=ObservedResult(reference,((R,:,:,[2,3]),(X,:,:,[2,3])))
    products=observation_product((a,b,ref),R;band=(10.,100.),reference_id=ref.gridpoint.id)
    @test [p.coordinates.samples for p in products]==[[3],[2],[2,3]]
    @test [p.coordinates.frequencies for p in products]==[[100.],[10.],[10.,100.]]
    @test a.quantities[1].coordinates.samples==[1,3]
    @test_throws ArgumentError observation_product((a,b,ref),R;band=:absent)
    conflicting=merge(only(b.errors),(settings=merge(only(b.errors).settings,(indices=[1,3],)),))
    invalid=ObservedResult(b.gridpoint,b.quantities,[conflicting],b.timings)
    @test_throws ArgumentError observation_product((a,invalid,ref),R;band=(10.,100.))
    empty=ObservedResult(a,((R,:,:,1:1),))
    @test isempty(first(observation_product((empty,b,ref),R;band=(10.,100.))).coordinates.samples)
    @test_throws ArgumentError observation_product((empty,),R;band=(10.,100.))
    other_reference=retain_gridpoint(base,gridpoint_id())
    additional=compare(other_reference,candidates[1],[R];bands=((10.,100.),))
    ambiguous=ObservedResult(a.gridpoint,a.quantities,[a.errors;additional],a.timings)
    @test_throws ArgumentError observation_product((ambiguous,),R;band=(10.,100.))
    @test only(observation_product((ambiguous,),R;band=(10.,100.),reference_id=ref.gridpoint.id)).coordinates.samples==[3]
end

@testitem "Makie / observed public forwarding preserves native array dispatch" tags=[:visual] begin
    using CairoMakie
    extension=Base.get_extension(LineCableModels,:LineCableModelsMakieExt)
    @test LineCableModels.plot===LineCableModels.PlotBuilder.plot
    for T in (ObservedResult,Vector{ObservedResult},Tuple{ObservedResult,ObservedResult},Vector{LineParameters})
        @test which(Makie.plot,Tuple{T}).module===extension
    end
    @test which(Makie.plot,Tuple{Vector{Float64}}).module!==extension
    source=LineParameters(fill(1.0+2im,1,1,3),fill(1e-6+2e-6im,1,1,3),[1.,10.,100.])
    observed=ObservedResult(source,(R,X);clip=false)
    for input in (observed,[observed],(observed,observed))
        p=Makie.plot(input;ydata=(R,),backend=:cairo,display_plot=false,controls=false)
        @test p isa UIPlot
        @test length(p.axes)==1
    end
end

@testitem "Makie / public conveniences re-express retained display units" tags=[:visual] begin
    using CairoMakie, Measurements
    using LineCableModels.Engine: retain_gridpoint,completed_formulation,compare
    using LineCableModels.Grammar: gridpoint_id,observation_product
    using LineCableModels.ReportBuilder: ReportArtifact,BenchmarkTableDefinition
    U=LineCableModels.Units
    options=(backend=:cairo,display_plot=false,controls=false,open_export=false)
    source_id=gridpoint_id().source_id
    z=reshape(complex.(1.:12.,21.:32.),2,2,3)
    f=[10.,100.,1000.]
    point(index;values=z,frequency=f)=retain_gridpoint(LineParameters(copy(values),values.*1e-6,frequency),
        gridpoint_id(;source_id,problem_index=index,formulation_index=index);fields=completed_formulation(Formulation()))
    a,b=point(1),point(2;values=2z)
    reference=point(3;values=.5z)
    request=@observe R[:,:,:]
    curves(p)=filter(item -> item isa Makie.Lines,first(p.axes).scene.plots)
    data(p)=[copy(c[1][]) for c in curves(p)]
    raw=LineCableModels.plot(a;ydata=(request,),length_unit=:base,options...)
    retained=ObservedResult(a,(request,);complete_pairs=true,length_unit=:base)
    explicit=LineCableModels.plot(retained;ydata=(request,),options...)
    @test data(raw)==data(explicit)
    @test keys(raw.addon_state.panel_data)==keys(explicit.addon_state.panel_data)
    for inputs in ([a,b],(a,b),ParametricResult(nothing,[a,b],(problems=[:unused],formulations=[:a,:b]),ComputationDetails()))
        p=LineCableModels.plot(inputs;ydata=(request,),length_unit=:base,reference,options...)
        @test last.(curves(p)[1][1][])≈real.(z[1,1,:])
        @test last.(curves(p)[2][1][])≈2real.(z[1,1,:])
        @test last.(curves(p)[3][1][])≈.5real.(z[1,1,:])
        @test all(o -> observation_product(o,request).unit==U.units(:base,:ohm;per=(:base,:meter)),p.addon_state.observed)
    end
    positional=LineCableModels.plot(a,b,(request,);length_unit=:base,options...)
    @test length(curves(positional))==2
    native=Makie.plot(a;ydata=(request,),length_unit=:base,options...)
    @test data(native)==data(raw)
    named=LineCableModels.plot((first_case=a,second_case=b);ydata=(request,),reference,length_unit=:base,options...)
    @test Set(values(named.addon_state.labels))==Set(("first_case","second_case",description(Formulation().methods.earth_impedance;compact=true)*" (reference)"))
    rich_labels=(rich("First",font=:bold),rich("Second",font=:italic))
    rich_page=LineCableModels.plot((a,b);ydata=(request,),reference,series_labels=rich_labels,options...)
    @test rich_page.addon_state.labels[:result_1]===rich_labels[1]
    @test endswith(rich_page.addon_state.labels[:result_3]," (reference)")
    kilo=ObservedResult(a;length_unit=:kilo)
    original=deepcopy(kilo)
    for (units,scale) in ((nothing,1000.),(:kilo,1000.),(:base,1.))
        p=units===nothing ? LineCableModels.plot(kilo;ydata=(R,),options...) :
            LineCableModels.plot(kilo;ydata=(R,),length_unit=units,options...)
        @test last.(only(curves(p))[1][])≈scale.*real.(z[1,1,:])
        repeated=LineCableModels.plot(first(p.addon_state.observed);ydata=(R,),length_unit=:base,options...)
        @test data(repeated)==data(raw)
    end
    single=ObservedResult(kilo,(R,);length_unit=:base)
    @test length(single.quantities)==1
    @test data(LineCableModels.plot(single;ydata=(R,),length_unit=:base,options...))==data(raw)
    freq=LineCableModels.plot(single;ydata=(R,),freq_unit=:kilo,options...)
    @test first.(only(curves(freq))[1][])≈f./1000
    @test first.(only(curves(LineCableModels.plot(first(freq.addon_state.observed);frequency_unit=:base,options...)))[1][])≈f
    @test data(LineCableModels.plot(kilo;ydata=(R,),units=(U.units(:base,:ohm;per=(:base,:meter)),),options...))==data(raw)
    @test last.(only(curves(LineCableModels.plot(single;ydata=(R,),quantity_units=:milli,options...)))[1][])≈1000real.(z[1,1,:])
    mixed=LineCableModels.plot((kilo,b);ydata=(R,),length_unit=:base,options...)
    @test last.(curves(mixed)[1][1][])≈real.(z[1,1,:])
    mixed_omitted=LineCableModels.plot((kilo,b);ydata=(R,),options...)
    @test last.(curves(mixed_omitted)[1][1][])≈1000real.(z[1,1,:])

    comparisons=compare(reference,a,[R];bands=(:all,))
    observed=ObservedResult(a;comparisons,length_unit=:kilo,timings=(seconds=.1,))
    ref=ObservedResult(reference;length_unit=:kilo)
    artifact=report(BenchmarkTableDefinition(),[observed];reference=ref)
    saved_tables=deepcopy(artifact.tables)
    # Reopening these raw sources after detachment would return invalid data.
    Z(a).=NaN;Z(reference).=NaN
    for candidates in (observed,[observed])
        report_handle=ReportArtifact(candidates,ref,artifact.tables,nothing,nothing)
        p=LineCableModels.plot(report_handle;ydata=(R,),length_unit=:base,options...)
        @test last.(curves(p)[1][1][])≈real.(z[1,1,:])
        @test last.(curves(p)[2][1][])≈.5real.(z[1,1,:])
        @test first(p.addon_state.observed).errors==observed.errors
        @test first(p.addon_state.observed).timings==observed.timings
        axisscale!(p,:y,:log10);resetview!(p)
        mktempdir() do directory
            @test isfile(export_svg(p;path=joinpath(directory,"retained.svg"),open_file=false))
        end
        @test length(curves(LineCableModels.plot(report_handle;ydata=(R,),reference=nothing,length_unit=:base,options...)))==1
        @test length(curves(LineCableModels.plot(report_handle;ydata=(R,),reference=single,length_unit=:base,options...)))==2
    end
    @test isequal(artifact.tables,saved_tables)
    @test isequal(kilo.quantities,original.quantities)
    @test kilo.gridpoint==original.gridpoint
    for kwargs in ((clip=false,),(atol=0.,),(frequencies=[1.,2.,3.],),
            (freq_unit=:kilo,frequency_unit=:base,), (units=(U.units(:base,:farad),),))
        @test_throws ArgumentError LineCableModels.plot(single;ydata=(R,),kwargs...,options...)
    end
    @test_throws ArgumentError LineCableModels.plot(single;ydata=(X,),options...)
    @test_throws ArgumentError LineCableModels.plot(single,(R,);ydata=(R,),options...)
    @test_throws ArgumentError LineCableModels.plot(b;ydata=(R,),freq_unit=:base,frequency_unit=:kilo,options...)
end
