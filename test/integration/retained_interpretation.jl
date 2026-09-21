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
