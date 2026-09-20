@testitem "ObservedResult / grouping requires physical identity and owned assumptions" tags=[:visual] begin
    using CairoMakie,Measurements
    using LineCableModels.Engine: retain_gridpoint,compare
    using LineCableModels.Grammar: gridpoint_id,observation_groups,observation_labels
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition,tabulate
    source_id=gridpoint_id().source_id
    base=LineParameters(fill(1.0+2im,1,1,3),fill(3.0+4im,1,1,3),[1.,10.,100.])
    function point(problem,formula;rho=100.,control=20.,value=base)
        retain_gridpoint(value,gridpoint_id(;source_id,problem_index=problem,formulation_index=formula);
            fields=(inputs=(radius=.01,earth_resistivity=rho),
                selections=(Z=((:internal,:fixture,(accuracy=control,)),),Y=((:shunt,Symbol("formula",formula),(;)),))))
    end
    first_point=point(1,1)
    irrelevant=point(1,2)
    different_input=point(2,1;rho=200.)
    different_control=point(1,3;control=30.)
    observed=observables([first_point,irrelevant,different_input,different_control])
    groups=observation_groups(observed;request=R)
    @test getproperty.(groups,:members)==[[1,2],[3],[4]]
    @test groups[1].identities==getproperty.([observed[1].gridpoint,observed[2].gridpoint],:id)
    @test length(observation_groups(observed;request=G))==4
    @test occursin("earth_resistivity=200",observation_labels(observed)[3])
    @test all(point -> point.gridpoint.inputs.radius==.01,observed)
    @test length(tabulate(observed))==4
    opts=(backend=:cairo,display_plot=false,controls=false,open_export=false)
    page=LineCableModels.plot(observed;ydata=(R,),opts...)
    @test page.addon_state.displayed_indices==[1,3,4]
    @test only(page.addon_state.display_groups).groups==groups
    @test length(page.addon_state.observed)==4
    @test length(filter(plot -> plot isa Makie.Lines,only(page.axes).scene.plots))==3
    @test_throws ArgumentError LineCableModels.plot(observed;ydata=(L,),opts...)
    @test_throws ArgumentError LineCableModels.plot(observed;clip=false,opts...)
    reference=retain_gridpoint(base,gridpoint_id())
    errors=compare(reference,[first_point,irrelevant,different_input],[R];bands=(:all,(10.,100.)))
    report_points=observables([first_point,irrelevant,different_input];comparisons=errors)
    artifact=report(BenchmarkTableDefinition(),report_points;reference=ObservedResult(reference))
    @test length(artifact.tables.features)==2
    @test all(feature -> length(feature.relative.formula)==1,artifact.tables.features)
    overlay=LineCableModels.plot(artifact;ydata=(R,),problem=2,band=(10.,100.),opts...)
    @test length(overlay.addon_state.observed)==2
    @test first(overlay.addon_state.observed).gridpoint.id.problem_index==2
    curves=filter(plot -> plot isa Makie.Lines,only(overlay.axes).scene.plots)
    @test length(curves)==2
    @test all(curve -> length(curve[1][])==2,curves)
    # Equal marginals do not establish the same dependency interpretation.
    a=measurement(1.,.1);b=measurement(1.,.1)
    uncertain(x)=LineParameters(fill(complex(x,2x),1,1,3),base.Y,[1.,10.,100.])
    independent=observables([point(1,1;value=uncertain(a)),point(1,2;value=uncertain(b))])
    @test length(observation_groups(independent;request=R))==2
    correlated=observables([point(1,1;value=uncertain(a)),point(1,2;value=uncertain(a))])
    @test length(observation_groups(correlated;request=R))==1
end
