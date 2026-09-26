@testitem "ModalAnalysis / structural state and forward imaginary branch" tags=[:unit] begin
    import LineCableModels.Engine as E
    impedance=reshape(ComplexF64[-im],1,1,1)
    phase=LineParameters(impedance,copy(impedance),[50.0])
    modal=compute(ModalAnalysisProblem(phase),ModalAnalysisFormulation(:default))
    @test gamma(modal)[1,1]≈im
    @test_throws DimensionMismatch E.ModalDomain(operators(modal),zeros(ComplexF64,2,1))
    @test_throws ArgumentError E.ModalDomain((Tv=operators(modal).Tv,
        Ti=operators(modal).Ti),zeros(ComplexF64,1,1))
    @test_throws DimensionMismatch ModalOperators(operators(modal).Tv,
        zeros(ComplexF64,2,2,1))
    large=reshape(ComplexF64[1,0,0,1],2,2,1)
    @test_throws DimensionMismatch LineParameters(modal.domain,
        SeriesImpedance(large),ShuntAdmittance(copy(large)),[50.0],ComputationDetails())
    @test_throws ArgumentError PropagationParameters(phase,-1.0,ComputationDetails())
    @test_throws ArgumentError PropagationParameters(modal,-1.0,ComputationDetails())
    @test_throws ArgumentError PropagationParameters(modal,10.0,
        ComputationDetails((segment=(line_length=20.0,),)))
    @test_throws ArgumentError PropagationParameters(modal,10.0,ComputationDetails())
    bound=PropagationParameters(modal;line_length=10.0)
    @test PropagationParameters(bound.parameters,bound.line_length,bound.details) isa
        PropagationParameters
    @test line_length(PropagationParameters(bound;line_length=20.0))==20.0
end

@testitem "ModalAnalysis / derived request acquisition and retained identity" tags=[:unit] begin
    import LineCableModels.Grammar as G
    import LineCableModels.Engine as E
    import LineCableModels.Units as U
    impedance=reshape(ComplexF64[2+im],1,1,1)
    admittance=reshape(ComplexF64[1e-6im],1,1,1)
    base=LineParameters(impedance,admittance,[50.0];
        details=ComputationDetails((inputs=(system=(line_length=10.0,),),)))
    point=details(base).data.gridpoint
    first_phase=E.retain_gridpoint(base,point;
        fields=E.completed_formulation(Formulation()))
    second_phase=E.retain_gridpoint(base,point;
        fields=E.completed_formulation(Formulation(earth_impedance=:unified)))
    distinct_phase=E.retain_gridpoint(base,point;
        fields=E.completed_formulation(Formulation(earth_impedance=:carson1926)))
    first_segment=PropagationParameters(compute(ModalAnalysisProblem(first_phase),
        ModalAnalysisFormulation(:default)))
    second_segment=PropagationParameters(compute(ModalAnalysisProblem(second_phase),
        ModalAnalysisFormulation(:default)))
    @test observe(first_segment,H,abs)≈abs.(H(first_segment))
    @test G.observation_request(first_segment,(H,abs)).quantity==U.quantity(H,abs)
    @test observe(first_segment,gamma,real)==real.(gamma(first_segment))
    @test observe(first_segment,H,imag)==imag.(H(first_segment))
    requests=((H,abs),(gamma,real),(H,imag))
    first_observed=ObservedResult(first_segment,requests;complete_pairs=true)
    second_observed=ObservedResult(second_segment,requests;complete_pairs=true)
    @test G.observation_quantity(first_phase,G.observation_requests(first_phase,(LineCableModels.Engine.Z,)).retained[1]).assumptions==
        G.observation_quantity(second_phase,G.observation_requests(second_phase,(LineCableModels.Engine.Z,)).retained[1]).assumptions
    @test length(G.observation_groups([first_observed,second_observed];request=(H,abs)))==1
    distinct_observed=ObservedResult(PropagationParameters(compute(
        ModalAnalysisProblem(distinct_phase),ModalAnalysisFormulation(:default))),requests;
        complete_pairs=true)
    @test length(G.observation_groups([first_observed,distinct_observed];request=(H,abs)))==2
    @test length(G.observation_labels([first_observed];request=(H,abs)))==1
    alternate_modal=compute(
        ModalAnalysisProblem(first_phase),ModalAnalysisFormulation(:default;
            options=(iteration=(max_iterations=101,),)))
    alternate=ObservedResult(PropagationParameters(alternate_modal),requests;
        complete_pairs=true)
    labels=G.observation_labels([first_observed,alternate];request=(H,abs))
    @test labels[1]!=labels[2]
    @test any(label -> occursin("max_iterations",label),labels)
    primary=ObservedResult(first_segment.parameters,(E.R,E.X))
    primary_alternate=ObservedResult(alternate_modal,(E.R,E.X))
    @test primary.quantities[1].assumptions!=primary_alternate.quantities[1].assumptions
    @test length(G.observation_groups([primary,primary_alternate];request=E.R))==2
    primary_labels=G.observation_labels([primary,primary_alternate];request=E.R)
    @test primary_labels[1]!=primary_labels[2]
    @test any(label -> occursin("max_iterations",label),primary_labels)
    descriptions=details(first_segment.parameters).data.formulation_fields
    @test descriptions.all==descriptions.Z==descriptions.Y
    @test descriptions.all !== descriptions.Z && descriptions.all !== descriptions.Y &&
        descriptions.Z !== descriptions.Y
    candidate_source=G.gridpoint_id().source_id
    candidates=[E.retain_gridpoint(result,G.gridpoint_id(;source_id=candidate_source,
        formulation_index=index)) for (index,result) in enumerate((first_segment.parameters,alternate_modal))]
    reference=E.retain_gridpoint(first_segment.parameters,G.gridpoint_id())
    comparisons=E.compare(reference,candidates,[E.R];bands=(:all,))
    @test isequal(comparisons[1].absolute,comparisons[2].absolute)
    @test isequal(comparisons[1].relative,comparisons[2].relative)
    @test comparisons[1].assumptions!=comparisons[2].assumptions
    compared=observables(candidates,(E.R,E.X);comparisons)
    @test length(G.observation_groups(compared;request=E.R,band=:all,
        normalization=:reference_rms,reference=details(reference).data.gridpoint))==2
    artifact=report(LineCableModels.ReportBuilder.BenchmarkTableDefinition((E.R,);bands=(:all,)),
        compared;reference=ObservedResult(reference,(E.R,E.X)))
    @test size(only(artifact.tables.features).absolute,1)==2
    admittance_cases=[E.retain_gridpoint(base,point;
        fields=E.completed_formulation(Formulation(earth_admittance=selection)))
        for selection in (:unified,:ideal)]
    phase_admittance=[ObservedResult(phase,(E.R,E.X)) for phase in admittance_cases]
    @test only(unique(G.observation_labels(phase_admittance;request=E.R))) isa String
    modal_admittance=[ObservedResult(compute(ModalAnalysisProblem(phase),
        ModalAnalysisFormulation(:default)),(E.R,E.X)) for phase in admittance_cases]
    @test length(G.observation_groups(modal_admittance;request=E.R))==2
    admittance_labels=G.observation_labels(modal_admittance;request=E.R)
    @test admittance_labels[1]!=admittance_labels[2]
    @test any(label -> occursin("Ideal",label),admittance_labels)
    modal_impedance=[ObservedResult(compute(ModalAnalysisProblem(phase),
        ModalAnalysisFormulation(:default)),(E.G,E.B)) for phase in (second_phase,distinct_phase)]
    @test length(G.observation_groups(modal_impedance;request=E.G))==2
    impedance_labels=G.observation_labels(modal_impedance;request=E.G)
    @test impedance_labels[1]!=impedance_labels[2]
    phase_Zc=Base.Fix2(Zc,(domain=PhaseDomain,))
    target=only(G.unit_targets(((phase_Zc,abs),),:pul;
        overrides=Dict(Zc=>:kilo)))
    @test target==U.units(:kilo,:ohm)
    mktempdir() do folder
        path=joinpath(folder,"modal-derived.json")
        LineCableModels.save(first_observed,path)
        restored=LineCableModels.import_data(Val(:observed),path)
        @test observe(restored,(gamma,real))≈observe(first_observed,(gamma,real))
        @test observe(restored,(H,imag))≈observe(first_observed,(H,imag))
    end
end

@testitem "ModalAnalysis / custom passive controls survive retained labels" tags=[:unit] begin
    import LineCableModels.ModalAnalysis as MA
    import LineCableModels.Engine as E
    import LineCableModels.Grammar as G

    struct ScaleModal <: AbstractFormulation
        scale::Float64
    end
    LineCableModels.formula_id(::ScaleModal)=:scale_modal
    Base.NamedTuple(value::ScaleModal)=(identifier=:scale_modal,
        parameters=(scale=value.scale,),options=(;))
    LineCableModels.formulation_options(::ScaleModal)=FormulationOptions()
    LineCableModels.description(value::ScaleModal;compact::Bool=false)=
        "Basis scale $(value.scale)"
    LineCableModels.description(::Type{ScaleModal},::Val{:scale},value::Real;
        compact::Bool=false)="scale=$(value)"
    E.initialize_buffers(::ScaleModal,::Type,input,invariants,buffers)=buffers
    function MA.decompose!(selected::ScaleModal,workspace,parameters,options)
        workspace.Tv[1,1,1]=selected.scale
        workspace.Ti[1,1,1]=1
        workspace.roots[1,1]=sqrt(workspace.input.Z[1,1,1]*workspace.input.Y[1,1,1])
        return workspace
    end

    z=reshape(ComplexF64[2+im],1,1,1)
    y=reshape(ComplexF64[1e-6im],1,1,1)
    base=LineParameters(z,y,[50.0])
    phase=E.retain_gridpoint(base,details(base).data.gridpoint;
        fields=E.completed_formulation(Formulation()))
    first_modal=compute(ModalAnalysisProblem(phase),ModalAnalysisFormulation(ScaleModal(1.0)))
    second_modal=compute(ModalAnalysisProblem(phase),ModalAnalysisFormulation(ScaleModal(2.0)))
    @test Tv(first_modal)!=Tv(second_modal)
    first=ObservedResult(first_modal,((Tv,abs),E.R,E.X);complete_pairs=true)
    second=ObservedResult(second_modal,((Tv,abs),E.R,E.X);complete_pairs=true)
    @test length(G.observation_groups([first,second];request=(Tv,abs)))==2
    @test first.quantities[1].assumptions!=second.quantities[1].assumptions
    fields=details(first_modal).data.formulation_fields.all
    modal_field=only(filter(field -> field.meaning==(:transformation,),fields))
    @test modal_field.value=="Basis scale 1.0"
    @test any(control -> control.scope==(:parameters,:scale),modal_field.control_fields)
    for request in ((Tv,abs),E.R)
        labels=G.observation_labels([first,second];request)
        @test labels[1]!=labels[2]
        @test all(label -> occursin("scale=",label),labels)
    end
end
