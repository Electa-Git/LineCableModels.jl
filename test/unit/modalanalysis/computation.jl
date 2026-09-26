@testitem "ModalAnalysis / retained directions, finite responses, transport, and diagnostics" tags=[:unit] begin
    using Test
    using LinearAlgebra
    using Serialization
    import LineCableModels.ModalAnalysis as MA
    import LineCableModels.Engine as E
    import LineCableModels.Grammar as G

    struct FixedModal <: AbstractFormulation
        voltage::Array{ComplexF64,3}
        current::Array{ComplexF64,3}
        roots::Matrix{ComplexF64}
        allocations::Base.RefValue{Int}
        calculations::Base.RefValue{Int}
        common::Base.RefValue{Any}
    end
    FixedModal(voltage,current,roots,allocations,calculations) =
        FixedModal(voltage,current,roots,allocations,calculations,Ref{Any}(nothing))
    LineCableModels.formula_id(::FixedModal)=:fixed_modal
    Base.NamedTuple(::FixedModal)=(identifier=:fixed_modal,)
    G.formulation_options(::FixedModal)=FormulationOptions()
    function E.initialize_buffers(selected::FixedModal,::Type{T},input,invariants,buffers) where {T}
        selected.allocations[]+=1
        @test size(buffers.admittance_impedance_product)==(invariants.n,invariants.n)
        @test eltype(buffers.admittance_impedance_product)==T
        selected.common[]=buffers
        return merge(buffers,(fixed=zeros(ComplexF64,invariants.n,invariants.n),))
    end
    function MA.decompose!(selected::FixedModal,workspace,parameters,options)
        selected.calculations[]+=1
        @test all(key -> getproperty(workspace.buffers,key) ===
            getproperty(selected.common[],key),keys(selected.common[]))
        copyto!(workspace.Tv,selected.voltage)
        copyto!(workspace.Ti,selected.current)
        copyto!(workspace.roots,selected.roots)
        return workspace
    end

    Tv0=ComplexF64[1+0.2im 0.3-0.4im;0.2+0.1im 1.4-0.2im]
    Ti0=ComplexF64[1.2+0.3im 0.2+0.1im;0.4-0.2im 1.1+0.5im]
    z=ComplexF64[2+3im,4+1im]
    y=ComplexF64[0.1+0.2im,0.2+0.3im]
    roots=sqrt.(z.*y)
    Zp=Tv0*Diagonal(z)/Ti0
    Yp=Ti0*Diagonal(y)/Tv0
    details0=ComputationDetails((inputs=(system=(line_length=600.0,),),))
    phase=LineParameters(PhaseDomain,SeriesImpedance(reshape(Zp,2,2,1)),
        ShuntAdmittance(reshape(Yp,2,2,1)),[50.0],details0)
    selected=FixedModal(reshape(Tv0,2,2,1),reshape(Ti0,2,2,1),
        reshape(roots,2,1),Ref(0),Ref(0))
    mf=ModalAnalysisFormulation(selected)
    modal=compute(ModalAnalysisProblem(phase),mf)
    @test selected.allocations[]==1
    @test selected.calculations[]==1
    @test Tv(modal)[:,:,1]==Tv0
    @test Ti(modal)[:,:,1]==Ti0
    @test gamma(modal)[:,1]==roots
    @test modal.Z.values[:,:,1]≈Diagonal(z)
    @test modal.Y.values[:,:,1]≈Diagonal(y)
    @test Z(transform(PhaseDomain,modal))[:,:,1]≈Zp
    @test Y(transform(PhaseDomain,modal))[:,:,1]≈Yp
    @test Zc(modal)[:,1]≈z./roots
    @test Yc(modal)[:,1]≈y./roots
    @test Zc(modal,PhaseDomain)[:,:,1]≈Tv0*Diagonal(z./roots)/Ti0
    @test Yc(modal,PhaseDomain)[:,:,1]≈Ti0*Diagonal(y./roots)/Tv0

    segment=PropagationParameters(modal)
    @test line_length(segment)==600.0
    short=PropagationParameters(segment;line_length=150.0)
    zero_segment=PropagationParameters(segment;line_length=0.0)
    @test H(short)[:,1]≈exp.(-roots.*150.0)
    @test H(zero_segment)==ones(ComplexF64,2,1)
    @test H(zero_segment,PhaseDomain;field=:voltage)[:,:,1]≈Matrix{ComplexF64}(I,2,2)
    @test H(short,PhaseDomain;field=:voltage)[:,:,1]≈Tv0*Diagonal(H(short)[:,1])/Tv0
    @test H(short,PhaseDomain;field=:current)[:,:,1]≈Ti0*Diagonal(H(short)[:,1])/Ti0
    @test H(short,PhaseDomain;field=:voltage)[:,:,1]*Zc(short,PhaseDomain)[:,:,1]≈
        Zc(short,PhaseDomain)[:,:,1]*H(short,PhaseDomain;field=:current)[:,:,1]
    zero_roots=LineParameters(E.ModalDomain(operators(modal),zeros(ComplexF64,size(gamma(modal)))),
        modal.Z,modal.Y,modal.f,modal.details)
    first_zero=PropagationParameters(zero_roots;line_length=100.0)
    second_zero=PropagationParameters(first_zero;line_length=200.0)
    @test H(first_zero)==H(second_zero)
    @test details(first_zero).data.gridpoint!=details(second_zero).data.gridpoint
    @test selected.calculations[]==1

    total=LineParameters(PhaseDomain,SeriesImpedance(reshape(Zp.*600,2,2,1);basis=:total),
        ShuntAdmittance(reshape(Yp.*600,2,2,1);basis=:total),[50.0],details0)
    total_selected=FixedModal(selected.voltage,selected.current,selected.roots.*600,
        Ref(0),Ref(0))
    total_modal=compute(ModalAnalysisProblem(total),ModalAnalysisFormulation(total_selected))
    @test gamma(total_modal)[:,1]≈roots.*600
    @test Zc(total_modal)≈Zc(modal)
    @test gamma(PropagationParameters(total_modal))≈gamma(segment)
    @test H(PropagationParameters(total_modal))≈H(segment)

    phase_results=ParametricResult(Combinatorial(mf),[phase,phase])
    problems=Gridspace{ModalAnalysisProblem}(phase_results)
    @test length(problems)==2
    run=compute(problems,mf)
    @test length(run)==2
    @test eltype(run)<:LineParameters
    @test size.(gamma.(run))==[(2,1),(2,1)]
    @test !Base.mightalias(Tv(modal),Tv(run[1]))
    @test !Base.mightalias(Ti(modal),Ti(run[1]))
    @test !Base.mightalias(gamma(modal),gamma(run[1]))
    @test !Base.mightalias(Z(modal),Z(run[1]))
    @test !Base.mightalias(details(modal).data.modal.diagnostics.z_coupling,
        details(run[1]).data.modal.diagnostics.z_coupling)
    @test length(Tv.(run))==2
    @test length(Ti.(run))==2
    @test Zc.(run,Ref(PhaseDomain))[1]==Zc(run[1],PhaseDomain)
    @test length(PropagationParameters.(run))==2
    @test length(collect(Gridspace{PropagationParameters}(run)))==2
    segments=collect(Gridspace{PropagationParameters}(run))
    finite_observed=ObservedResult(segments[1],(H,))
    @test finite_observed.gridpoint.source_gridpoint==
        details(segments[1]).data.source_gridpoint
    @test finite_observed.gridpoint.source_gridpoint==
        details(run[1]).data.gridpoint
    @test finite_observed.gridpoint.source_gridpoint!=
        details(phase).data.gridpoint
    rebound=PropagationParameters(segments[1];line_length=75.0)
    rebound_observed=ObservedResult(rebound,(H,))
    @test rebound_observed.gridpoint.source_gridpoint==
        details(segments[1]).data.gridpoint
    @test H.(segments)[1]==H(segments[1])
    @test H.(segments,Ref(PhaseDomain);field=:current)[2]==
        H(segments[2],PhaseDomain;field=:current)
    segment_lengths=PropagationParameters(modal;line_length=Grid((100.0,200.0)))
    @test segment_lengths isa Gridspace{PropagationParameters}
    @test line_length.(collect(segment_lengths))==[100.0,200.0]
    zipped=Gridspace{PropagationParameters}((parameters,length) ->
        PropagationParameters(parameters;line_length=length),
        (run,Grid((100.0,200.0)));combine=:zip)
    @test line_length.(collect(zipped))==[100.0,200.0]
    @test selected.allocations[]==3
    @test selected.calculations[]==3
    variants=ModalAnalysisFormulation(Grid((selected,selected)))
    product=compute(problems,variants)
    @test length(product)==4
    @test product[1,1]===product[1]
    @test product[2,2]===product[4]
    @test all(value -> value isa LineParameters,product)

    bad=FixedModal(copy(selected.voltage),copy(selected.current),copy(selected.roots),Ref(0),Ref(0))
    bad.voltage[1,2,1]+=0.25
    bad_form=ModalAnalysisFormulation(bad)
    noisy=@test_logs (:warn,r"off-diagonal tolerance") compute(
        ModalAnalysisProblem(phase),bad_form;
        options=(offdiagonal_tolerance=1e-12,))
    @test noisy.Z.values[1,2,1]!=0
    @test !isempty(details(noisy).data.modal.diagnostics.z_coupling)
    singular=FixedModal(zeros(ComplexF64,2,2,1),selected.current,selected.roots,Ref(0),Ref(0))
    @test_throws SingularException compute(ModalAnalysisProblem(phase),
        ModalAnalysisFormulation(singular))

    io=IOBuffer()
    Serialization.serialize(io,(modal,short))
    seekstart(io)
    loaded_modal,loaded_segment=Serialization.deserialize(io)
    @test Z(loaded_modal)==Z(modal)
    @test H(loaded_segment)==H(short)
    @test Z(transform(PhaseDomain,loaded_modal))[:,:,1]≈Zp
    @test selected.calculations[]==7
    no_length=LineParameters(reshape(Zp,2,2,1),reshape(Yp,2,2,1),[50.0])
    unbound=compute(ModalAnalysisProblem(no_length),mf)
    @test line_length(unbound)===nothing
    @test_throws ArgumentError PropagationParameters(unbound)
    @test line_length(PropagationParameters(unbound;line_length=25.0))==25.0
    @test_throws ArgumentError LineCableModels.ImportExport.serialize_value(modal)

    timed_phase=E.retain_gridpoint(phase,details(phase).data.gridpoint;
        fields=(timing=(wall_seconds=1.0,),))
    untimed=compute(ModalAnalysisProblem(timed_phase),mf)
    @test !haskey(details(untimed).data,:timing)
    modal_timed=compute(ModalAnalysisProblem(phase),mf;options=(timing=true,))
    @test haskey(details(modal_timed).data,:timing)
    @test details(modal_timed).data.gridpoint==details(modal).data.gridpoint
    @test details(modal_timed).data.source_gridpoint==details(modal).data.source_gridpoint
    @test !haskey(details(PropagationParameters(modal_timed)).data,:timing)

    voltage_H=Base.Fix2(H,(domain=PhaseDomain,field=:voltage))
    current_H=Base.Fix2(H,(domain=PhaseDomain,field=:current))
    observed=ObservedResult(short,(gamma,Zc,Yc,H,Tv,Ti,voltage_H,current_H,(H,abs));
        complete_pairs=true)
    @test observed.gridpoint.source_gridpoint==details(short).data.source_gridpoint
    @test length(observed.quantities)==18
    @test G.observation_product(observed,(Tv,real)).coordinates.column_labels==["1","2"]
    artifact=report(TableReportDefinition((voltage_H,current_H,Tv,Ti)),observed)
    @test keys(artifact.tables.H)==(:H_phase_voltage_real,:H_phase_voltage_imag,
        :H_phase_current_real,:H_phase_current_imag)
    @test size(artifact[(Tv,real)])==(1,5)
    mktempdir() do folder
        for extension in (".json",".jls")
            saved=joinpath(folder,"modal-observations"*extension)
            LineCableModels.save(observed,saved)
            restored=LineCableModels.import_data(Val(:observed),saved)
            @test length(restored.quantities)==18
            @test observe(restored,(voltage_H,real))≈observe(observed,(voltage_H,real))
            @test observe(restored,(current_H,real))≈observe(observed,(current_H,real))
            @test observe(restored,(H,abs))≈observe(observed,(H,abs))
            @test restored.gridpoint.id==observed.gridpoint.id
            @test restored.gridpoint.source_gridpoint==observed.gridpoint.source_gridpoint
            finite_path=joinpath(folder,"finite-observations"*extension)
            LineCableModels.save(finite_observed,finite_path)
            finite_restored=LineCableModels.import_data(Val(:observed),finite_path)
            @test finite_restored.gridpoint.source_gridpoint==
                finite_observed.gridpoint.source_gridpoint
            rebound_path=joinpath(folder,"rebound-observations"*extension)
            LineCableModels.save(rebound_observed,rebound_path)
            rebound_restored=LineCableModels.import_data(Val(:observed),rebound_path)
            @test rebound_restored.gridpoint.source_gridpoint==
                rebound_observed.gridpoint.source_gridpoint
        end
    end
end
