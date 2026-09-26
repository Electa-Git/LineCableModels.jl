@testitem "ModalAnalysis / repeated eigenvalues preserve the modal subspace" tags=[:unit] begin
    using LinearAlgebra

    frequencies=[50.0,55.0,60.0,65.0,70.0]
    impedance=zeros(ComplexF64,3,3,length(frequencies))
    admittance=similar(impedance)
    principal_axes=Matrix{Float64}[]
    for (index,angle) in enumerate(range(0.0,0.3;length=length(frequencies)))
        rotation=[cos(angle) 0 sin(angle);0 1 0;-sin(angle) 0 cos(angle)]
        push!(principal_axes,rotation)
        impedance[:,:,index]=(1+2im).*(rotation*Diagonal([1.0,1.0,5.0])*transpose(rotation))
        admittance[:,:,index]=(1e-8+4e-7im).*Matrix(I,3,3)
    end
    phase=LineParameters(impedance,admittance,frequencies)
    modal=compute(ModalAnalysisProblem(phase),ModalAnalysisFormulation(:default))
    maps=operators(modal)
    for index in eachindex(frequencies)
        vectors=maps.Ti[:,:,index]
        product=admittance[:,:,index]*impedance[:,:,index]
        eigenvalues=gamma(modal)[:,index].^2
        @test norm(product*vectors-vectors*Diagonal(eigenvalues))<=
            1e-8*norm(product)*norm(vectors)
        for mode in 1:3
            paired=impedance[:,:,index]*vectors[:,mode]
            @test maps.Tv[:,mode,index]≈paired/norm(paired)
        end
        target=(1e-8+4e-7im)*(1+2im)
        repeated=findall(value->abs(value-target)<=1e-8*abs(target),eigenvalues)
        @test length(repeated)==2
        actual_basis=vectors[:,repeated]
        actual_projector=actual_basis*pinv(actual_basis)
        expected_basis=principal_axes[index][:,1:2]
        @test actual_projector≈expected_basis*transpose(expected_basis) rtol=1e-8
    end
    rebuilt=LineCableModels.ModalAnalysis.transform(PhaseDomain,modal)
    @test rebuilt.Z.values≈impedance rtol=1e-12
    @test rebuilt.Y.values≈admittance rtol=1e-12
    @test phase.Z.values==impedance
    @test phase.Y.values==admittance

    for selector in (2,2:4,[5,2],:)
        indices=selector isa Integer ? (selector:selector) : selector
        selected=modal[selector]
        selected_maps=operators(selected)
        @test selected.f==frequencies[indices]
        @test selected.Z.values==modal.Z.values[:,:,indices]
        @test selected.Y.values==modal.Y.values[:,:,indices]
        @test selected_maps.Tv==maps.Tv[:,:,indices]
        @test selected_maps.Ti==maps.Ti[:,:,indices]
        @test gamma(selected)==gamma(modal)[:,indices]
        selected_phase=LineCableModels.ModalAnalysis.transform(PhaseDomain,selected)
        @test selected_phase.Z.values≈impedance[:,:,indices] rtol=1e-12
        @test selected_phase.Y.values≈admittance[:,:,indices] rtol=1e-12
        segment=PropagationParameters(modal;line_length=25.0)
        selected_segment=segment[selector]
        @test selected_segment.parameters.f==selected.f
        @test gamma(selected_segment)==gamma(selected)
        @test H(selected_segment)==exp.(-gamma(selected).*25.0)
    end
    selected=modal[2:4]
    voltage_before=copy(maps.Tv)
    current_before=copy(maps.Ti)
    operators(selected).Tv[1,1,1]+=1
    operators(selected).Ti[1,1,1]+=1
    @test maps.Tv==voltage_before
    @test maps.Ti==current_before
end

@testitem "ModalAnalysis / limited iteration reports matched and unrecovered results" tags=[:unit] begin
    using LinearAlgebra
    frequencies=[50.0,100.0,200.0,400.0]
    impedance=zeros(ComplexF64,2,2,length(frequencies))
    admittance=similar(impedance)
    for (index,angle) in enumerate((0.0,0.4,0.8,1.2))
        rotation=[cos(angle) -sin(angle);sin(angle) cos(angle)]
        impedance[:,:,index]=rotation*Diagonal(ComplexF64[
            1+index*im,3+2index*im])*transpose(rotation)
        admittance[:,:,index]=rotation*Diagonal(ComplexF64[
            1e-8+index*1e-7im,2e-8+index*3e-7im])*transpose(rotation)
    end
    phase=LineParameters(impedance,admittance,frequencies)
    matched=compute(ModalAnalysisProblem(phase),ModalAnalysisFormulation(:default;
        options=(iteration=(max_iterations=1,),)))
    diagnostic=details(matched).data.modal.diagnostics
    @test !isempty(diagnostic.fallback_frequencies)
    @test diagnostic.fallback_frequencies==diagnostic.missed_frequencies
    @test size(diagnostic.eigen_residual)==(2,length(frequencies))
    @test all(isfinite,diagnostic.eigen_residual)
    @test all(value -> value===nothing || value isa Bool,diagnostic.converged)
    @test all(value -> value===nothing || value>=0,diagnostic.iterations)
    for index in eachindex(frequencies)
        vectors=Ti(matched)[:,:,index]
        product=admittance[:,:,index]*impedance[:,:,index]
        eigenvalues=gamma(matched)[:,index].^2
        @test norm(product*vectors-vectors*Diagonal(eigenvalues))<=
            1e-9*norm(product)*norm(vectors)
        @test all(isfinite,matched.Z.values[:,:,index])
        @test all(isfinite,matched.Y.values[:,:,index])
    end
    rebuilt=LineCableModels.ModalAnalysis.transform(PhaseDomain,matched)
    @test rebuilt.Z.values≈impedance rtol=1e-12
    @test rebuilt.Y.values≈admittance rtol=1e-12
    unassisted=compute(ModalAnalysisProblem(phase),ModalAnalysisFormulation(:default;
        options=(iteration=(max_iterations=1,fallback=:none),)))
    @test !isempty(details(unassisted).data.modal.diagnostics.missed_frequencies)
    @test isempty(details(unassisted).data.modal.diagnostics.fallback_frequencies)
    @test all(isfinite,unassisted.Z.values)
    @test all(isfinite,unassisted.Y.values)
    selected=unassisted[[4,2]]
    @test details(selected).data.modal.diagnostics.missed_frequencies==[1,2]
    @test details(selected).data.modal.diagnostics.z_coupling==
        details(unassisted).data.modal.diagnostics.z_coupling[[4,2]]
    @test details(selected).data.modal.diagnostics.eigen_residual==
        details(unassisted).data.modal.diagnostics.eigen_residual[:,[4,2]]
end
