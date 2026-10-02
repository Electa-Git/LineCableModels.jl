@testitem "ModalAnalysis / uncertain matrices / inferred round trip" tags=[:unit] setup=[
    UseEngineSupport,
    TestNumerics
] begin
    using Measurements
    using Test

    impedance=Array{Complex{Measurement{Float64}}, 3}(undef, 2, 2, 1)
    admittance=similar(impedance)
    for column in 1:2, row in 1:2

        diagonal=row==column
        impedance[row, column, 1]=complex(
            measurement(diagonal ? 3.0 : 1.0, 1e-3),
            measurement(diagonal ? 5.0 : 0.5, 1e-3)
        )
        admittance[row, column, 1]=complex(
            measurement(diagonal ? 3e-9 : 1e-9, 1e-12),
            measurement(diagonal ? 5e-9 : 0.5e-9, 1e-12)
        )
    end
    parameters=LineParameters(PhaseDomain, impedance, admittance, [50.0])
    modal=@inferred compute(
        ModalAnalysisProblem(parameters),
        ModalAnalysisFormulation(:default)
    )
    rebuilt=LineCableModels.ModalAnalysis.transform(PhaseDomain,modal)

    @test eltype(modal) === Complex{Measurement{Float64}}
    @test TestNumerics.isapprox_scaled(
        Measurements.value.(real.(rebuilt.Z.values)),
        Measurements.value.(real.(impedance))
    )
    @test any(!iszero, Measurements.uncertainty.(real.(modal.Z.values)))
    @test any(!iszero, Measurements.uncertainty.(real.(gamma(modal))))
    segment=PropagationParameters(modal;line_length=10.0)
    @test any(!iszero, Measurements.uncertainty.(real.(H(segment))))
end

@testitem "ModalAnalysis / tracked eigensystems / numerical invariants" tags=[:unit] setup=[
    UseEngineSupport,
    TestNumerics,
    TestAssertions
] begin
    using LinearAlgebra
    const ModalAnalysis=LineCableModels.ModalAnalysis

    descriptions=(default = "Chrysochos et al. Levenberg–Marquardt modal transformation (2014)",)

    frequencies=[50.0, 100.0]
    impedance=zeros(ComplexF64, 2, 2, 2)
    admittance=zeros(ComplexF64, 2, 2, 2)
    impedance[:, :, 1]=[2+4im 0.2+0.1im; 0.2+0.1im 3+5im]
    impedance[:, :, 2]=[2.2+8im 0.21+0.2im; 0.21+0.2im 3.2+10im]
    admittance[:, :, 1]=[4+8im -0.5-1im; -0.5-1im 5+9im] .* 1e-9
    admittance[:, :, 2]=[4.2+16im -0.55-2im; -0.55-2im 5.2+18im] .* 1e-9
    parameters=LineParameters(PhaseDomain, impedance, admittance, frequencies)

    for identifier in keys(descriptions)
        formulation=ModalAnalysisFormulation(identifier)
        @test formula_id(formulation) === :modal
        @test formula_id(formulation.formula) === :chrysochos2014
        @test description(formulation.formula) == descriptions[identifier]
        transformed=@inferred compute(
            ModalAnalysisProblem(parameters),
            formulation
        )
        maps=operators(transformed)
        @test size(maps.Tv) == (2, 2, 2)
        @test size(maps.Ti) == (2, 2, 2)
        @test domain(transformed) === ModalDomain
        @test TestAssertions.all_finite(maps.Tv)
        @test TestAssertions.all_finite(maps.Ti)
        @test TestAssertions.all_finite(transformed.Z.values)
        @test TestAssertions.all_finite(transformed.Y.values)
        for frequency_index in eachindex(frequencies)
            basis_matrix=maps.Ti[:, :, frequency_index]
            @test abs(det(basis_matrix)) > sqrt(eps(Float64))
            for matrix in (transformed.Z.values, transformed.Y.values)
                slice=@view matrix[:, :, frequency_index]
                @test norm(slice-Diagonal(diag(slice))) <=
                      1e-6*max(norm(slice), eps(Float64))
            end
        end
        rebuilt=LineCableModels.ModalAnalysis.transform(PhaseDomain,transformed)
        @test rebuilt.Z.values ≈ impedance rtol = 1e-6
        @test rebuilt.Y.values ≈ admittance rtol = 1e-6
    end

    modal=compute(
        ModalAnalysisProblem(parameters),
        ModalAnalysisFormulation(:default)
    )
    roots=ModalAnalysis.gamma(modal)
    @test size(roots)==(2,length(frequencies))
    @test all(isfinite,roots)
    @test size(Zc(modal))==size(roots)
    @test size(Yc(modal))==size(roots)
    @test size(Zc(modal,PhaseDomain))==size(impedance)
    @test size(Yc(modal,PhaseDomain))==size(impedance)
    assignment_work=ModalAnalysis._assignment_workspace(ComplexF64,3)
    @test ModalAnalysis.hungarian_assignment!(
        [4.0 1.0 3.0; 2.0 0.0 5.0; 3.0 2.0 2.0],assignment_work) ==
          [2, 1, 3]

    smooth_frequencies=collect(10.0 .^ range(1, 4; length = 17))
    smooth_impedance=zeros(ComplexF64, 2, 2, length(smooth_frequencies))
    smooth_admittance=similar(smooth_impedance)
    for frequency_index in eachindex(smooth_frequencies)
        angle=(frequency_index-1)/(length(smooth_frequencies)-1)*(pi/3)
        basis_matrix=[cos(angle) -sin(angle); sin(angle) cos(angle)]
        modal_impedance=Diagonal(ComplexF64[
        2 + 0.01frequency_index + 3im,
        4 + 0.02frequency_index + 5im
])
        modal_admittance=Diagonal(ComplexF64[
            (4 + 0.01frequency_index) + 8im,
            (8 + 0.02frequency_index) + 12im
        ] .* 1e-9)
        smooth_impedance[:, :, frequency_index]=basis_matrix*modal_impedance*transpose(basis_matrix)
        smooth_admittance[:, :, frequency_index]=basis_matrix*modal_admittance*transpose(basis_matrix)
    end
    smooth_parameters=LineParameters(
        PhaseDomain,
        smooth_impedance,
        smooth_admittance,
        smooth_frequencies
    )
    for identifier in keys(descriptions)
        transformed=compute(
            ModalAnalysisProblem(smooth_parameters),
            ModalAnalysisFormulation(identifier)
        )
        maps=operators(transformed)
        for frequency_index in 2:length(smooth_frequencies), mode in 1:2

            previous=@view maps.Ti[:, mode, frequency_index - 1]
            current=@view maps.Ti[:, mode, frequency_index]
            overlap=abs(dot(previous, current))/(norm(previous)*norm(current))
            @test overlap > 0.99
        end
    end

    for controls in ((max_iterations = 0,), (convergence = -1,))
        @test_throws ArgumentError ModalAnalysisFormulation(
            formula(:default; options = (iteration = controls,)))
    end
end

@testitem "ModalAnalysis / independent maps preserve ordered nonreciprocal entries" tags=[:unit] setup=[FormulaFixtures] begin
    const TR = LineCableModels.ModalAnalysis
    const FM = LineCableModels.FormulaMethod
    Z = reshape(ComplexF64[2+3im 0.2+0.1im; 0.7+0.3im 4+5im], 2, 2, 1)
    Y = reshape(ComplexF64[2+6im -0.4-0.5im; -0.2-0.1im 3+8im] .* 1e-9, 2, 2, 1)
    A = ComplexF64[1 0.3; 0.1im 2]
    B = ComplexF64[2 0.2im; 0.4 1]
    selected=FormulaFixtures.FixedModalMaps(reshape(copy(A),2,2,1),reshape(copy(B),2,2,1))
    phase = LineParameters(PhaseDomain, Z, Y, [50.0])
    modal = compute(ModalAnalysisProblem(phase),
        ModalAnalysisFormulation(selected);
        options = (offdiagonal_tolerance = 2.0,))
    @test modal.Z.values[:, :, 1] ≈ A \ (Z[:, :, 1] * B)
    @test modal.Y.values[:, :, 1] ≈ B \ (Y[:, :, 1] * A)
    @test details(modal).data.modal.requested.identifier === :FixedModalMaps
    @test details(modal).data.modal.effective.identifier === :FixedModalMaps
    @test details(modal).data.modal.identifier === :FixedModalMaps
    rebuilt = LineCableModels.ModalAnalysis.transform(PhaseDomain,modal)
    @test rebuilt.Z.values ≈ Z
    @test rebuilt.Y.values ≈ Y
    @test phase.Z.values == Z && phase.Y.values == Y
end

@testitem "ModalAnalysis / current commuting modes / eigenvalues and reconstruction" tags=[:unit] begin
    using LinearAlgebra
    f=[10.0,100.0,1000.0];theta=pi/6
    q=[cos(theta) -sin(theta);sin(theta) cos(theta)]
    zs=[q*Diagonal([1+s*1e-3,3+s*2e-3])*transpose(q) for s in 2pi*im.*f]
    ys=[q*Diagonal([1e-6+s*1e-9,3e-6+s*2e-9])*transpose(q) for s in 2pi*im.*f]
    source=LineParameters(cat(zs...;dims=3),cat(ys...;dims=3),f)
    modes=compute(ModalAnalysisProblem(source),ModalAnalysisFormulation())
    restored=LineCableModels.ModalAnalysis.transform(PhaseDomain,modes)
    @test Z(restored) ≈ Z(source) rtol=1e-10
    @test Y(restored) ≈ Y(source) rtol=1e-10
    for k in eachindex(f)
        s=2pi*im*f[k]
        expected=[(1+s*1e-3)*(1e-6+s*1e-9),(3+s*2e-3)*(3e-6+s*2e-9)]
        actual=eigvals(Z(modes)[:,:,k]*Y(modes)[:,:,k])
        @test sort(actual;by=abs) ≈ sort(expected;by=abs) rtol=1e-10
    end
end
