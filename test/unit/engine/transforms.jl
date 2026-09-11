@testitem "Transforms / uncertain matrices / inferred round trip" tags=[:unit] setup=[
    EngineTestSupport,
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
        ModalTransformationProblem(parameters),
        ModalTransformationFormulation(:default)
    )
    @test_throws ArgumentError compute(ModalTransformationProblem(modal);
        options = (offdiagonal_tolerance = 1e-6,))
    inverse_problem = ModalTransformationProblem(modal)
    @test (@inferred computation_options(typeof(inverse_problem), (;))) == (;)
    @test_throws ArgumentError computation_options(typeof(inverse_problem), (unknown=true,))
    rebuilt=@inferred compute(inverse_problem)

    @test eltype(modal) === Complex{Measurement{Float64}}
    @test TestNumerics.isapprox_scaled(
        Measurements.value.(real.(rebuilt.Z.values)),
        Measurements.value.(real.(impedance))
    )
    @test any(!iszero, Measurements.uncertainty.(real.(modal.Z.values)))
end

@testitem "Transforms / tracked eigensystems / numerical invariants" tags=[:unit] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestNumerics,
    TestAssertions
] begin
    using LinearAlgebra
    const Transforms=LineCableModels.Transforms

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
        formulation=ModalTransformationFormulation(identifier)
        @test formula_id(formulation) === identifier
        @test description(formulation) == descriptions[identifier]
        transformed=@inferred compute(
            ModalTransformationProblem(parameters),
            formulation
        )
        maps=operators(transformed)
        @test size(maps.voltage) == (2, 2, 2)
        @test size(maps.current) == (2, 2, 2)
        @test domain(transformed) === ModalDomain
        @test TestAssertions.all_finite(maps.voltage)
        @test TestAssertions.all_finite(maps.current)
        @test TestAssertions.all_finite(transformed.Z.values)
        @test TestAssertions.all_finite(transformed.Y.values)
        for frequency_index in eachindex(frequencies)
            basis_matrix=transpose(maps.voltage[:, :, frequency_index])
            @test abs(det(basis_matrix)) > sqrt(eps(Float64))
            for matrix in (transformed.Z.values, transformed.Y.values)
                slice=@view matrix[:, :, frequency_index]
                @test norm(slice-Diagonal(diag(slice))) <=
                      1e-6*max(norm(slice), eps(Float64))
            end
        end
        rebuilt=@inferred compute(ModalTransformationProblem(transformed))
        @test rebuilt.Z.values ≈ impedance rtol = 1e-6
        @test rebuilt.Y.values ≈ admittance rtol = 1e-6
    end

    modal=compute(
        ModalTransformationProblem(parameters),
        ModalTransformationFormulation(:default)
    )
    gamma=Transforms.gamma(modal)
    @test size(gamma) == size(impedance)
    @test all(index -> isdiag(gamma[:, :, index]), eachindex(frequencies))

    modal_products=Transforms.modal_quantities(modal)
    @test length(modal_products) == 6
    @test all(values -> size(values) == size(impedance), modal_products)
    @test Transforms._hungarian([4.0 1.0 3.0; 2.0 0.0 5.0; 3.0 2.0 2.0]) ==
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
            ModalTransformationProblem(smooth_parameters),
            ModalTransformationFormulation(identifier)
        )
        maps=operators(transformed)
        for frequency_index in 2:length(smooth_frequencies), mode in 1:2

            previous=@view maps.voltage[mode, :, frequency_index - 1]
            current=@view maps.voltage[mode, :, frequency_index]
            overlap=abs(dot(previous, current))/(norm(previous)*norm(current))
            @test overlap > 0.99
        end
    end

    for controls in ((max_iterations = 0,), (convergence = -1,))
        @test_throws ArgumentError ModalTransformationFormulation(
            formula(:default; options = (iteration = controls,)))
    end
end

@testitem "Engine / reduction / reorder, Kron, and bundle invariants" tags=[:unit] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestNumerics
] begin
    using LinearAlgebra
    const Engine=LineCableModels.Engine

    phase_map=[2, 0, 1, 2, 0, 1]
    @test Engine.reorder_indices(phase_map) == [1, 3, 4, 6, 2, 5]
    source=reshape(ComplexF64.(1:36), 6, 6)
    reordered, reordered_map=Engine.reorder_M(source, phase_map)
    @test reordered_map == [2, 1, 2, 1, 0, 0]
    @test reordered == source[[1, 3, 4, 6, 2, 5], [1, 3, 4, 6, 2, 5]]
    @test_throws ArgumentError Engine.reorder_M(ones(2, 3), [1, 2])

    matrix=ComplexF64[4 1 2; 1 5 3; 2 3 8]
    reduction_map=[1, 2, 0]
    expected=matrix[1:2, 1:2]-matrix[1:2, 3:3]*
                              inv(matrix[3:3, 3:3])*matrix[3:3, 1:2]
    @test TestNumerics.isapprox_scaled(kronify(matrix, reduction_map), expected)
    destination=zeros(ComplexF64, 2, 2)
    @test Engine.kronify!(matrix, reduction_map, destination) === nothing
    @test TestNumerics.isapprox_scaled(destination, expected)

    bundled, merged_map=Engine.merge_bundles!(copy(matrix), [1, 1, 0])
    @test merged_map == [1, 0, 0]
    change_of_basis=Matrix{ComplexF64}(I, 3, 3)
    change_of_basis[1, 2]=-1
    @test bundled == transpose(change_of_basis) * matrix * change_of_basis

    unconnected=ComplexF64[4 1 2; 1 5 3; 2 3 8]
    unchanged, unconnected_map=Engine.merge_bundles!(copy(unconnected), [0, 0, 1])
    @test unchanged == unconnected
    @test unconnected_map == [0, 0, 1]

    mixed=ComplexF64[6 1 2 3; 1 7 4 5; 2 4 8 6; 3 5 6 9]
    mixed_basis=Matrix{ComplexF64}(I, 4, 4)
    mixed_basis[1, 2]=-1
    mixed_result, mixed_map=Engine.merge_bundles!(copy(mixed), [2, 2, 0, 0])
    @test mixed_result == transpose(mixed_basis)*mixed*mixed_basis
    @test mixed_map == [2, 0, 0, 0]
    @test_throws ArgumentError Engine.merge_bundles!(ones(2, 3), [1, 1])
end

@testitem "Transforms / independent maps preserve ordered nonreciprocal entries" tags=[:unit] begin
    const TR = LineCableModels.Transforms
    const FM = LineCableModels.FormulaMethod
    Z = reshape(ComplexF64[2+3im 0.2+0.1im; 0.7+0.3im 4+5im], 2, 2, 1)
    Y = reshape(ComplexF64[2+6im -0.4-0.5im; -0.2-0.1im 3+8im] .* 1e-9, 2, 2, 1)
    A = ComplexF64[1 0.3; 0.1im 2]
    B = ComplexF64[2 0.2im; 0.4 1]
    replacement = (parameters,
        physical,
        options,
        workspace) -> TR.ModalOperators(reshape(copy(A), 2, 2, 1), reshape(copy(B), 2, 2, 1))
    @eval LineCableModels.computation_options(
        ::FM{:default, typeof(TR.modal_operators)}, ::$(typeof(replacement))) = (;)
    phase = LineParameters(PhaseDomain, Z, Y, [50.0])
    modal = compute(ModalTransformationProblem(phase),
        ModalTransformationFormulation(formula(:default; hooks = (contribution = replacement,)));
        options = (offdiagonal_tolerance = 2.0,))
    @test modal.Z.values[:, :, 1] ≈ A * Z[:, :, 1] / B
    @test modal.Y.values[:, :, 1] ≈ B * Y[:, :, 1] / A
    @test details(modal).modal.modified
    @test isempty(details(modal).modal.options)
    rebuilt = compute(ModalTransformationProblem(modal))
    @test rebuilt.Z.values ≈ Z
    @test rebuilt.Y.values ≈ Y
    @test phase.Z.values == Z && phase.Y.values == Y
end
