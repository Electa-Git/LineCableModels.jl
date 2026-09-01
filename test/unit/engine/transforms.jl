@testitem "Transforms / Fortescue / unitary modal invariants" tags=[:unit] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestNumerics
] begin
    using LinearAlgebra
    const Transforms=LineCableModels.Transforms

    frequency_values=[50.0, 500.0]
    impedance_slice=ComplexF64[3+9im 1+2im 1+2im
                               1+2im 3+9im 1+2im
                               1+2im 1+2im 3+9im]
    admittance_slice=ComplexF64[6+12im -2-3im -2-3im
                                -2-3im 6+12im -2-3im
                                -2-3im -2-3im 6+12im] .* 1e-9
    impedance=cat(impedance_slice, 2 .* impedance_slice; dims = 3)
    admittance=cat(admittance_slice, 2 .* admittance_slice; dims = 3)
    parameters=LineParameters(
        PhaseDomain,
        impedance,
        admittance,
        frequency_values;
        details = (source = :test,)
    )

    transform=Transforms.modal_basis(Val(:Fortescue), 3)
    @test TestNumerics.isapprox_scaled(
        transform * transform',
        Matrix{ComplexF64}(I, 3, 3)
    )
    @test_throws ArgumentError Transforms.modal_basis(Val(:Fortescue), 0)

    formulation=ModalTransformationFormulation(:Fortescue)
    selected=@inferred ModalTransformationFormulation(
        formula(:Chrysochos2014; tolerance=1e-7)
    )
    default_formulation=@inferred ModalTransformationFormulation()
    @test Transforms.DEFAULT === :Chrysochos2014
    @test formula_id(default_formulation) === :Chrysochos2014
    @test formula_id(Transforms.Formula(:default)) === :Chrysochos2014
    @test :default ∉ Transforms.formulas()
    @test description(formulation) ==
          "Fortescue symmetrical-component transformation (1918)"
    @test formula_id(formulation) === :Fortescue
    @test formula_id(selected) === :Chrysochos2014
    @test LineCableModels.Transforms.assumptions(selected).tolerance == 1e-7
    @test Transforms.formulas() == (
        :Chrysochos2014,
        :Fan2009,
        :Fortescue,
        :Wedepohl1996
    )
    modal=@inferred compute(ModalTransformationProblem(parameters), formulation)
    explicit=compute(
        LineCableModelsModal(),
        ModalTransformationProblem(parameters),
        formulation
    )
    @test Z(explicit) == Z(modal)
    @test Y(explicit) == Y(modal)
    maps=operators(modal)
    @test maps.voltage === maps.current
    @test maps.voltage[:, :, 1] == transform
    @test maps.voltage[:, :, 2] == transform
    @test domain(modal) === ModalDomain
    @test basis(modal) === :pul
    @test frequencies(modal) == frequency_values
    @test details(modal) == (source = :test,)
    for frequency_index in eachindex(frequency_values)
        for matrix in (modal.Z.values, modal.Y.values)
            slice=@view matrix[:, :, frequency_index]
            @test norm(slice-Diagonal(diag(slice))) <=
                  1e-12*max(norm(slice), eps(Float64))
        end
        @test tr(modal.Z.values[:, :, frequency_index]) ≈
              tr(impedance[:, :, frequency_index])
        @test tr(modal.Y.values[:, :, frequency_index]) ≈
              tr(admittance[:, :, frequency_index])
    end
    rebuilt=@inferred compute(ModalTransformationProblem(modal))
    @test domain(rebuilt) === PhaseDomain
    @test TestNumerics.isapprox_scaled(Z(rebuilt), impedance)
    @test TestNumerics.isapprox_scaled(Y(rebuilt), admittance)
    @test @allocated(compute(ModalTransformationProblem(parameters), formulation)) <=
          8_192
    @test @allocated(compute(ModalTransformationProblem(modal))) <= 8_192
    sliced=modal[1]
    @test size(operators(sliced).voltage) == (3, 3, 1)
    @test details(sliced) == details(modal)
    @test_throws MethodError LineParameters(
        ModalDomain,
        impedance,
        admittance,
        frequency_values
    )
    @test_throws ArgumentError ModalTransformationFormulation(
        :Fortescue;
        unsupported = true
    )
    @test_throws ArgumentError ModalTransformationFormulation(
        formula(:Fortescue; order=:before)
    )
    @test_throws DomainError ModalTransformationFormulation(
        :Fortescue;
        tolerance = -1.0
    )

    custom_route=(
        source, _)->begin
        source_maps=operators(modal)
        ModalOperators(copy(source_maps.voltage), copy(source_maps.current))
    end
    experiment=ModalTransformationFormulation(Transforms.Formula(
        :Experiment,
        custom_route,
        (tolerance = 1e-4,)
    ))
    experimental=compute(ModalTransformationProblem(parameters), experiment)
    @test formula_id(experimental.domain.formula) === :Experiment
    @test TestNumerics.isapprox_scaled(Z(experimental), Z(modal))
    @test TestNumerics.isapprox_scaled(Y(experimental), Y(modal))

    override=ModalTransformationFormulation(
        :Fortescue;
        route = custom_route
    )
    overridden=compute(ModalTransformationProblem(parameters), override)
    @test formula_id(overridden.domain.formula) === :Fortescue
    @test Z(overridden) == Z(experimental)
    @test Y(overridden) == Y(experimental)
end

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
        ModalTransformationFormulation(:Fortescue)
    )
    rebuilt=@inferred compute(ModalTransformationProblem(modal))

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

    descriptions=(
        Chrysochos2014 =
            "Chrysochos et al. Levenberg–Marquardt modal transformation (2014)",
        Fan2009 = "Fan et al. eigenvector-tracking transformation (2009)",
        Wedepohl1996 = "Wedepohl et al. Newton–Raphson modal transformation (1996)"
    )

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
        ModalTransformationFormulation(:Chrysochos2014)
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
            2+0.01frequency_index+3im,
            4+0.02frequency_index+5im
        ])
        modal_admittance=Diagonal(ComplexF64[
            (4+0.01frequency_index)+8im,
            (8+0.02frequency_index)+12im
        ].*1e-9)
        smooth_impedance[:, :, frequency_index]=
            basis_matrix*modal_impedance*transpose(basis_matrix)
        smooth_admittance[:, :, frequency_index]=
            basis_matrix*modal_admittance*transpose(basis_matrix)
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
            previous=@view maps.voltage[mode, :, frequency_index-1]
            current=@view maps.voltage[mode, :, frequency_index]
            overlap=abs(dot(previous, current))/(norm(previous)*norm(current))
            @test overlap > 0.99
        end
    end

    @test_throws DomainError compute(
        ModalTransformationProblem(parameters),
        ModalTransformationFormulation(:Chrysochos2014; max_iterations=0)
    )
    @test_throws DomainError compute(
        ModalTransformationProblem(parameters),
        ModalTransformationFormulation(:Chrysochos2014; convergence=-1)
    )
    @test_throws DomainError compute(
        ModalTransformationProblem(parameters),
        ModalTransformationFormulation(:Fan2009; history_weight=-1)
    )
    @test_throws DomainError compute(
        ModalTransformationProblem(parameters),
        ModalTransformationFormulation(:Fan2009; coalescence_tolerance=-1)
    )
    @test_throws DomainError compute(
        ModalTransformationProblem(parameters),
        ModalTransformationFormulation(:Wedepohl1996; max_iterations=0)
    )
    @test_throws DomainError compute(
        ModalTransformationProblem(parameters),
        ModalTransformationFormulation(:Wedepohl1996; convergence=-1)
    )
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
