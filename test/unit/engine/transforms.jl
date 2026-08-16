@testitem "Engine / Fortescue / unitary modal invariants" tags=[:unit] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestNumerics
] begin
    using LinearAlgebra
    const Transforms=LineCableModels.Engine.Transforms

    frequency_values=[50.0, 500.0]
    impedance_slice=ComplexF64[3+9im 1+2im 1+2im
                               1+2im 3+9im 1+2im
                               1+2im 1+2im 3+9im]
    admittance_slice=ComplexF64[6+12im -2-3im -2-3im
                                -2-3im 6+12im -2-3im
                                -2-3im -2-3im 6+12im] .* 1e-9
    impedance=cat(impedance_slice, 2 .* impedance_slice; dims = 3)
    admittance=cat(admittance_slice, 2 .* admittance_slice; dims = 3)
    parameters=LineParameters(PhaseDomain, impedance, admittance, frequency_values)

    transform=Transforms.fortescue_F(3)
    @test TestNumerics.isapprox_scaled(
        transform * transform',
        Matrix{ComplexF64}(I, 3, 3)
    )
    @test_throws ArgumentError Transforms.fortescue_F(0)

    returned_transform, modal=Fortescue()(parameters)
    @test returned_transform == transform
    @test domain(modal) === ModalDomain
    @test basis(modal) === :per_length
    @test frequencies(modal) == frequency_values
    for frequency_index in eachindex(frequency_values)
        @test isdiag(modal.Z.values[:, :, frequency_index])
        @test isdiag(modal.Y.values[:, :, frequency_index])
        @test tr(modal.Z.values[:, :, frequency_index]) ≈
              tr(impedance[:, :, frequency_index])
        @test tr(modal.Y.values[:, :, frequency_index]) ≈
              tr(admittance[:, :, frequency_index])
    end
    @test_throws ErrorException Fortescue()(modal)
end

@testitem "Engine / Levenberg / eigen and modal residuals" tags=[:unit] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestNumerics,
    TestAssertions
] begin
    using LinearAlgebra
    const Transforms=LineCableModels.Engine.Transforms

    @test get_description(Transforms.Levenberg()) ==
          "Levenberg–Marquardt (frequency-tracked eigen decomposition)"
    imaginary=ComplexF64[1 + 2im, -3 + 4im]
    @test Transforms.imag!(imaginary) === imaginary
    @test imaginary == ComplexF64[2, 4]

    frequencies=[50.0, 100.0]
    impedance=zeros(ComplexF64, 2, 2, 2)
    admittance=zeros(ComplexF64, 2, 2, 2)
    impedance[:, :, 1]=[2+4im 0.2+0.1im; 0.2+0.1im 3+5im]
    impedance[:, :, 2]=[2.2+8im 0.21+0.2im; 0.21+0.2im 3.2+10im]
    admittance[:, :, 1]=[4+8im -0.5-1im; -0.5-1im 5+9im] .* 1e-9
    admittance[:, :, 2]=[4.2+16im -0.55-2im; -0.55-2im 5.2+18im] .* 1e-9
    parameters=LineParameters(PhaseDomain, impedance, admittance, frequencies)

    transformation, modal=Transforms.Levenberg(tol = 1e-7)(parameters)
    @test size(transformation) == (2, 2, 2)
    @test domain(modal) === ModalDomain
    @test TestAssertions.all_finite(transformation)
    @test TestAssertions.all_finite(modal.Z.values)
    @test TestAssertions.all_finite(modal.Y.values)

    for frequency_index in eachindex(frequencies)
        basis_matrix=transformation[:, :, frequency_index]
        @test abs(det(basis_matrix)) > sqrt(eps(Float64))
        @test maximum(abs,
            modal.Z.values[:, :, frequency_index] -
            Diagonal(diag(modal.Z.values[:, :, frequency_index]))) < 1e-5
        @test maximum(abs,
            modal.Y.values[:, :, frequency_index] -
            Diagonal(diag(modal.Y.values[:, :, frequency_index]))) < 1e-12
    end

    gamma=Transforms._calc_gamma(transformation, impedance, admittance)
    @test size(gamma) == size(impedance)
    @test all(index -> isdiag(gamma[:, :, index]), eachindex(frequencies))

    modal_products=Transforms._calc_modal_quantities(
        transformation,
        impedance,
        admittance
    )
    @test length(modal_products) == 6
    @test all(values -> size(values) == size(impedance), modal_products)
    @test_throws DimensionMismatch Transforms.rot!(ones(ComplexF64, 2, 3))
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
    @test_throws ArgumentError Engine.merge_bundles!(ones(2, 3), [1, 1])
end
