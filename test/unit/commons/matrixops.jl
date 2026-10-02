@testitem "Commons / matrix reductions / reorder, Kron, bundle and transposition invariants" tags=[:unit] setup=[
    UseEngineSupport,
    TestNumerics
] begin
    using LinearAlgebra
    const Commons=LineCableModels.Commons

    phase_map=[2, 0, 1, 2, 0, 1]
    @test Commons.reorder_indices(phase_map) == [1, 3, 4, 6, 2, 5]
    matrix=ComplexF64[4 1 2; 1 5 3; 2 3 8]
    reduction_map=[1, 2, 0]
    expected=matrix[1:2, 1:2]-matrix[1:2, 3:3]*
                              inv(matrix[3:3, 3:3])*matrix[3:3, 1:2]
    @test TestNumerics.isapprox_scaled(kron_reduce(matrix, reduction_map), expected)
    destination=zeros(ComplexF64, 2, 2)
    @test Commons.kron_reduce!(matrix, reduction_map, destination) === nothing
    @test TestNumerics.isapprox_scaled(destination, expected)

    # Complex asymmetric entries and ordered, noncontiguous indices expose
    # accidental conjugation, symmetry assumptions, and reordered terminals.
    for T in (Float32, Float64, BigFloat)
        ordered = Complex{T}[8+im 1-2im 2+im 3; 2+3im 9-im 1 2im;
                            1 3+im 10+2im 2; 4im 1-im 3+im 11]
        keep, eliminate = [4, 1], [3, 2]
        reduced = zeros(Complex{T}, 2, 2)
        actual = Commons.kron_reduce!(ordered, keep, eliminate, reduced,
            similar(reduced), similar(reduced), similar(reduced))
        expected_ordered = ordered[keep, keep] - ordered[keep, eliminate] *
            (ordered[eliminate, eliminate] \ ordered[eliminate, keep])
        @test actual === reduced
        @test actual ≈ expected_ordered
        @test kron_reduce(ordered, [2, 0, 3, 0]) ≈ ordered[[1, 3], [1, 3]] -
            ordered[[1, 3], [2, 4]] * (ordered[[2, 4], [2, 4]] \ ordered[[2, 4], [1, 3]])
        @test kron_reduce(ordered, [1, 2, 3, 4]) == ordered
        aliased = copy(ordered)
        @test Commons.kron_reduce!(aliased, [1, 2, 3, 4], aliased) === nothing
        @test aliased == ordered
    end
    @test_throws SingularException kron_reduce(ComplexF64[1 2; 3 0], [1, 0])
    @test_throws DimensionMismatch kron_reduce(matrix, [1, 0])

    bundled, merged_map=Commons.merge_bundles!(copy(matrix), [1, 1, 0])
    @test merged_map == [1, 0, 0]
    change_of_basis=Matrix{ComplexF64}(I, 3, 3)
    change_of_basis[1, 2]=-1
    @test bundled == transpose(change_of_basis) * matrix * change_of_basis

    unconnected=ComplexF64[4 1 2; 1 5 3; 2 3 8]
    unchanged, unconnected_map=Commons.merge_bundles!(copy(unconnected), [0, 0, 1])
    @test unchanged == unconnected
    @test unconnected_map == [0, 0, 1]

    mixed=ComplexF64[6 1 2 3; 1 7 4 5; 2 4 8 6; 3 5 6 9]
    mixed_basis=Matrix{ComplexF64}(I, 4, 4)
    mixed_basis[1, 2]=-1
    mixed_result, mixed_map=Commons.merge_bundles!(copy(mixed), [2, 2, 0, 0])
    @test mixed_result == transpose(mixed_basis)*mixed*mixed_basis
    @test mixed_map == [2, 0, 0, 0]
    @test_throws ArgumentError Commons.merge_bundles!(ones(2, 3), [1, 1])
    @test Commons.bundle_operations([2, 0, 1, 2, 1, 2]) == [(1, 4), (3, 5), (1, 6)]

    # Each option combination fixes the retained rows of the reordered terminals
    # [1, 3, 4, 6, 2, 5] with phases [2, 1, 2, 1, 0, 0].
    for (reduce_bundle, kron_reduction, keep, eliminate, retained, bundles) in (
            (true, true, [1, 2], [3, 4, 5, 6], [2, 1], [(1, 3), (2, 4)]),
            (true, false, [1, 2, 5, 6], [3, 4], [2, 1, -1, -1], [(1, 3), (2, 4)]),
            (false, true, [1, 2, 3, 4], [5, 6], [2, 1, 2, 1], Tuple{Int, Int}[]),
            (false, false, collect(1:6), Int[], [2, 1, 2, 1, 0, 0], Tuple{Int, Int}[]))
        plan=Commons.ReductionPlan(phase_map; reduce_bundle, kron_reduction,
            ideal_transposition = false)
        @test plan.permutation == [1, 3, 4, 6, 2, 5]
        @test (plan.keep, plan.eliminate, plan.phase_map, plan.bundles) ==
              (keep, eliminate, retained, bundles)
        @test plan.indices == plan.permutation[keep]
        buffers=Commons.ReductionBuffers{ComplexF64}(plan)
        @test size(buffers.potential) == (length(keep), length(keep))
        @test size(buffers.factor) == (length(eliminate), length(eliminate))
    end

    matrix=[Float64(2i+3j+i*j) for i in 1:3, j in 1:3]
    circulant=copy(matrix)
    @test Commons.ideal_transposition!(circulant) === circulant
    @test all(circulant[i, j] == circulant[mod1(i + 1, 3), mod1(j + 1, 3)]
    for i in 1:3, j in 1:3)
    @test sum(circulant) ≈ sum(matrix)
    @test_throws DimensionMismatch Commons.ideal_transposition!(ones(2, 3))
end

@testitem "Commons / matrix reductions / passive network constrained solves" tags=[:unit] begin
    using LinearAlgebra
    const Commons=LineCableModels.Commons
    incidence=[1.0 0 0 1 1 0;0 1 0 -1 0 1;0 0 1 0 -1 -1]
    for f in (10.0,100.0,1000.0)
        s=2pi*im*f
        r=collect(1:6).*1e-3;l=collect(7:12).*1e-6
        g=collect(1:2:11).*1e-9;c=collect(2:2:12).*1e-10
        z=inv(incidence*Diagonal(inv.(r.+s.*l))*transpose(incidence))
        y=incidence*Diagonal(g.+s.*c)*transpose(incidence)
        p=s*inv(y)
        # Set the eliminated conductor's voltage or potential to zero in the
        # original equations. Solve for the constrained currents or charges
        # without constructing a Schur complement.
        for matrix in (z,p)
            excitation=Matrix{ComplexF64}(I,3,3)[:,1:2]
            currents=matrix\excitation
            expected=inv(currents[1:2,:])
            @test kron_reduce(matrix,[1,2,0]) ≈ expected rtol=1e-10 atol=0
        end
        function reduced(phase_map; reduce_bundle)
            plan=Commons.ReductionPlan(phase_map; reduce_bundle, kron_reduction=true,
                ideal_transposition=false)
            Zr=zeros(ComplexF64,2,2);Yr=similar(Zr)
            Commons.reduce_line_matrices!(Zr,Yr,z,p,s,plan,Commons.ReductionBuffers{ComplexF64}(plan))
            return Zr,Yr
        end
        Zr,Yr=reduced([1,2,0];reduce_bundle=false)
        @test Zr ≈ inv((z\Matrix{ComplexF64}(I,3,3)[:,1:2])[1:2,:]) rtol=1e-10
        @test Yr ≈ s*((p\Matrix{ComplexF64}(I,3,3)[:,1:2])[1:2,:]) rtol=1e-10
        Zb,Yb=reduced([1,1,2];reduce_bundle=true)
        equal_potentials=[1.0 0;1 0;0 1]
        @test Zb ≈ inv(transpose(equal_potentials)*(z\equal_potentials)) rtol=1e-10
        @test Yb ≈ s*transpose(equal_potentials)*(p\equal_potentials) rtol=1e-10
    end
end

@testitem "Commons / matrix reductions / admittance inversion diagnostics" tags=[:unit] begin
    using LinearAlgebra
    const Commons=LineCableModels.Commons
    plan=Commons.ReductionPlan([1, 2]; reduce_bundle=false, kron_reduction=false,
        ideal_transposition=false)
    buffers=Commons.ReductionBuffers{ComplexF64}(plan)
    Zr, Yr=zeros(ComplexF64, 2, 2), zeros(ComplexF64, 2, 2)
    s=2pi*50im
    for potential in (ComplexF64[2.0 0.25; 0.25 1.5], ComplexF64[3.0+0.1im 0.5; 0.5 2.0+0.2im])
        record=Commons.reduce_line_matrices!(Zr, Yr, potential, potential, s, plan, buffers,
            Val(true))
        @test Zr == potential
        @test potential * Yr ≈ s * I
        @test record.residual ≤ 100eps(Float64)
        @test Commons.reduce_line_matrices!(Zr, Yr, potential, potential, s, plan,
            buffers) === nothing
        @test potential * Yr ≈ s * I
    end
    # An overflowing condition estimate is not a failed solve. The physical
    # inverse is finite and remains the result of the original LU operation.
    ill_conditioned=ComplexF64[1e-200 0; 0 1e200]
    record=@test_logs (:warn, r"condition estimate is not finite") begin
        Commons.reduce_line_matrices!(Zr, Yr, ill_conditioned, ill_conditioned, one(s), plan,
            buffers, Val(true))
    end
    @test Yr == inv(ill_conditioned)
    @test isinf(record.condition_number)
    @test record.residual ≤ 100eps(Float64)
    @test_throws ArgumentError Commons.reduce_line_matrices!(Zr, Yr, ill_conditioned,
        fill(ComplexF64(NaN), 2, 2), s, plan, buffers, Val(true))
    @test_throws SingularException Commons.reduce_line_matrices!(Zr, Yr, ill_conditioned,
        zeros(ComplexF64, 2, 2), s, plan, buffers, Val(true))
end
