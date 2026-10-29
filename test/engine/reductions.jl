using LinearAlgebra: I, lu!

@testset "Engine reductions" begin
    Engine = LineCableModels.Engine
    # A dense complex matrix, diagonally dominant so that every diagonal block is regular.
    matrix = [complex(sin(i + 2j), cos(3i - j)) for i in 1:5, j in 1:5] + 6I
    kept = [4, 1, 3]
    eliminated = [5, 2]
    reduced = similar(matrix, 3, 3)
    factor = similar(matrix, 2, 2)
    coupling = similar(matrix, 3, 2)
    right_hand_side = similar(matrix, 2, 3)
    schur = matrix[kept, kept] -
            matrix[kept, eliminated] * inv(matrix[eliminated, eliminated]) *
            matrix[eliminated, kept]
    @test isapprox(
        Engine.eliminate!(reduced, matrix, kept, eliminated, factor, coupling,
            right_hand_side),
        schur; rtol = 1e-13)
    # Without an eliminated row, the elimination copies the kept block in the order of
    # `kept`.
    @test Engine.eliminate!(similar(matrix, 3, 3), matrix, kept, Int[],
        similar(matrix, 0, 0), similar(matrix, 3, 0), similar(matrix, 0, 3)) ==
          matrix[kept, kept]
    # A bundle of the conductors 1, 3 and 4 changes the basis of the matrix by `basis`.
    pairs = [(1, 3), (1, 4)]
    basis = Matrix{ComplexF64}(I, 5, 5)
    basis[3, 1] = basis[4, 1] = -1
    @test isapprox(Engine.bundle!(copy(matrix), pairs), basis * matrix * transpose(basis);
        rtol = 1e-14)
    mean = sum(matrix[circshift(1:5, k), circshift(1:5, k)] for k in 0:4) / 5
    @test isapprox(Engine.transposition!(copy(matrix)), mean; rtol = 1e-14)
    @test_throws DimensionMismatch Engine.transposition!(similar(matrix, 2, 3))
    # After the first calls above, the bundle and the transposition do not allocate. The
    # elimination allocates what its one LU decomposition allocates: the pivot vector.
    target = copy(matrix)
    @test (@allocated Engine.bundle!(target, pairs)) == 0
    @test (@allocated Engine.transposition!(target)) == 0
    block = matrix[eliminated, eliminated]
    lu!(copy(block))
    pivots = @allocated lu!(block)
    @test (@allocated Engine.eliminate!(
        reduced, matrix, kept, eliminated, factor, coupling,
        right_hand_side)) == pivots
end
