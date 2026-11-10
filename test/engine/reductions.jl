using LinearAlgebra: BlasInt, I, SingularException, ldiv!, lu!

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
    pivots = Vector{BlasInt}(undef, 3)
    schur = matrix[kept, kept] -
            matrix[kept, eliminated] * inv(matrix[eliminated, eliminated]) *
            matrix[eliminated, kept]
    @test isapprox(
        Engine.eliminate!(reduced, matrix, kept, eliminated, factor, coupling,
            right_hand_side, pivots),
        schur; rtol = 1e-13)
    # Without an eliminated row, the elimination copies the kept block in the order of
    # `kept`.
    @test Engine.eliminate!(similar(matrix, 3, 3), matrix, kept, Int[],
        similar(matrix, 0, 0), similar(matrix, 3, 0), similar(matrix, 0, 3), pivots) ==
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
    # After the first calls above, the bundle, the transposition and the elimination do
    # not allocate. The LU decomposition of the elimination writes into the pivots of the
    # caller.
    target = copy(matrix)
    @test (@allocated Engine.bundle!(target, pairs)) == 0
    @test (@allocated Engine.transposition!(target)) == 0
    @test (@allocated Engine.eliminate!(
        reduced, matrix, kept, eliminated, factor, coupling,
        right_hand_side, pivots)) == 0
end

@testset "factorize! by the element type" begin
    Engine = LineCableModels.Engine
    pivots = Vector{BlasInt}(undef, 12)
    # Complex matrices of orders 2 to 12 with a zero first entry, so that the factorization
    # interchanges rows. On a BLAS element type the factors, the interchanges and the solves
    # are those of `lu!`, bit for bit.
    @testset "order $n" for n in 2:12
        A = [complex(sin(i + 2j), cos(3i - j)) for i in 1:n, j in 1:n] + n * I
        A[1, 1] = 0
        B = [complex(cos(i - j), sin(i * j)) for i in 1:n, j in 1:3]
        expected = lu!(copy(A))
        factors = Engine.factorize!(copy(A), pivots)
        @test reinterpret(UInt64, factors.factors) == reinterpret(UInt64, expected.factors)
        @test factors.ipiv == expected.ipiv
        @test reinterpret(UInt64, ldiv!(factors, copy(B))) ==
              reinterpret(UInt64, ldiv!(expected, copy(B)))
    end
    # After a first call, the factorization of a BLAS element type and the solve with it
    # allocate 0 bytes. The function and its arguments pass through an anonymous
    # function, where their types are concrete.
    A = [complex(sin(i + 2j), cos(3i - j)) for i in 1:6, j in 1:6] + 6I
    B = [complex(cos(i - j), sin(i * j)) for i in 1:6, j in 1:3]
    ldiv!(Engine.factorize!(copy(A), pivots), copy(B))
    bytes = only(map((f, a, p, b) -> @allocated(ldiv!(f(a, p), b)),
        (Engine.factorize!,), (copy(A),), (pivots,), (copy(B),)))
    @test bytes == 0
    # Any other element type takes the generic LU of `lu!`.
    dual = map(z -> complex(ForwardDiff.Dual(real(z), 1.0), ForwardDiff.Dual(imag(z), 0.0)), A)
    @test Engine.factorize!(copy(dual), pivots).factors == lu!(copy(dual)).factors
    uncertain = map(z -> complex(measurement(real(z), 1e-3), measurement(imag(z), 1e-3)), A)
    @test Engine.factorize!(copy(uncertain), pivots).factors == lu!(copy(uncertain)).factors
    # `factorize!` refuses an exactly singular matrix by the same exception as `lu!`.
    @test_throws SingularException(1) Engine.factorize!(zeros(ComplexF64, 3, 3), pivots)
end
