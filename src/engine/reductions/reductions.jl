"""
$(TYPEDSIGNATURES)

Return the plan of the bundle and Kron reductions of a system under the formulation's
options. The plan lists the conductor indices of the reductions and the labels of the
conductors and of the kept rows. The rows of the reduced matrices start with the first
conductor of each conductor index `k ≥ 1`, in the order of the system. The other conductors
of each index follow, and the conductors of index zero come last. With `bundle`, the
reductions merge each of the other conductors into the first of its index and remove its
row. With `kron`, they remove the rows of index zero.

# Arguments

- `system`: cable designs, placements, and conductor connections.
- `options`: formulation options that select the bundle and Kron reductions.

# Returns

- `label`: label of each conductor, `cable:<i>:<terminal>`.
- `bundles`: conductor pairs `(first, other)` to merge, empty when `bundle` is false.
- `kept`: conductor indices to keep, in result order.
- `eliminated`: conductor indices to remove.
- `terminals`: labels of the kept rows, with a bundle's conductor labels joined by `+`.
"""
function reductions(system::CableSystem, options::FormulationOptions)
    connection = [row.connection for row in system.terminals]
    label = ["cable:$(row.cable):$(system.designs[row.cable].terminals[row.terminal].name)"
             for row in system.terminals]
    firsts = [i
              for i in eachindex(connection)
              if connection[i] > 0 && connection[i] ∉ view(connection, 1:(i - 1))]
    others = [i
              for j in firsts
              for i in eachindex(connection)
              if i > j && connection[i] == connection[j]]
    grounded = [i for i in eachindex(connection) if iszero(connection[i])]
    permutation = [firsts; others; grounded]
    merged = length(firsts) .+ eachindex(others)
    zeroed = (length(firsts) + length(others)) .+ eachindex(grounded)
    bundles = [(k, row)
               for (row, i) in zip(merged, others)
               for k in eachindex(firsts) if connection[firsts[k]] == connection[i]]
    eliminated = [options.bundle ? merged : Int[]; options.kron ? zeroed : Int[]]
    kept = setdiff(eachindex(permutation), eliminated)
    terminals = [options.bundle && connection[i] > 0 ?
                 join(label[connection .== connection[i]], "+") : label[i]
                 for i in permutation[kept]]
    bundles = [(permutation[first], permutation[other]) for (first, other) in bundles]
    kept = permutation[kept]
    eliminated = permutation[eliminated]
    return (; label, bundles = options.bundle ? bundles : Tuple{Int, Int}[], kept,
        eliminated, terminals)
end

"""
$(TYPEDSIGNATURES)

Merge each other conductor of a bundle into the first conductor of its index, in the rows
and columns of `matrix`. The conductors of a bundle share one voltage and carry its current
together. For the pairs ``(i, j)`` of `pairs`, with ``i`` the first conductor of a bundle
and ``j`` another one, the change of basis is

```math
M' = T M T^{\\mathsf{T}}, \\qquad T = I - \\sum_{(i, j)} e_j e_i^{\\mathsf{T}},
```

where ``e_i`` is the unit vector of row ``i``. It subtracts column ``i`` from column ``j``,
then row ``i`` from row ``j``. Row ``j`` then gives the voltage of ``j`` less that of ``i``,
zero inside a bundle. The Kron elimination of ``j`` by [`eliminate!`](@ref) completes the
merge.

# Arguments

- `matrix`: square matrix of impedances \\[Ω/m\\] or potential coefficients \\[m/F\\], one
  row per conductor.
- `pairs`: the rows `(first, other)` of each conductor that a bundle merges into the first
  conductor of its index.

# Returns

- `matrix`, with transformed columns and rows, in its original unit \\[Ω/m or m/F\\].
"""
function bundle!(matrix::AbstractMatrix, pairs::AbstractVector{Tuple{Int, Int}})
    for (first, other) in pairs
        for row in axes(matrix, 1)
            matrix[row, other] -= matrix[row, first]
        end
    end
    for (first, other) in pairs
        for column in axes(matrix, 2)
            matrix[other, column] -= matrix[first, column]
        end
    end
    return matrix
end

"""
$(TYPEDSIGNATURES)

Factor the square matrix `A` in place into ``P A = L U``. ``L`` is unit lower triangular and
``U`` upper triangular. ``P`` holds the row interchanges that move the largest entry of each
column to the diagonal. For a BLAS element type, LAPACK's `getrf!` writes the interchanges
into the first `size(A, 1)` entries of `pivots`, and the factorization does not allocate.
Any other element type takes Julia's generic LU of `lu!`, which allocates its own
interchanges.

# Arguments

- `A`: square matrix, overwritten by its factors.
- `pivots`: integer vector of at least `size(A, 1)` entries that receives the row
  interchanges.

# Returns

- The `LU` factorization, its factors in `A`.

# Errors

- Throws `SingularException` for an exactly singular matrix, as `lu!` does.
"""
function factorize!(
        A::StridedMatrix{T}, pivots::AbstractVector{BlasInt}) where {T <: BlasFloat}
    ipiv = view(pivots, 1:checksquare(A))
    _, _, info = LAPACK.getrf!(A, ipiv)
    info > 0 && throw(SingularException(info))
    return LU{T, typeof(A), typeof(ipiv)}(A, ipiv, info)
end

factorize!(A::AbstractMatrix, ::AbstractVector{BlasInt}) = lu!(A)

"""
$(TYPEDSIGNATURES)

Write into `reduced` the Kron elimination of `matrix`, the Schur complement of the block of
eliminated rows and columns, retained on the kept rows and columns:

```math
M_{\\mathrm{r}} = M_{kk} - M_{ke} M_{ee}^{-1} M_{ek},
```

where ``k`` indexes the kept rows and columns in the order of `kept`, and ``e`` the
eliminated ones in the order of `eliminated`. An eliminated conductor is at zero voltage:
grounded, or merged into a bundle by [`bundle!`](@ref). The LU decomposition of ``M_{ee}``
by [`factorize!`](@ref) solves for ``M_{ee}^{-1} M_{ek}``. Without an eliminated row,
``M_{\\mathrm{r}} = M_{kk}``.

# Arguments

- `reduced`: destination matrix with `length(kept)` rows and columns \\[Ω/m or m/F\\].
- `matrix`: square matrix of impedances \\[Ω/m\\] or potential coefficients \\[m/F\\], one
  row per conductor.
- `kept`: rows of `matrix` that the elimination keeps.
- `eliminated`: rows of `matrix` that the elimination removes.
- `factor`: array for the LU factors of ``M_{ee}``, square with `length(eliminated)` rows.
  Before factorization its unit is \\[Ω/m or m/F\\].
- `coupling`: matrix of `length(kept)` rows and `length(eliminated)` columns that receives
  ``M_{ke}`` \\[Ω/m or m/F\\].
- `rhs`: matrix of `length(eliminated)` rows and `length(kept)` columns that receives
  ``M_{ee}^{-1} M_{ek}`` \\[dimensionless\\].
- `pivots`: integer vector of at least `length(eliminated)` entries that receives the row
  interchanges of the LU decomposition.

# Returns

- `reduced`, in the unit of `matrix` \\[Ω/m or m/F\\].
"""
function eliminate!(
        reduced::AbstractMatrix{Complex{T}}, matrix::AbstractMatrix{Complex{T}},
        kept::AbstractVector{Int}, eliminated::AbstractVector{Int},
        factor::AbstractMatrix{Complex{T}}, coupling::AbstractMatrix{Complex{T}},
        rhs::AbstractMatrix{Complex{T}}, pivots::AbstractVector{BlasInt}) where {T <: Real}
    retained = length(kept)
    removed = length(eliminated)
    for column in 1:retained, row in 1:retained

        reduced[row, column] = matrix[kept[row], kept[column]]
    end
    isempty(eliminated) && return reduced
    for column in 1:removed, row in 1:removed

        factor[row, column] = matrix[eliminated[row], eliminated[column]]
    end
    for column in 1:removed, row in 1:retained

        coupling[row, column] = matrix[kept[row], eliminated[column]]
    end
    for column in 1:retained, row in 1:removed

        rhs[row, column] = matrix[eliminated[row], kept[column]]
    end
    ldiv!(factorize!(factor, pivots), rhs)
    return mul!(reduced, coupling, rhs, -one(Complex{T}), one(Complex{T}))
end

"""
$(TYPEDSIGNATURES)

Allocate the arrays of the reductions of a computation: the Kron elimination by
[`eliminate!`](@ref) and the inversion of the reduced potential coefficients after it. With
`m = length(kept)` and `e = length(eliminated)`, `reduced` has `m` rows and columns,
`factor` `e`, `coupling` `m` rows and `e` columns, `rhs` `e` rows and `m` columns, and
`identity` is the identity matrix of order `m`, each of `Complex{T}`. The integer vector
`pivots` of `max(m, e)` entries serves each LU decomposition of a frequency, one after the
other.

# Arguments

- `eliminate!`: the reduction whose arrays these are.
- `T`: real number type of the computation.
- `kept`: conductors that the reductions keep.
- `eliminated`: conductors that the Kron elimination removes.

# Returns

- A named tuple of `reduced`, `factor`, `coupling`, `rhs`, `identity` and `pivots`.
"""
function initialize_buffers(::typeof(eliminate!), ::Type{T}, kept::AbstractVector{Int},
        eliminated::AbstractVector{Int}) where {T <: Real}
    m = length(kept)
    e = length(eliminated)
    return (reduced = Matrix{Complex{T}}(undef, m, m),
        factor = Matrix{Complex{T}}(undef, e, e),
        coupling = Matrix{Complex{T}}(undef, m, e),
        rhs = Matrix{Complex{T}}(undef, e, m),
        identity = Matrix{Complex{T}}(I, m, m),
        pivots = Vector{BlasInt}(undef, max(m, e)))
end

"""
$(TYPEDSIGNATURES)

Replace `matrix` by its ideal transposition, the mean over the ``n`` cyclic permutations of
its rows and columns:

```math
\\bar{M}_{ij} = \\frac{1}{n} \\sum_{k = 0}^{n - 1} M_{i + k, j + k},
```

with the indices taken cyclically, modulo ``n``. Each cyclic diagonal, the entries
``(i, j)`` of one value of ``j - i`` modulo ``n``, takes the mean of its entries, so the
result is a circulant matrix.

# Arguments

- `matrix`: square matrix of impedances \\[Ω/m\\] or admittances \\[S/m\\], one row per
  conductor.

# Returns

- `matrix`, with circulant entries in its original unit \\[Ω/m or S/m\\].

# Errors

- Throws `DimensionMismatch` when `matrix` is not square.
"""
function transposition!(matrix::AbstractMatrix)
    n = checksquare(matrix)
    for offset in 0:(n - 1)
        total = zero(eltype(matrix))
        for row in 1:n
            total += matrix[row, 1 + mod(row - 1 + offset, n)]
        end
        for row in 1:n
            matrix[row, 1 + mod(row - 1 + offset, n)] = total / n
        end
    end
    return matrix
end
