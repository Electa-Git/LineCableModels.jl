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
zero inside a bundle, and the Kron elimination of ``j`` (`eliminate!`) completes the merge.

# Arguments

- `matrix`: square matrix of impedances \\[Ω/m\\] or potential coefficients \\[m/F\\], one
  row per conductor.
- `pairs`: the rows `(first, other)` of each conductor that a bundle merges into the first
  conductor of its index.

# Returns

- `matrix`, with the merged columns and rows.
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

Write into `reduced` the Kron elimination of `matrix`, the Schur complement of the block of
its eliminated rows and columns:

```math
M_{\\mathrm{r}} = M_{kk} - M_{ke} M_{ee}^{-1} M_{ek},
```

where ``k`` indexes the kept rows and columns in the order of `kept`, and ``e`` the
eliminated ones in the order of `eliminated`. An eliminated conductor is at zero voltage:
grounded, or merged into a bundle by `bundle!`. The LU decomposition of ``M_{ee}`` solves
for ``M_{ee}^{-1} M_{ek}``. Without an eliminated row, ``M_{\\mathrm{r}} = M_{kk}``.

# Arguments

- `reduced`: matrix that receives the result, of `length(kept)` rows and columns.
- `matrix`: square matrix of impedances \\[Ω/m\\] or potential coefficients \\[m/F\\], one
  row per conductor.
- `kept`: rows of `matrix` that the elimination keeps.
- `eliminated`: rows of `matrix` that the elimination removes.
- `factor`: matrix of `length(eliminated)` rows and columns that receives the LU factors of
  ``M_{ee}``.
- `coupling`: matrix of `length(kept)` rows and `length(eliminated)` columns that receives
  ``M_{ke}``.
- `right_hand_side`: matrix of `length(eliminated)` rows and `length(kept)` columns that
  receives ``M_{ee}^{-1} M_{ek}``.

# Returns

- `reduced`, in the unit of `matrix`.
"""
function eliminate!(
        reduced::AbstractMatrix{Complex{T}}, matrix::AbstractMatrix{Complex{T}},
        kept::AbstractVector{Int}, eliminated::AbstractVector{Int},
        factor::AbstractMatrix{Complex{T}}, coupling::AbstractMatrix{Complex{T}},
        right_hand_side::AbstractMatrix{Complex{T}}) where {T <: Real}
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

        right_hand_side[row, column] = matrix[eliminated[row], kept[column]]
    end
    ldiv!(lu!(factor), right_hand_side)
    return mul!(reduced, coupling, right_hand_side, -one(Complex{T}), one(Complex{T}))
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

- `matrix`, circulant.

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
