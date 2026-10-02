# Matrix reductions shared by the line-parameter backends. Primitive matrices are
# ordered by terminal; a `ReductionPlan` fixes the index operations once.

"""
$(TYPEDSIGNATURES)

Replace each cyclic diagonal of the square `matrix` by its mean, in place, and
return `matrix`. The result is the ideally transposed matrix: every phase occupies
every position for an equal share of the line length.
"""
function ideal_transposition!(matrix::AbstractMatrix)
    n = checksquare(matrix)
    coefficients = similar(diag(matrix))
    @inbounds for offset in 0:(n - 1)
        total = zero(eltype(matrix))
        for row in 1:n
            total += matrix[row, 1 + mod(row - 1 + offset, n)]
        end
        coefficients[offset + 1] = total / n
    end
    @inbounds for row in 1:n, column in 1:n

        matrix[row, column] = coefficients[mod1(column - row + 1, n)]
    end
    return matrix
end

"""
$(TYPEDSIGNATURES)

Return the terminal permutation that places the first conductor of each active
phase in encounter order, then the remaining conductors of each phase in the same
phase order, then the conductors whose `map` entry is zero.
"""
function reorder_indices(map::AbstractVector{<:Integer})
    n = length(map)
    phases = Int[]                     # encounter order of active phase IDs
    firsts = Int[]
    sizehint!(firsts, n)
    eliminated = Int[]                 # phase-zero conductors
    sizehint!(eliminated, n)
    tails = Dict{Int, Vector{Int}}()   # phase => remaining indices

    seen = Set{Int}()
    @inbounds for (i, p) in pairs(map)
        if p > 0
            if !(p in seen)
                push!(seen, p)
                push!(phases, p)
                push!(firsts, i)
            else
                push!(get!(tails, p, Int[]), i)
            end
        else
            push!(eliminated, i)
        end
    end

    perm = Vector{Int}(undef, n)
    k = 1
    @inbounds begin
        for i in firsts
            perm[k] = i
            k += 1
        end
        for p in phases
            if haskey(tails, p)
                for i in tails[p]
                    perm[k] = i
                    k += 1
                end
            end
        end
        for i in eliminated
            perm[k] = i
            k += 1
        end
    end
    return perm
end

"""
$(TYPEDSIGNATURES)

Eliminate matrix rows and columns whose `phase_map` entry is zero.

For retained indices `1` and eliminated indices `2`, calculate the Schur
complement

```math
M_{\\mathrm{red}} = M_{11} - M_{12}M_{22}^{-1}M_{21}.
```

# Arguments

- `M`: square complex matrix.
- `phase_map`: active-phase assignment for each row and column. Nonzero IDs
  identify active phases. Zero marks a grounded or eliminated conductor.

# Returns

- The reduced matrix.
"""
function kron_reduce(
        M::Matrix{Complex{T}},
        phase_map::Vector{Int}
) where {T <: Real}
    retained = count(!=(0), phase_map)
    reduced = similar(M, retained, retained)
    kron_reduce!(M, phase_map, reduced)
    return reduced
end

"""
$(TYPEDSIGNATURES)

Write the Kron-reduced matrix from [`kron_reduce`](@ref) into `Mred`.

# Arguments

- `M`: square complex matrix.
- `phase_map`: active-phase assignment for each row and column. Nonzero IDs
  identify active phases. Zero marks a grounded or eliminated conductor.
- `Mred`: destination matrix.

# Returns

- `nothing`.
"""
function kron_reduce!(
        M::Matrix{Complex{T}},
        phase_map::Vector{Int},
        Mred::Matrix{Complex{T}}
) where {T <: Real}
    checksquare(M) == length(phase_map) || throw(DimensionMismatch(
        "phase map must contain one entry per matrix row"))
    keep = findall(!=(0), phase_map)
    eliminate = findall(==(0), phase_map)
    # This entry point also permits Mred to alias M. The reusable-buffer entry
    # point below requires independent scratch, as in the computation workspace.
    source = Base.unalias(Mred, M)
    kron_reduce!(source, keep, eliminate, Mred,
        similar(M, length(eliminate), length(eliminate)),
        similar(M, length(keep), length(eliminate)),
        similar(M, length(eliminate), length(keep)))
    return nothing
end

"""
$(TYPEDSIGNATURES)

Write the Schur complement of `matrix` that retains the indices `keep` and
eliminates the indices `eliminate` into `reduced`, and return `reduced`.

`factor`, `coupling` and `right_hand_side` are scratch blocks of sizes
`length(eliminate)²`, `length(keep)×length(eliminate)` and
`length(eliminate)×length(keep)`. The destination and the scratch blocks must not
alias `matrix` or one another. Without eliminated indices, `reduced` receives
`matrix[keep, keep]`.
"""
function kron_reduce!(
        matrix::AbstractMatrix{Complex{T}},
        keep::AbstractVector{Int},
        eliminate::AbstractVector{Int},
        reduced::AbstractMatrix{Complex{T}},
        factor::AbstractMatrix{Complex{T}},
        coupling::AbstractMatrix{Complex{T}},
        right_hand_side::AbstractMatrix{Complex{T}}
) where {T <: Real}
    retained = length(keep)
    removed = length(eliminate)
    size(reduced) == (retained, retained) || throw(DimensionMismatch(
        "reduced matrix storage must be $retained×$retained"
    ))
    size(factor) == (removed, removed) || throw(DimensionMismatch(
        "Kron factor storage must be $removed×$removed"
    ))
    size(coupling) == (retained, removed) || throw(DimensionMismatch(
        "Kron coupling storage must be $retained×$removed"
    ))
    size(right_hand_side) == (removed, retained) || throw(DimensionMismatch(
        "Kron right-hand-side storage must be $removed×$retained"
    ))
    @inbounds for column in 1:retained, row in 1:retained

        reduced[row, column] = matrix[keep[row], keep[column]]
    end
    iszero(removed) && return reduced
    @inbounds for column in 1:removed, row in 1:removed

        factor[row, column] = matrix[eliminate[row], eliminate[column]]
    end
    @inbounds for column in 1:removed, row in 1:retained

        coupling[row, column] = matrix[keep[row], eliminate[column]]
    end
    @inbounds for column in 1:retained, row in 1:removed

        right_hand_side[row, column] = matrix[eliminate[row], keep[column]]
    end
    factorization = lu!(factor)
    ldiv!(factorization, right_hand_side)
    mul!(
        reduced,
        coupling,
        right_hand_side,
        -one(Complex{T}),
        one(Complex{T})
    )
    return reduced
end

"""
$(TYPEDSIGNATURES)

Return the pairs `(first, duplicate)` that merge each further conductor of an
active phase into the first conductor of that phase, in the order of `phases`.
Zero entries take part in no pair.
"""
function bundle_operations(phases::AbstractVector{<:Integer})
    first_index = Dict{Int, Int}()
    operations = Tuple{Int, Int}[]
    @inbounds for (index, phase) in pairs(phases)
        phase > 0 || continue
        first = get(first_index, phase, 0)
        if iszero(first)
            first_index[phase] = index
        else
            push!(operations, (first, index))
        end
    end
    return operations
end

"""
$(TYPEDSIGNATURES)

Apply the bundle change of basis of the pairs from [`bundle_operations`](@ref) to
`matrix` in place: each duplicate column, then each duplicate row, loses its first
conductor's column or row. Return `matrix`.
"""
function merge_bundles!(
        matrix::AbstractMatrix{T},
        operations::AbstractVector{<:Tuple{Int, Int}}
) where {T}
    @inbounds for (first, duplicate) in operations
        base_column = @view matrix[:, first]
        column = @view matrix[:, duplicate]
        if matrix isa StridedMatrix{T} &&
           T <: Union{Float32, Float64, ComplexF32, ComplexF64}
            axpy!(-one(T), base_column, column)
        else
            column .-= base_column
        end
    end
    @inbounds for (first, duplicate) in operations
        base_row = @view matrix[first, :]
        row = @view matrix[duplicate, :]
        if matrix isa StridedMatrix{T} &&
           T <: Union{Float32, Float64, ComplexF32, ComplexF64}
            axpy!(-one(T), base_row, row)
        else
            row .-= base_row
        end
    end
    return matrix
end

"""
$(TYPEDSIGNATURES)

Apply the bundle change of basis to conductors assigned to the same active
phase. A zero assignment denotes an independent conductor selected for
grounded or eliminated-conductor reduction and is never interpreted as a bundle
identity.
"""
function merge_bundles!(
        matrix::AbstractMatrix{T},
        phases::AbstractVector{<:Integer}
) where {T}
    n = size(matrix, 1)
    size(matrix, 2) == n == length(phases) ||
        throw(ArgumentError("shape mismatch"))
    operations = bundle_operations(phases)
    merge_bundles!(matrix, operations)
    reduced = copy(phases)
    @inbounds for (_, duplicate) in operations
        reduced[duplicate] = 0
    end
    return matrix, reduced
end

"""
$(TYPEDEF)

Fix the index operations that reduce primitive matrices, ordered by terminal, to the
matrices of the retained phases: terminal reorder, bundle change of basis, Kron
elimination and ideal transposition.

$(TYPEDFIELDS)
"""
struct ReductionPlan
    "Reordering of the primitive terminals, from [`reorder_indices`](@ref)."
    permutation::Vector{Int}
    "Bundle change-of-basis pairs in reordered indices, empty without bundle reduction."
    bundles::Vector{Tuple{Int, Int}}
    "Reordered indices retained by the Kron elimination."
    keep::Vector{Int}
    "Reordered indices eliminated by the Kron elimination."
    eliminate::Vector{Int}
    "Whether the retained matrices are ideally transposed."
    transposition::Bool
    "Phase assignment of each retained row."
    phase_map::Vector{Int}
    "Primitive terminal of each retained row."
    indices::Vector{Int}
end

"""
$(TYPEDSIGNATURES)

Build the reduction plan of `phase_map`, the active-phase assignment of each
primitive terminal. Nonzero IDs identify active phases, and zero marks a grounded
or eliminated conductor.

# Keywords

- `reduce_bundle`: merge the conductors of each active phase into one.
- `kron_reduction`: eliminate the conductors with phase zero.
- `ideal_transposition`: average the retained matrices over cyclic transposition.
"""
function ReductionPlan(phase_map::AbstractVector{<:Integer}; reduce_bundle::Bool,
        kron_reduction::Bool, ideal_transposition::Bool)
    permutation = reorder_indices(phase_map)
    ordered = phase_map[permutation]
    reduced = copy(ordered)
    seen = Set{Int}()
    @inbounds for (index, phase) in pairs(ordered)
        if phase > 0 && phase in seen
            reduced[index] = 0
        elseif phase > 0
            push!(seen, phase)
        end
    end
    # Without Kron reduction, bundle reduction still eliminates the merged
    # duplicates and keeps the conductors with phase zero, marked -1.
    retained = if reduce_bundle
        kron_reduction ? reduced :
        [ordered[index] == 0 ? -1 : reduced[index] for index in eachindex(reduced)]
    else
        kron_reduction ? ordered : nothing
    end
    bundles = reduce_bundle ? bundle_operations(ordered) : Tuple{Int, Int}[]
    retained === nothing && return ReductionPlan(permutation, bundles,
        collect(eachindex(ordered)), Int[], ideal_transposition, ordered, copy(permutation))
    keep = findall(!=(0), retained)
    return ReductionPlan(permutation, bundles, keep, findall(==(0), retained),
        ideal_transposition, retained[keep], permutation[keep])
end

"""
$(TYPEDEF)

Hold the storage of [`reduce_line_matrices!`](@ref) for one [`ReductionPlan`](@ref)
and one element type.

$(TYPEDFIELDS)
"""
struct ReductionBuffers{T <: Number}
    "Reordered primitive matrix, reused for `Z` and `P`."
    ordered::Matrix{T}
    "Retained potential coefficients, factorized in place."
    potential::Matrix{T}
    "Kron elimination block of the eliminated indices."
    factor::Matrix{T}
    "Kron elimination block coupling retained to eliminated indices."
    coupling::Matrix{T}
    "Kron elimination right-hand side."
    right_hand_side::Matrix{T}
    "Identity of the retained size."
    identity::Matrix{T}
end

"""
$(TYPEDSIGNATURES)

Allocate the storage of [`reduce_line_matrices!`](@ref) for `plan` with elements of
type `T`.
"""
function ReductionBuffers{T}(plan::ReductionPlan) where {T <: Number}
    n, retained, removed = length(plan.permutation), length(plan.keep), length(plan.eliminate)
    return ReductionBuffers{T}(Matrix{T}(undef, n, n), Matrix{T}(undef, retained, retained),
        Matrix{T}(undef, removed, removed), Matrix{T}(undef, retained, removed),
        Matrix{T}(undef, removed, retained), Matrix{T}(I, retained, retained))
end

"""
$(TYPEDSIGNATURES)

Reduce one frequency of primitive line matrices to the retained phases of `plan`.

`Zprimitive` \\[Ω/m\\] and `Pprimitive` \\[m/F\\] are ordered by terminal. Both are
reordered, merged by bundle and Kron-reduced. With ideal transposition, the
retained `Z` and `P` are averaged. `Z` receives the retained series impedance and
`Y` receives the shunt admittance ``Y = s P^{-1}`` \\[S/m\\] of the retained `P`,
with `s = jω`. `buffers` comes from [`ReductionBuffers`](@ref) for `plan`.

With `diagnostics = Val(true)`, nonfinite `P` or `Y` throws `ArgumentError`, and the
return value is the infinity-norm residual of ``P P^{-1} - I`` and the 2-norm
condition number of `P`. A nonfinite condition estimate or a residual above
``\\max(\\sqrt{ε}, 32nε\\max(1, κ))`` produces a warning. Otherwise the return value
is `nothing`.
"""
function reduce_line_matrices!(
        Z::AbstractMatrix{T},
        Y::AbstractMatrix{T},
        Zprimitive::AbstractMatrix,
        Pprimitive::AbstractMatrix,
        s::Number,
        plan::ReductionPlan,
        buffers::ReductionBuffers{T},
        diagnostics::Union{Val{true}, Val{false}} = Val(false)
) where {T <: Number}
    function retain!(retained, primitive)
        ordered = buffers.ordered
        @inbounds for column in eachindex(plan.permutation), row in eachindex(plan.permutation)
            ordered[row, column] = primitive[plan.permutation[row], plan.permutation[column]]
        end
        merge_bundles!(ordered, plan.bundles)
        kron_reduce!(ordered, plan.keep, plan.eliminate, retained,
            buffers.factor, buffers.coupling, buffers.right_hand_side)
        plan.transposition && ideal_transposition!(retained)
        return retained
    end
    retain!(Z, Zprimitive)
    P = retain!(buffers.potential, Pprimitive)
    if diagnostics === Val(false)
        ldiv!(Y, lu!(P), buffers.identity)
        Y .*= s
        return nothing
    end
    all(isfinite, P) || throw(ArgumentError(
        "retained potential coefficients contain nonfinite values"))
    R = real(T)
    condition_number = convert(R, cond(P))
    isfinite(condition_number) ||
        @warn "Potential-coefficient condition estimate is not finite" condition_number
    # The reordered storage is free after the Kron elimination of P.
    coefficients = copyto!(view(buffers.ordered, axes(P)...), P)
    ldiv!(Y, lu!(P), buffers.identity)
    all(isfinite, Y) || throw(ArgumentError("computed admittance contains nonfinite values"))
    residual = convert(R, norm(coefficients * Y - buffers.identity, Inf))
    tolerance = max(sqrt(eps(R)), convert(R, 32size(Y, 1) * eps(R) * max(one(R), condition_number)))
    isfinite(residual) && residual <= tolerance ||
        @warn "Potential-coefficient inversion residual target was not met" residual tolerance condition_number
    Y .*= s
    return (; residual, condition_number)
end
