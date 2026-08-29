using LinearAlgebra: axpy!

function reorder_indices(map::AbstractVector{<:Integer})
    n = length(map)
    phases = Int[]                     # encounter order of phases > 0
    firsts = Int[]
    sizehint!(firsts, n)
    zeros = Int[]
    sizehint!(zeros, n)
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
            push!(zeros, i)
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
        for i in zeros
            perm[k] = i
            k += 1
        end
    end
    return perm
end

# Non-mutating reorder (2D)
function reorder_M(M::AbstractMatrix, map::AbstractVector{<:Integer})
    n = size(M, 1)
    n == size(M, 2) == length(map) || throw(ArgumentError("shape mismatch"))
    perm = reorder_indices(map)
    return M[perm, perm], map[perm]
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

- `M`: Square complex matrix.
- `phase_map`: Phase assignment for each row and column. Zero marks an
  eliminated conductor.

# Returns

- The reduced matrix.
"""
function kronify(
        M::Matrix{Complex{T}},
        phase_map::Vector{Int}
) where {T <: Real}
    keep = findall(!=(0), phase_map)
    eliminate = findall(==(0), phase_map)

    M11 = M[keep, keep]
    M12 = M[keep, eliminate]
    M21 = M[eliminate, keep]
    M22 = M[eliminate, eliminate]

    return M11 - (M12 * inv(M22)) * M21
end

"""
$(TYPEDSIGNATURES)

Write the Kron-reduced matrix from [`kronify`](@ref) into `Mred`.

# Arguments

- `M`: Square complex matrix.
- `phase_map`: Phase assignment for each row and column.
- `Mred`: Destination matrix.

# Returns

- `nothing`.
"""
function kronify!(
        M::Matrix{Complex{T}},
        phase_map::Vector{Int},
        Mred::Matrix{Complex{T}}
) where {T <: Real}
    keep = findall(!=(0), phase_map)
    eliminate = findall(==(0), phase_map)

    M11 = M[keep, keep]
    M12 = M[keep, eliminate]
    M21 = M[eliminate, keep]
    M22 = M[eliminate, eliminate]
    @views @inbounds Mred .= M11 - (M12 * inv(M22)) * M21
    return nothing
end

function kronify!(
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

Apply the bundle change of basis to conductors assigned to the same positive
phase. Phase zero denotes independent conductors selected for elimination and
is never interpreted as a bundle identity.
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
