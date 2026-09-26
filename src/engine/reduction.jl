using LinearAlgebra: axpy!

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

- `M`: Square complex matrix.
- `phase_map`: Active-phase assignment for each row and column. Nonzero IDs
  identify active phases; zero marks a grounded/eliminated conductor.

# Returns

- The reduced matrix.
"""
function kronify(
        M::Matrix{Complex{T}},
        phase_map::Vector{Int}
) where {T <: Real}
    retained = count(!=(0), phase_map)
    reduced = similar(M, retained, retained)
    kronify!(M, phase_map, reduced)
    return reduced
end

"""
$(TYPEDSIGNATURES)

Write the Kron-reduced matrix from [`kronify`](@ref) into `Mred`.

# Arguments

- `M`: Square complex matrix.
- `phase_map`: Active-phase assignment for each row and column. Nonzero IDs
  identify active phases; zero marks a grounded/eliminated conductor.
- `Mred`: Destination matrix.

# Returns

- `nothing`.
"""
function kronify!(
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
    kronify!(source, keep, eliminate, Mred,
        similar(M, length(eliminate), length(eliminate)),
        similar(M, length(keep), length(eliminate)),
        similar(M, length(eliminate), length(keep)))
    return nothing
end

# The destination and three scratch blocks must not alias the source or one
# another. Callers with aliased storage use the phase-map entry point above.
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

Apply the bundle change of basis to conductors assigned to the same active
phase. A zero assignment denotes an independent conductor selected for
grounded/eliminated-conductor reduction and is never interpreted as a bundle
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

function _structural_reduction(
        matrix::Matrix{Complex{T}},
        permutation::Vector{Int},
        bundle_pairs::Vector{Tuple{Int, Int}},
        kron_map,
        options::FormulationOptions
) where {T <: Real}
    transformed = matrix[permutation, permutation]
    options.data.reduce_bundle && merge_bundles!(transformed, bundle_pairs)
    if kron_map !== nothing
        eliminate = findall(==(0), kron_map)
        if !isempty(eliminate)
            transformed = kronify!(
                transformed,
                findall(!=(0), kron_map),
                eliminate,
                Matrix{Complex{T}}(
                    undef,
                    count(!=(0), kron_map),
                    count(!=(0), kron_map)
                ),
                Matrix{Complex{T}}(
                    undef,
                    length(eliminate),
                    length(eliminate)
                ),
                Matrix{Complex{T}}(
                    undef,
                    count(!=(0), kron_map),
                    length(eliminate)
                ),
                Matrix{Complex{T}}(
                    undef,
                    length(eliminate),
                    count(!=(0), kron_map)
                )
            )
        end
    end
    options.data.ideal_transposition && ideal_transposition!(transformed)
    return transformed
end

"""
$(TYPEDSIGNATURES)

Apply the engine-owned terminal reorder, bundle merge, Kron elimination, and
ideal-transposition operations to primitive impedance and potential-coefficient
frequency scans.

# Arguments

- `Z_primitive`: Primitive series impedance \\[Ω/m\\], ordered by terminal.
- `P_primitive`: Primitive inverse-admittance coefficient \\[m/S\\], ordered
  by terminal. Its inverse is shunt admittance in S/m.
- `phase_map`: Active-phase assignment aligned with primitive terminal order.
  Nonzero IDs identify retained phases; zero identifies a grounded/eliminated
  conductor.
- `options`: Normalized shared line-parameter formulation options.

# Returns

- A named tuple containing reduced `Z` \\[Ω/m\\], reduced `P` \\[m/S\\], and the
  retained active-phase assignment.
"""
function reduce_primitive_matrices(
        Z_primitive::Array{Complex{T}, 3},
        P_primitive::Array{Complex{T}, 3},
        phase_map::Vector{Int},
        options::FormulationOptions
) where {T <: Real}
    size(Z_primitive) == size(P_primitive) || throw(DimensionMismatch(
        "primitive Z and P scans must have identical dimensions",
    ))
    n = size(Z_primitive, 1)
    size(Z_primitive, 2) == n == length(phase_map) || throw(DimensionMismatch(
        "primitive matrices and phase_map must describe the same terminals",
    ))
    all(name -> haskey(options.data, name),
        (:reduce_bundle, :kron_reduction, :ideal_transposition)) ||
        throw(ArgumentError("shared line-parameter options are incomplete"))

    permutation, reordered, kron_map = _reduction_map(
        phase_map, (options = options,)
    )
    bundle_pairs = bundle_operations(reordered)
    retained_map = kron_map === nothing ? reordered : kron_map[findall(!=(0), kron_map)]
    output_size = length(retained_map)
    frequency_count = size(Z_primitive, 3)
    Z = Array{Complex{T}, 3}(undef, output_size, output_size, frequency_count)
    P = similar(Z)
    for frequency in 1:frequency_count
        @views Z[:, :, frequency] .= _structural_reduction(
            Matrix(Z_primitive[:, :, frequency]),
            permutation,
            bundle_pairs,
            kron_map,
            options
        )
        @views P[:, :, frequency] .= _structural_reduction(
            Matrix(P_primitive[:, :, frequency]),
            permutation,
            bundle_pairs,
            kron_map,
            options
        )
    end
    return (; Z, P, phase_map = retained_map,
        indices = kron_map === nothing ? permutation : permutation[findall(!=(0), kron_map)])
end
