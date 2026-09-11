# Shared numerical operations for frequency-tracked eigensystems.

function _input(
        lp::LineParameters{Tc, U, PhaseDomain, Basis}
) where {Tc <: Complex, U <: Real, Basis}
    n, columns, nfrequencies = size(lp.Z.values)
    n == columns || throw(DimensionMismatch("Z must be square"))
    size(lp.Y.values) == (n, n, nfrequencies) || throw(
        DimensionMismatch("Y must be n×n×nfreq")
    )
    length(lp.f) == nfrequencies || throw(
        DimensionMismatch("f must contain one value per frequency slice")
    )
    nfrequencies > 0 || throw(ArgumentError(
        "modal transformation requires at least one frequency"
    ))
    impedance = nominal(lp.Z.values)
    admittance = nominal(lp.Y.values)
    frequencies = nominal(lp.f)
    all(isfinite, impedance) || throw(DomainError(
        impedance,
        "phase-domain impedance must be finite"
    ))
    all(isfinite, admittance) || throw(DomainError(
        admittance,
        "phase-domain admittance must be finite"
    ))
    all(isfinite, frequencies) || throw(DomainError(
        frequencies,
        "frequencies must be finite"
    ))
    return impedance, admittance, frequencies
end

@inline function _product!(destination, admittance, impedance, frequency::Integer)
    mul!(
        destination,
        @view(admittance[:, :, frequency]),
        @view(impedance[:, :, frequency])
    )
    return destination
end

function _maps(current::AbstractArray{T, 3}) where {T <: Complex}
    n, columns, nfrequencies = size(current)
    n == columns || throw(DimensionMismatch(
        "current eigenvectors must be an n×n×nfreq tensor"
    ))
    voltage = similar(current)
    inverse = similar(current)
    factor = Matrix{T}(undef, n, n)
    identity = Matrix{T}(I, n, n)
    @inbounds for frequency in 1:nfrequencies
        copyto!(
            @view(voltage[:, :, frequency]),
            transpose(@view(current[:, :, frequency]))
        )
        copyto!(factor, @view(current[:, :, frequency]))
        copyto!(@view(inverse[:, :, frequency]), identity)
        ldiv!(lu!(factor), @view(inverse[:, :, frequency]))
    end
    return ModalOperators(voltage, inverse)
end

@inline function _unit!(vector::AbstractVector)
    scale = norm(vector)
    isfinite(scale) && !iszero(scale) || return false
    vector ./= scale
    return true
end

function _orient!(vector::AbstractVector)
    _unit!(vector) || return false
    pivot = firstindex(vector)
    magnitude = abs(vector[pivot])
    @inbounds for index in Iterators.drop(eachindex(vector), 1)
        candidate = abs(vector[index])
        if candidate > magnitude
            pivot = index
            magnitude = candidate
        end
    end
    iszero(magnitude) && return false
    vector .*= conj(vector[pivot]) / magnitude
    return true
end

function _align!(vector::AbstractVector, reference::AbstractVector)
    _unit!(vector) || return false
    overlap = dot(reference, vector)
    scale = abs(overlap)
    if scale > sqrt(eps(typeof(real(scale))))
        vector .*= conj(overlap) / scale
        return true
    end
    return _orient!(vector)
end

function _transpose!(vector::AbstractVector{T}) where {T <: Complex}
    squared = zero(T)
    magnitude = zero(typeof(real(zero(T))))
    @inbounds for value in vector
        squared += value * value
        magnitude += abs2(value)
    end
    abs(squared) > sqrt(eps(typeof(magnitude))) * max(magnitude, one(magnitude)) ||
        return false
    vector ./= sqrt(squared)
    return true
end

function _seed(matrix::AbstractMatrix{T}) where {T <: Complex}
    decomposition = eigen(matrix)
    values = Vector{T}(decomposition.values)
    vectors = Matrix{T}(decomposition.vectors)
    @inbounds for mode in axes(vectors, 2)
        _orient!(@view(vectors[:, mode])) || throw(ArgumentError(
            "eigendecomposition produced a zero eigenvector"
        ))
    end
    return values, vectors
end

# Minimum-cost square assignment by the O(n³) Hungarian algorithm.
function _hungarian(cost::AbstractMatrix{R}) where {R <: Real}
    n = checksquare(cost)
    n == 0 && return Int[]
    u = zeros(R, n + 1)
    v = zeros(R, n + 1)
    matching = zeros(Int, n + 1)
    way = zeros(Int, n + 1)
    minimums = Vector{R}(undef, n + 1)
    used = falses(n + 1)

    @inbounds for row in 1:n
        matching[1] = row
        fill!(minimums, R(Inf))
        fill!(used, false)
        column = 1
        while true
            used[column] = true
            matched_row = matching[column]
            delta = R(Inf)
            next_column = 0
            for candidate in 2:(n + 1)
                used[candidate] && continue
                reduced = cost[matched_row, candidate - 1] -
                          u[matched_row + 1] - v[candidate]
                if reduced < minimums[candidate]
                    minimums[candidate] = reduced
                    way[candidate] = column
                end
                if minimums[candidate] < delta
                    delta = minimums[candidate]
                    next_column = candidate
                end
            end
            isfinite(delta) || throw(ArgumentError(
                "modal assignment cost must be finite"
            ))
            for candidate in 1:(n + 1)
                if used[candidate]
                    u[matching[candidate] + 1] += delta
                    v[candidate] -= delta
                elseif candidate > 1
                    minimums[candidate] -= delta
                end
            end
            column = next_column
            iszero(matching[column]) && break
        end
        while true
            previous = way[column]
            matching[column] = matching[previous]
            column = previous
            column == 1 && break
        end
    end

    assignment = Vector{Int}(undef, n)
    @inbounds for column in 1:n
        assignment[matching[column + 1]] = column
    end
    return assignment
end

function _match!(
        values::AbstractVector{T},
        vectors::AbstractMatrix{T},
        previous_values::AbstractVector{T},
        previous_vectors::AbstractMatrix{T}
) where {T <: Complex}
    n = length(values)
    length(previous_values) == n || throw(DimensionMismatch(
        "eigenvalue sequences must have equal length"
    ))
    size(vectors) == size(previous_vectors) == (n, n) || throw(
        DimensionMismatch("eigenvector matrices must be n×n")
    )
    R = typeof(real(zero(T)))
    cost = Matrix{R}(undef, n, n)
    @inbounds for previous in 1:n, current in 1:n
        denominator = norm(@view(previous_vectors[:, previous])) *
                      norm(@view(vectors[:, current]))
        overlap = iszero(denominator) ? zero(R) :
                  abs(dot(
            @view(previous_vectors[:, previous]),
            @view(vectors[:, current])
        )) / denominator
        cost[previous, current] = one(R) - overlap
    end
    assignment = _hungarian(cost)
    ordered_values = copy(values)
    ordered_vectors = copy(vectors)
    @inbounds for mode in 1:n
        source = assignment[mode]
        values[mode] = ordered_values[source]
        copyto!(@view(vectors[:, mode]), @view(ordered_vectors[:, source]))
    end
    return assignment
end

function _valid(
        matrix::AbstractMatrix{T},
        values::AbstractVector{T},
        vectors::AbstractMatrix{T},
        tolerance::Real
) where {T <: Complex}
    all(isfinite, values) && all(isfinite, vectors) || return false
    R = typeof(real(zero(T)))
    condition = cond(vectors)
    isfinite(condition) && condition <= inv(sqrt(eps(R))) || return false
    scale = max(norm(matrix, Inf), eps(R))
    limit = max(convert(R, tolerance), sqrt(eps(R))) * scale
    residual = Vector{T}(undef, size(matrix, 1))
    @inbounds for mode in eachindex(values)
        mul!(residual, matrix, @view(vectors[:, mode]))
        residual .-= values[mode] .* @view(vectors[:, mode])
        norm(residual, Inf) <= limit || return false
    end
    return true
end

function _fallback(
        matrix::AbstractMatrix{T},
        previous_values::AbstractVector{T},
        previous_vectors::AbstractMatrix{T}
) where {T <: Complex}
    values, vectors = _seed(matrix)
    _match!(values, vectors, previous_values, previous_vectors)
    @inbounds for mode in eachindex(values)
        _align!(@view(vectors[:, mode]), @view(previous_vectors[:, mode]))
    end
    return values, vectors
end
