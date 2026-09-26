# Shared numerical operations for frequency-tracked eigensystems.

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

function normalize_bilinear!(vector::AbstractVector{T}) where {T <: Complex}
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
function _assignment_workspace(::Type{T},n) where {T<:Complex}
    R=typeof(real(zero(T)))
    return (cost=Matrix{R}(undef,n,n),u=zeros(R,n+1),v=zeros(R,n+1),
        matching=zeros(Int,n+1),way=zeros(Int,n+1),minimums=Vector{R}(undef,n+1),
        used=falses(n+1),assignment=Vector{Int}(undef,n),
        ordered_values=Vector{T}(undef,n),ordered_vectors=Matrix{T}(undef,n,n),
        residual=Vector{T}(undef,n))
end

function hungarian_assignment!(cost::AbstractMatrix{R},work) where {R <: Real}
    n = checksquare(cost)
    n == 0 && return Int[]
    u,v,matching,way,minimums,used=work.u,work.v,work.matching,
        work.way,work.minimums,work.used
    fill!(u,zero(R));fill!(v,zero(R))
    fill!(matching,0);fill!(way,0)

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

    assignment = work.assignment
    @inbounds for column in 1:n
        assignment[matching[column + 1]] = column
    end
    return assignment
end
function _match!(
        values::AbstractVector{T},
        vectors::AbstractMatrix{T},
        previous_values::AbstractVector{T},
        previous_vectors::AbstractMatrix{T},work
) where {T <: Complex}
    n = length(values)
    length(previous_values) == n || throw(DimensionMismatch(
        "eigenvalue sequences must have equal length"
    ))
    size(vectors) == size(previous_vectors) == (n, n) || throw(
        DimensionMismatch("eigenvector matrices must be n×n")
    )
    R = typeof(real(zero(T)))
    cost = work.cost
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
    assignment = hungarian_assignment!(cost,work)
    ordered_values = work.ordered_values
    ordered_vectors = work.ordered_vectors
    copyto!(ordered_values,values)
    copyto!(ordered_vectors,vectors)
    @inbounds for mode in 1:n
        source = assignment[mode]
        values[mode] = ordered_values[source]
        copyto!(@view(vectors[:, mode]), @view(ordered_vectors[:, source]))
    end
    return assignment
end
function check_eigenpairs!(
        matrix::AbstractMatrix{T},
        values::AbstractVector{T},
        vectors::AbstractMatrix{T},
        tolerance::Real,residual::AbstractVector{T}
) where {T <: Complex}
    all(isfinite, values) && all(isfinite, vectors) || return false
    R = typeof(real(zero(T)))
    condition = cond(vectors)
    isfinite(condition) && condition <= inv(sqrt(eps(R))) || return false
    scale = max(norm(matrix, Inf), eps(R))
    limit = max(convert(R, tolerance), sqrt(eps(R))) * scale
    @inbounds for mode in eachindex(values)
        mul!(residual, matrix, @view(vectors[:, mode]))
        residual .-= values[mode] .* @view(vectors[:, mode])
        norm(residual, Inf) <= limit || return false
    end
    return true
end
function recompute_matched_eigenpairs!(
        matrix::AbstractMatrix{T},
        previous_values::AbstractVector{T},
        previous_vectors::AbstractMatrix{T},work
) where {T <: Complex}
    values, vectors = _seed(matrix)
    _match!(values, vectors, previous_values, previous_vectors,work)
    @inbounds for mode in eachindex(values)
        _align!(@view(vectors[:, mode]), @view(previous_vectors[:, mode]))
    end
    return values, vectors
end
