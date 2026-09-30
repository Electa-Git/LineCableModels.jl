# Shared numerical operations for frequency-tracked eigensystems.

# Minimize the current column's imaginary norm, preserve the voltage/current
# pairing, then resolve the remaining sign against the previous frequency.
function orient_modes!(voltage::AbstractArray{T, 3}, current::AbstractArray{T, 3},
        rotate::Bool) where {T <: Complex}
    R=typeof(real(zero(T)))
    for frequency in axes(current, 3), mode in axes(current, 2)

        vector=@view current[:, mode, frequency]
        phase=one(T)
        scale=maximum(value -> max(abs(real(value)), abs(imag(value))), vector)
        if rotate && !iszero(scale)
            real_sum=zero(R)
            imaginary_sum=zero(R)
            for value in vector
                scaled=value/scale
                real_sum += real(scaled)^2-imag(scaled)^2
                imaginary_sum += 2real(scaled)*imag(scaled)
            end
            phase=cis(-atan(imaginary_sum, real_sum)/2)
        end
        if frequency>firstindex(current, 3) && !iszero(scale)
            previous=@view current[:, mode, frequency - 1]
            previous_scale=maximum(value -> max(abs(real(value)), abs(imag(value))), previous)
            if !iszero(previous_scale)
                overlap=zero(T)
                for row in eachindex(vector)
                    overlap += conj(previous[row]/previous_scale)*(vector[row]/scale)
                end
                real(overlap*phase)<zero(R) && (phase=-phase)
            end
        end
        if phase!=one(T)
            vector .*= phase
            @views voltage[:, mode, frequency] .*= phase
        end
    end
    return nothing
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
        component_magnitude = abs(vector[index])
        if component_magnitude > magnitude
            pivot = index
            magnitude = component_magnitude
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
            for column_slot in 2:(n + 1)
                used[column_slot] && continue
                reduced = cost[matched_row, column_slot - 1] -
                          u[matched_row + 1] - v[column_slot]
                if reduced < minimums[column_slot]
                    minimums[column_slot] = reduced
                    way[column_slot] = column
                end
                if minimums[column_slot] < delta
                    delta = minimums[column_slot]
                    next_column = column_slot
                end
            end
            isfinite(delta) || throw(ArgumentError(
                "modal assignment cost must be finite"
            ))
            for column_slot in 1:(n + 1)
                if used[column_slot]
                    u[matching[column_slot] + 1] += delta
                    v[column_slot] -= delta
                elseif column_slot > 1
                    minimums[column_slot] -= delta
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

# Complex eigenpair equations with the bilinear normalization tᵀt = 1.
function eigenpair_residual!(residual, x, matrix)
    n = size(matrix, 1)
    vector = @view x[1:n]
    mul!(@view(residual[1:n]), matrix, vector)
    constraint = -one(eltype(x))
    @inbounds for row in 1:n
        residual[row] -= x[end] * vector[row]
        constraint += vector[row] * vector[row]
    end
    residual[end] = constraint
    return residual
end

function eigenpair_jacobian!(jacobian, x, matrix)
    n = size(matrix, 1)
    copyto!(@view(jacobian[1:n, 1:n]), matrix)
    @inbounds for row in 1:n
        jacobian[row, row] -= x[end]
        jacobian[row, end] = -x[row]
        jacobian[end, row] = 2x[row]
    end
    jacobian[end, end] = zero(eltype(jacobian))
    return jacobian
end

# Greedy square assignment: choose the smallest remaining entry.
function greedy_assignment!(assignment, cost)
    for _ in eachindex(assignment)
        entry = argmin(cost)
        assignment[entry[1]] = entry[2]
        @views cost[entry[1], :] .= Inf
        @views cost[:, entry[2]] .= Inf
    end
    return assignment
end
