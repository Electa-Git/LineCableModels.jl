"""
$(TYPEDSIGNATURES)

Align modal `parameters` with `reference` at the frequency index returned by
`ordering_frequency(selected, frequencies(reference))`. Both results must have
the same frequency and phase axes and coefficient basis.

Maximize the sum of normalized absolute current-column overlaps using an optimal
assignment. Apply that permutation and one complex phase factor per mode to both
voltage and current bases, with the corresponding similarity transformation of
the modal coefficient matrices. Hold the assignment and phase factors fixed over
frequency. Each overlap is dimensionless and lies in `[0, 1]`.

Return a named tuple containing the aligned `parameters`, the `permutation` from
reference modes to input modes, their `overlap` and the ordering `frequency` index.
The caller selects any acceptance threshold. A zero overlap uses a unit phase factor.

The optional `buffers` come from
`initialize_buffers(align_modes, T, parameters, (;), (;))`, with real scalar type
`T`. Returned arrays share those buffers and are overwritten by the next call.
Omitting `buffers` allocates independent arrays.
"""
function align_modes(selected::AbstractFormulation,
        parameters::LineParameters{<:Complex,<:Real,<:ModalDomain},
        reference::LineParameters{<:Complex,<:Real,<:ModalDomain};
        buffers=initialize_buffers(align_modes,
            typeof(real(zero(eltype(Ti(parameters))))), parameters, (;), (;)))
    frequencies(parameters) == frequencies(reference) || throw(DimensionMismatch(
        "modal alignment requires equal frequency axes"))
    basis(parameters) === basis(reference) || throw(ArgumentError(
        "modal alignment requires equal coefficient bases"))
    current, nominal_current = Ti(parameters), Ti(reference)
    size(current) == size(nominal_current) || throw(DimensionMismatch(
        "modal alignment requires equal phase and mode axes"))
    frequency = ordering_frequency(selected, frequencies(reference))
    for nominal_mode in axes(current,2), mode in axes(current,2)
        left = @view nominal_current[:,nominal_mode,frequency]
        right = @view current[:,mode,frequency]
        denominator = norm(left)*norm(right)
        isfinite(denominator) && denominator > 0 || throw(ArgumentError(
            "modal alignment requires finite nonzero current columns"))
        buffers.cost[nominal_mode,mode] = 1-clamp(abs(dot(left,right))/denominator,0,1)
    end
    permutation = hungarian_assignment!(buffers.cost,buffers)
    for mode in eachindex(permutation)
        buffers.overlap[mode] = 1-buffers.cost[mode,permutation[mode]]
        overlap = dot(@view(nominal_current[:,mode,frequency]),
            @view(current[:,permutation[mode],frequency]))
        buffers.phases[mode] = iszero(overlap) ? one(overlap) : conj(overlap)/abs(overlap)
    end
    voltage = Tv(parameters)
    roots = gamma(parameters)
    for k in axes(current,3), mode in eachindex(permutation)
        original = permutation[mode]
        factor = buffers.phases[mode]
        buffers.gamma[mode,k] = roots[original,k]
        for phase in axes(current,1)
            buffers.Ti[phase,mode,k] = current[phase,original,k]*factor
            buffers.Tv[phase,mode,k] = voltage[phase,original,k]*factor
        end
        for row in eachindex(permutation)
            ratio = factor/buffers.phases[row]
            buffers.Z[row,mode,k] = parameters.Z[permutation[row],original,k]*ratio
            buffers.Y[row,mode,k] = parameters.Y[permutation[row],original,k]*ratio
        end
    end
    domain = ModalDomain(ModalOperators(buffers.Tv,buffers.Ti),buffers.gamma)
    aligned = LineParameters(domain,buffers.Z,buffers.Y,frequencies(parameters);
        basis=basis(parameters),details=details(parameters))
    return (parameters=aligned, permutation, overlap=buffers.overlap, frequency)
end

"""
$(TYPEDSIGNATURES)

Extend `buffers` with arrays for optimal modal assignment and the aligned
coefficient matrices, propagation roots and paired bases of `input`.
"""
function initialize_buffers(::typeof(align_modes), ::Type{T},
        input::LineParameters{<:Complex,<:Real,<:ModalDomain}, plan, buffers) where {T<:Real}
    n = size(input.Z,1)
    return merge(buffers, (cost=Matrix{T}(undef,n,n), u=zeros(T,n+1), v=zeros(T,n+1),
        matching=zeros(Int,n+1), way=zeros(Int,n+1), minimums=Vector{T}(undef,n+1),
        used=falses(n+1), assignment=Vector{Int}(undef,n), overlap=Vector{T}(undef,n),
        phases=Vector{Complex{T}}(undef,n), Z=similar(input.Z.values), Y=similar(input.Y.values),
        gamma=similar(gamma(input)), Tv=similar(Tv(input)), Ti=similar(Ti(input))))
end
