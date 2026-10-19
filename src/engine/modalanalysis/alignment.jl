"""
$(TYPEDSIGNATURES)

Align the modal `parameters` with `reference` independently at each frequency.
Both results must have the same frequency and phase axes and coefficient basis.

Cluster reference modes by connected components of squared propagation roots.
Reference roots share an edge when their relative separation is below `sqrt(eps(T))`,
where `T` is the real scalar type of the assignment costs. Equal zero roots share a cluster.

Maximize total normalized current-column projection onto the reference cluster
spans with an optimal assignment. For a singleton this projection is the ordinary
absolute overlap. Apply the unitary polar factor of each assigned-to-reference
cross-overlap matrix to both paired bases and by similarity to the modal coefficients.
Keep the input roots in their assigned order.

# Returns

- `parameters`: Aligned line parameters with paired voltage and current bases.
- `permutation`: Reference-to-input mode indices, with axes mode × frequency.
- `overlap`: Assigned normalized projection norms \\[dimensionless\\], with axes mode × frequency.
- `cluster_sizes`: Reference cluster sizes for each mode, with axes mode × frequency.

# Notes

The caller selects any acceptance threshold. The optional `buffers` come from
`initialize_buffers(align_modes, T, parameters, (;), (;))`, with real scalar type
`T`. Returned arrays share those buffers, and the next call overwrites them.
Omitting `buffers` allocates independent arrays.
"""
function align_modes(parameters::LineParameters{C,R,D,Basis},
        reference::LineParameters{<:Complex,<:Real,<:ModalDomain};
        buffers=initialize_buffers(align_modes,
            typeof(real(zero(eltype(Ti(parameters))))), parameters, (;), (;))) where {C,R,D<:ModalDomain,Basis}
    validate(parameters, align_modes, reference)
    current, nominal_current = Ti(parameters), Ti(reference)
    voltage = Tv(parameters)
    roots, nominal_roots = gamma(parameters), gamma(reference)
    n = size(current,2)
    tolerance = sqrt(eps(eltype(buffers.cost)))
    for frequency in axes(current,3)
        for mode in 1:n
            buffers.cluster_ids[mode] = mode
        end
        for left in 1:n, right in left+1:n
            a, b = nominal_roots[left,frequency]^2, nominal_roots[right,frequency]^2
            scale = max(abs(a),abs(b))
            separation = iszero(scale) ? zero(scale) : abs(a-b)/scale
            if separation < tolerance
                kept = min(buffers.cluster_ids[left],buffers.cluster_ids[right])
                merged = max(buffers.cluster_ids[left],buffers.cluster_ids[right])
                for mode in 1:n
                    buffers.cluster_ids[mode] == merged && (buffers.cluster_ids[mode] = kept)
                end
            end
        end
        for cluster in 1:n
            buffers.cluster_ids[cluster] == cluster || continue
            count = 0
            for mode in 1:n
                if buffers.cluster_ids[mode] == cluster
                    count += 1
                    buffers.members[count] = mode
                    for phase in 1:n
                        buffers.columns[phase,count] = nominal_current[phase,mode,frequency]
                    end
                end
            end
            orthogonal = qr!(@view(buffers.columns[:,1:count]))
            for mode in 1:n
                right = @view current[:,mode,frequency]
                copyto!(buffers.projection,right)
                lmul!(adjoint(orthogonal.Q),buffers.projection)
                overlap = clamp(norm(@view(buffers.projection[1:count]))/norm(right),0,1)
                for member in 1:count
                    buffers.cost[buffers.members[member],mode] = 1-overlap
                end
            end
            for member in 1:count
                buffers.cluster_sizes[buffers.members[member],frequency] = count
            end
        end
        permutation = hungarian_assignment!(buffers.cost,buffers)
        for mode in 1:n
            buffers.permutation[mode,frequency] = permutation[mode]
            buffers.overlap[mode,frequency] = 1-buffers.cost[mode,permutation[mode]]
            buffers.gamma[mode,frequency] = roots[permutation[mode],frequency]
        end
        fill!(buffers.rotation,zero(eltype(buffers.rotation)))
        for cluster in 1:n
            buffers.cluster_ids[cluster] == cluster || continue
            count = 0
            for mode in 1:n
                if buffers.cluster_ids[mode] == cluster
                    count += 1
                    buffers.members[count] = mode
                end
            end
            for row in 1:count, column in 1:count
                buffers.cross[row,column] = dot(
                    @view(current[:,permutation[buffers.members[row]],frequency]),
                    @view(nominal_current[:,buffers.members[column],frequency]))
            end
            polar = svd!(@view(buffers.cross[1:count,1:count]))
            for row in 1:count, column in 1:count
                value = zero(eltype(buffers.rotation))
                for mode in 1:count
                    value += polar.U[row,mode]*polar.Vt[mode,column]
                end
                buffers.rotation[buffers.members[row],buffers.members[column]] = value
            end
        end
        for (source,destination) in ((current,buffers.Ti),(voltage,buffers.Tv))
            for mode in 1:n, phase in 1:n
                buffers.columns[phase,mode] = source[phase,permutation[mode],frequency]
            end
            mul!(@view(destination[:,:,frequency]),buffers.columns,buffers.rotation)
        end
        for (source,destination) in ((parameters.Z,buffers.Z),(parameters.Y,buffers.Y))
            for column in 1:n, row in 1:n
                buffers.coefficients[row,column] = source[permutation[row],permutation[column],frequency]
            end
            mul!(buffers.product,buffers.coefficients,buffers.rotation)
            mul!(@view(destination[:,:,frequency]),adjoint(buffers.rotation),buffers.product)
        end
    end
    domain = ModalDomain(ModalOperators(buffers.Tv,buffers.Ti),buffers.gamma)
    aligned = LineParameters(domain,
        SeriesImpedance{eltype(buffers.Z),Basis}(buffers.Z),
        ShuntAdmittance{eltype(buffers.Y),Basis}(buffers.Y),
        frequencies(parameters),details(parameters))
    return (parameters=aligned, permutation=buffers.permutation, overlap=buffers.overlap,
        cluster_sizes=buffers.cluster_sizes)
end

"""
$(TYPEDSIGNATURES)

Validate equal modal axes and coefficient bases and finite nonzero current-column
norm products before alignment. Return `parameters`.
"""
function validate(parameters::LineParameters{<:Complex,<:Real,<:ModalDomain},
        ::typeof(align_modes), reference::LineParameters{<:Complex,<:Real,<:ModalDomain})
    frequencies(parameters) == frequencies(reference) || throw(DimensionMismatch(
        "modal alignment requires equal frequency axes"))
    basis(parameters) === basis(reference) || throw(ArgumentError(
        "modal alignment requires equal coefficient bases"))
    current, nominal_current = Ti(parameters), Ti(reference)
    size(current) == size(nominal_current) || throw(DimensionMismatch(
        "modal alignment requires equal phase and mode axes"))
    for frequency in axes(current,3), nominal_mode in axes(current,2), mode in axes(current,2)
        denominator = norm(@view(nominal_current[:,nominal_mode,frequency])) *
            norm(@view(current[:,mode,frequency]))
        isfinite(denominator) && denominator > 0 || throw(ArgumentError(
            "modal alignment requires finite nonzero current columns"))
    end
    return parameters
end

"""
$(TYPEDSIGNATURES)

Extend `buffers` with arrays for optimal modal assignment and the aligned
coefficient matrices, propagation roots and paired bases of `input`.
"""
function initialize_buffers(::typeof(align_modes), ::Type{T},
        input::LineParameters{<:Complex,<:Real,<:ModalDomain}, plan, buffers) where {T<:Real}
    n, _, nf = size(input.Z)
    return merge(buffers, (cost=Matrix{T}(undef,n,n), u=zeros(T,n+1), v=zeros(T,n+1),
        matching=zeros(Int,n+1), way=zeros(Int,n+1), minimums=Vector{T}(undef,n+1),
        used=falses(n+1), assignment=Vector{Int}(undef,n),
        permutation=Matrix{Int}(undef,n,nf), overlap=Matrix{T}(undef,n,nf),
        cluster_ids=Vector{Int}(undef,n), members=Vector{Int}(undef,n),
        cluster_sizes=Matrix{Int}(undef,n,nf), projection=Vector{Complex{T}}(undef,n),
        columns=Matrix{Complex{T}}(undef,n,n), cross=Matrix{Complex{T}}(undef,n,n),
        rotation=Matrix{Complex{T}}(undef,n,n), coefficients=Matrix{Complex{T}}(undef,n,n),
        product=Matrix{Complex{T}}(undef,n,n), Z=similar(input.Z.values), Y=similar(input.Y.values),
        gamma=similar(gamma(input)), Tv=similar(Tv(input)), Ti=similar(Ti(input))))
end
