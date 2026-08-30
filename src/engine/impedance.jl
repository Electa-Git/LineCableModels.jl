function impedance!(
        destination::AbstractMatrix{Complex{T}},
        workspace::LineParametersWorkspace{T},
        frequency::Int,
        formulation::LineParametersFormulation
) where {T <: Real}
    input = workspace.input
    rho_cond = workspace.invariants.rho_cond
    indices = workspace.invariants.cable_indices
    stratified = media(formulation.methods.earth_impedance) === Val(:stratified)
    pairs = stratified ?
            workspace.invariants.earth_pairs : workspace.invariants.homogeneous_pairs
    earth_matrix = workspace.buffers.earth_matrix
    earth_media = stratified ? workspace.buffers.earth_layers :
                  workspace.buffers.earth_media
    capture = workspace.capture
    fill!(destination, zero(Complex{T}))
    earth!(
        earth_matrix, pairs, input, earth_media, frequency,
        formulation.methods.earth_impedance,
        _gamma(input.Γ, frequency),
        workspace.buffers.earth_impedance_segments,
        stratified ? earth_media.thickness : nothing
    )
    _stash!(_capture_target(capture, :Zg), frequency, earth_matrix)
    s = input.jω[frequency]

    @inbounds for cable in 1:input.n_cables
        conductors = indices[cable]
        count = length(conductors)
        for position in count:-1:1
            index = conductors[position]
            inner = input.r_in[index]
            outer = input.r_ext[index]
            rho = rho_cond[index]
            mu = input.mu_cond[index]
            interaction = formulation.methods.internal_impedance(
                inner, outer, rho, mu, s
            )
            outside = interaction(Val(:outer))
            inside = position < count ?
                     formulation.methods.internal_impedance(
                input.r_in[conductors[position + 1]],
                input.r_ext[conductors[position + 1]],
                rho_cond[conductors[position + 1]],
                input.mu_cond[conductors[position + 1]],
                s
            )(Val(:inner)) : zero(outside)
            mutual = interaction(Val(:mutual))
            insulation = formulation.methods.insulation_impedance(
                outer, input.r_ins_ext[index], input.mu_ins[index], s
            )
            loop = outside + inside + insulation
            if position > 1
                for row in 1:(position - 1), column in 1:(position - 1)

                    destination[conductors[row], conductors[column]] += loop - 2 * mutual
                end
                for row in 1:(position - 1)
                    destination[conductors[position], conductors[row]] += loop - mutual
                    destination[conductors[row], conductors[position]] += loop - mutual
                end
            end
            destination[index, index] += loop
        end
    end
    _stash!(_capture_target(capture, :Zin), frequency, destination)

    @inbounds for cable in 1:input.n_cables
        conductors = indices[cable]
        self = earth_matrix[cable, cable]
        for row in conductors, column in conductors

            destination[row, column] += self
        end
    end
    @inbounds for left in 1:(input.n_cables - 1)
        for right in (left + 1):input.n_cables
            mutual = earth_matrix[left, right]
            for row in indices[left], column in indices[right]

                destination[row, column] += mutual
                destination[column, row] += mutual
            end
        end
    end
    _stash!(_capture_target(capture, :Z), frequency, destination)
    return destination
end
