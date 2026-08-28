function compute_impedance_matrix!(
        destination::AbstractMatrix{Complex{T}},
        workspace::LineParametersWorkspace{T},
        frequency::Int,
        formulation::LineParametersFormulation
) where {T <: Real}
    input = workspace.normalized
    rho_cond = workspace.prepared.rho_cond
    earth = workspace.prepared.earth
    indices = workspace.prepared.cable_indices
    cables = workspace.prepared.cable_representatives
    earth_matrix = workspace.buffers.earth_matrix
    capture = workspace.capture
    fill!(destination, zero(Complex{T}))
    compute_earth_return_matrix!(
        earth_matrix, cables, input, earth, frequency,
        formulation.methods.earth_impedance
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
            outside = formulation.methods.internal_impedance(:outer, inner, outer, rho, mu, s)
            inside = position < count ?
                     formulation.methods.internal_impedance(
                :inner,
                input.r_in[conductors[position + 1]],
                input.r_ext[conductors[position + 1]],
                rho_cond[conductors[position + 1]],
                input.mu_cond[conductors[position + 1]],
                s
            ) : zero(outside)
            mutual = formulation.methods.internal_impedance(:mutual, inner, outer, rho, mu, s)
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
