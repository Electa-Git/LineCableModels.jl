function compute_admittance_matrix!(
        destination::AbstractMatrix{Complex{T}},
        workspace::LineParametersWorkspace{T},
        frequency::Int,
        formulation::LineParametersFormulation
) where {T <: Real}
    input = workspace.normalized
    earth = workspace.prepared.earth
    indices = workspace.prepared.cable_indices
    cables = workspace.prepared.cable_representatives
    earth_matrix = workspace.buffers.earth_matrix
    coefficients = workspace.buffers.coefficients
    tails = workspace.buffers.tails
    capture = workspace.capture
    fill!(destination, zero(Complex{T}))
    compute_earth_return_matrix!(
        earth_matrix, cables, input, earth, frequency,
        formulation.methods.earth_admittance
    )
    _stash!(_capture_target(capture, :Pg), frequency, earth_matrix)
    s = input.jω[frequency]

    @inbounds for cable in 1:input.n_cables
        conductors = indices[cable]
        count = length(conductors)
        for component in 1:count
            coefficients[component] = InsulationAdmittance.potential_coefficient(
                formulation.methods.insulation_admittance, input, conductors[component], s
            )
        end
        tails[count] = coefficients[count]
        for gap in (count - 1):-1:1
            tails[gap] = coefficients[gap] + tails[gap + 1]
        end
        for row in 1:count, column in 1:count

            destination[conductors[row], conductors[column]] += tails[max(row, column)]
        end
    end
    _stash!(_capture_target(capture, :Pin), frequency, destination)

    @inbounds for cable in 1:input.n_cables
        self = earth_matrix[cable, cable]
        for row in indices[cable], column in indices[cable]

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
    _stash!(_capture_target(capture, :P), frequency, destination)
    return destination
end
