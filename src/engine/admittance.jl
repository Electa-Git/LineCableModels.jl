function admittance!(
        destination::AbstractMatrix{Complex{T}},
        workspace::LineParametersWorkspace{T},
        frequency::Int,
        formulation::LineParametersFormulation
) where {T <: Real}
    input = workspace.input
    indices = workspace.invariants.cable_indices
    pairs = workspace.invariants.homogeneous_pairs
    earth_matrix = workspace.buffers.earth_matrix
    earth_media = workspace.buffers.earth_media
    coefficients = workspace.buffers.coefficients
    tails = workspace.buffers.tails
    capture = workspace.capture
    fill!(destination, zero(Complex{T}))
    earth!(
        earth_matrix, pairs, input, earth_media, frequency,
        formulation.methods.earth_admittance,
        _gamma(input.Γ, frequency),
        workspace.buffers.earth_segments
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
