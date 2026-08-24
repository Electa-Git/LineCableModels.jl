function compute_admittance_matrix!(
        destination::AbstractMatrix{Complex{T}},
        input::AnalyticalInput{T},
        earth,
        frequency::Int,
        formulation::AnalyticalFormulation,
        trace
) where {T <: Real}
    fill!(destination, zero(Complex{T}))
    indices, cables = _cable_indices(input)
    earth_matrix = Matrix{Complex{T}}(undef, input.n_cables, input.n_cables)
    compute_earth_return_matrix!(
        earth_matrix, cables, input, earth, frequency,
        formulation.earth_admittance
    )
    _stash!(_trace_target(trace, :Pg), frequency, earth_matrix)
    s = input.jω[frequency]

    @inbounds for cable in 1:input.n_cables
        conductors = indices[cable]
        count = length(conductors)
        coefficients = Vector{Complex{T}}(undef, count)
        for component in 1:count
            coefficients[component] = InsulationAdmittance.potential_coefficient(
                formulation.insulation_admittance, input, conductors[component], s
            )
        end
        tails = Vector{Complex{T}}(undef, count)
        tails[count] = coefficients[count]
        for gap in (count - 1):-1:1
            tails[gap] = coefficients[gap] + tails[gap + 1]
        end
        for row in 1:count, column in 1:count

            destination[conductors[row], conductors[column]] += tails[max(row, column)]
        end
    end
    _stash!(_trace_target(trace, :Pin), frequency, destination)

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
    _stash!(_trace_target(trace, :P), frequency, destination)
    return destination
end
