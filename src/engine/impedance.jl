"""
$(TYPEDSIGNATURES)

Assemble the earth-free primitive series-impedance matrix of independent
concentric cable assemblies.

The selected internal- and insulation-impedance formulas contribute their
outer, inner, transfer, and longitudinal-insulation terms directly to
`destination`. The matrix remains unreduced.

# Arguments

- `destination`: Reusable primitive series-impedance matrix [Ω/m].
- `input`: Concrete local cable arrays.
- `rho_cond`: Temperature-corrected conductor resistivities [Ω·m].
- `methods`: Resolved formulation methods.
- `s`: Complex angular frequency ``jω`` [1/s].

# Returns

- `destination`, overwritten with the assembled local matrix.
"""
function cable_impedance!(
        destination::AbstractMatrix{Complex{T}},
        input::LocalCableData{T},
        rho_cond::AbstractVector{T},
        methods::NamedTuple,
        s::Complex{T}
) where {T <: Real}
    fill!(destination, zero(Complex{T}))
    @inbounds for conductors in input.assemblies
        count = length(conductors)
        inside = zero(eltype(destination))
        for position in count:-1:1
            index = conductors[position]
            surfaces = InternalImpedance.surface_impedances(methods.internal_impedance,
                position > 1 ? Val((:outer,:transfer,:inner)) : Val((:outer,)),
                input.r_in[index],
                input.r_ext[index],
                rho_cond[index],
                input.mu_cond[index],
                s
            )
            outside = surfaces.outer
            transfer = position > 1 ? surfaces.transfer : zero(outside)
            insulation = methods.insulation_impedance(
                input.r_ext[index],
                input.r_ins_ext[index],
                input.mu_ins[index],
                s
            )
            loop = outside + inside + insulation
            if position > 1
                for row in 1:(position - 1), column in 1:(position - 1)

                    destination[conductors[row], conductors[column]] += loop - 2 * transfer
                end
                for row in 1:(position - 1)
                    destination[index, conductors[row]] += loop - transfer
                    destination[conductors[row], index] += loop - transfer
                end
            end
            destination[index, index] += loop
            # Reuse this wall's prepared state when its inner surface is needed
            # by the next contained conductor.
            position > 1 && (inside = surfaces.inner)
        end
    end
    return destination
end

function impedance!(
        destination::AbstractMatrix{Complex{T}},
        workspace::LineParametersWorkspace{T},
        frequency::Int,
        formulation::LineParametersFormulation
) where {T <: Real}
    input = workspace.input
    rho_cond = workspace.invariants.rho_cond
    indices = workspace.invariants.cable_indices
    earth_matrix = workspace.buffers.earth_matrix
    earth_media = workspace.buffers.earth_materials.earth_impedance
    capture = workspace.capture
    s = input.jω[frequency]
    cable_impedance!(
        destination,
        input.cable,
        rho_cond,
        formulation.methods,
        s
    )
    _stash!(_capture_target(capture, :Zin), frequency, destination)

    earth!(
        earth_matrix, workspace.invariants.earth_bindings.earth_impedance, earth_media, s,
        formulation.methods.earth_impedance,
        _gamma(input.Γ, frequency),
        workspace.buffers.earth_numerical.earth_impedance,
        get(earth_media, :thickness, nothing)
    )
    _stash!(_capture_target(capture, :Zg), frequency, earth_matrix)

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
                destination[column, row] += earth_matrix[right, left]
            end
        end
    end
    _stash!(_capture_target(capture, :Z), frequency, destination)
    return destination
end
