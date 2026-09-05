"""
$(TYPEDSIGNATURES)

Calculate the shunt admittance of one homogeneous coaxial annulus:

```math
y=\\frac{2\\pi\\kappa}{\\log(r_{ex}/r_{in})}.
```

# Returns

- Complex shunt admittance per unit length ``y`` \\[S/m\\].
"""
@inline function layer_admittance(r_in::T, r_ex::T, κ::Complex{T}) where {T <: Real}
    return 2 * (one(T) * π) * κ / log(r_ex / r_in)
end

"""
$(TYPEDSIGNATURES)

Calculate one homogeneous coaxial annulus's potential coefficient
``p=s/y`` from its registered constitutive relation's complex admittivity. A
zero inner radius or zero-thickness annulus contributes zero.

# Returns

- Complex potential coefficient per unit length ``p`` \\[m/F\\].
"""
@inline function potential_coefficient(
        r_in::T,
        r_ex::T,
        κ::Complex{T},
        s::Complex{T}
) where {T <: Real}
    if isapprox(r_in, zero(T); atol = eps(T)) ||
       isapprox(r_in, r_ex; atol = eps(T))
        return zero(Complex{T})
    end
    y = layer_admittance(r_in, r_ex, κ)
    return s / y
end

@inline function radial_coefficient(coefficients, layers::UnitRange{Int})
    coefficient = zero(eltype(coefficients))
    @inbounds for layer in layers
        coefficient += coefficients[layer]
    end
    return coefficient
end

"""
$(TYPEDSIGNATURES)

Evaluate the registered insulation and semicon constitutive relations and
store the potential coefficient of every physical dielectric layer.

# Returns

- `coefficients`, overwritten in physical radial-layer order [m/F].
"""
function dielectric!(
        coefficients::AbstractVector{Complex{T}},
        input::LocalCableData{T},
        methods::NamedTuple,
        frequency::T,
        temperature::T,
        s::Complex{T}
) where {T <: Real}
    @inbounds for layer in input.insulation_indices
        κ = constitutive(
            methods.insulation_admittance,
            input.dielectric_materials[layer],
            frequency,
            temperature
        )
        coefficients[layer] = potential_coefficient(
            input.r_layer_in[layer],
            input.r_layer_ext[layer],
            κ,
            s
        )
    end
    @inbounds for layer in input.semicon_indices
        κ = constitutive(
            methods.semicon_admittance,
            input.dielectric_materials[layer],
            frequency,
            temperature
        )
        coefficients[layer] = potential_coefficient(
            input.r_layer_in[layer],
            input.r_layer_ext[layer],
            κ,
            s
        )
    end
    return coefficients
end

"""
$(TYPEDSIGNATURES)

Assemble the unreduced cable-internal potential-coefficient matrix from the
series sum of physical radial dielectric layers.

# Returns

- `destination`, overwritten with the local potential coefficients [m/F].
"""
function cable_potential!(
        destination::AbstractMatrix{Complex{T}},
        input::LocalCableData{T},
        methods::NamedTuple,
        frequency::T,
        temperature::T,
        s::Complex{T},
        layer_coefficients::AbstractVector{Complex{T}},
        coefficients::AbstractVector{Complex{T}},
        tails::AbstractVector{Complex{T}}
) where {T <: Real}
    fill!(destination, zero(Complex{T}))
    dielectric!(layer_coefficients, input, methods, frequency, temperature, s)
    @inbounds for conductors in input.assemblies
        count = length(conductors)
        for component in 1:count
            coefficients[component] = radial_coefficient(
                layer_coefficients,
                input.dielectric_ranges[conductors[component]]
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
    return destination
end

"""
$(TYPEDSIGNATURES)

Assemble the earth-free nodal shunt-admittance matrix of independent
concentric assemblies.

For each conductor-owned radial interval, the physical layer coefficients
combine in series as ``p=\\sum_k p_k`` and contribute the branch admittance
``y=s/p``. An interval outside the last retained conductor terminates at the
external reference.

# Returns

- `destination`, overwritten with the local shunt admittance [S/m].
"""
function cable_admittance!(
        destination::AbstractMatrix{Complex{T}},
        input::LocalCableData{T},
        methods::NamedTuple,
        frequency::T,
        temperature::T,
        s::Complex{T},
        layer_coefficients::AbstractVector{Complex{T}}
) where {T <: Real}
    fill!(destination, zero(Complex{T}))
    dielectric!(layer_coefficients, input, methods, frequency, temperature, s)
    @inbounds for conductors in input.assemblies
        count = length(conductors)
        for position in 1:count
            index = conductors[position]
            coefficient = radial_coefficient(
                layer_coefficients,
                input.dielectric_ranges[index]
            )
            iszero(coefficient) && continue
            admittance = s / coefficient
            destination[index, index] += admittance
            if position < count
                outside = conductors[position + 1]
                destination[outside, outside] += admittance
                destination[index, outside] -= admittance
                destination[outside, index] -= admittance
            end
        end
    end
    return destination
end

function admittance!(
        destination::AbstractMatrix{Complex{T}},
        workspace::LineParametersWorkspace{T},
        frequency::Int,
        formulation::LineParametersFormulation
) where {T <: Real}
    input = workspace.input
    indices = workspace.invariants.cable_indices
    stratified = media(formulation.methods.earth_admittance) === Val(:stratified)
    pairs = stratified ?
            workspace.invariants.earth_pairs : workspace.invariants.homogeneous_pairs
    earth_matrix = workspace.buffers.earth_matrix
    earth_media = stratified ? workspace.buffers.earth_layers :
                  workspace.buffers.earth_media
    coefficients = workspace.buffers.coefficients
    tails = workspace.buffers.tails
    layer_coefficients = workspace.buffers.layer_coefficients
    capture = workspace.capture
    s = input.jω[frequency]
    cable_potential!(
        destination,
        input.cable,
        formulation.methods,
        input.freq[frequency],
        input.temperature,
        s,
        layer_coefficients,
        coefficients,
        tails
    )
    _stash!(_capture_target(capture, :Pin), frequency, destination)

    earth!(
        earth_matrix, pairs, input, earth_media, frequency,
        formulation.methods.earth_admittance,
        _gamma(input.Γ, frequency),
        workspace.buffers.earth_admittance_segments,
        stratified ? earth_media.thickness : nothing
    )
    _stash!(_capture_target(capture, :Pg), frequency, earth_matrix)

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
