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

- Complex charge-based potential coefficient ``p`` \\[m/F\\].

Here `r_in` and `r_ex` are the dielectric radii in meters, `κ` is complex
admittivity in S/m, and `s` is complex frequency in s⁻¹. For sinusoidal
evaluation, ``s=jω``. This definition gives ``y=s/p``. The FEM backend extracts
the inverse-admittance coefficient ``P=1/y`` in m/S and passes ``p=sP`` to
[`reduce_line_matrices!`](@ref LineCableModels.Commons.reduce_line_matrices!).
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

"""
$(TYPEDSIGNATURES)

Evaluate the selected dielectric law at frequency in Hz and temperature in °C.
`rho` is the resistivity at `temperature` \\[Ω·m\\]. By default, the selected
`temperature_dependence` law evaluates it once before the dielectric equation. That law
defaults to the linear resistivity law, and `nothing` retains the reference resistivity. A
computation that evaluates the law once for all frequencies passes `rho` instead.

The temporary material records the operating temperature when resistivity changes. The source
material remains unchanged. The lossless dielectric law ignores resistivity and loss
tangent. Returns complex admittivity in S/m.
"""
function constitutive(
        formula::Union{InsulationAdmittanceFormulation, SemiconAdmittanceFormulation},
        material::Material, frequency::Real, temperature::Real;
        temperature_dependence = TemperatureDependent.Formula(:default), workspace = nothing,
        rho = constitutive(temperature_dependence, material, temperature; workspace))
    if rho != material.rho
        material = Material(material.kind, rho, material.eps_r, material.mu_r,
            temperature, material.alpha; rho_thermal = material.rho_thermal,
            theta_max = material.theta_max, tan_delta = material.tan_delta,
            sigma_solar = material.sigma_solar)
    end
    return formula(material, frequency, temperature; workspace)
end

"""
$(TYPEDSIGNATURES)

Evaluate the dielectric law that the material's class selects. `laws` lists the selected
insulation and semicon laws, in that order. An `:insulator` takes the insulation law, and a
`:semicon` takes the semicon law. The keywords pass to the selected law. Returns complex
admittivity in S/m.
"""
function constitutive(laws::Tuple{I, S}, material::Material{<:Real, :insulator},
        frequency::Real, temperature::Real; kwargs...
) where {I <: InsulationAdmittanceFormulation, S <: SemiconAdmittanceFormulation}
    return constitutive(first(laws), material, frequency, temperature; kwargs...)
end

function constitutive(laws::Tuple{I, S}, material::Material{<:Real, :semicon},
        frequency::Real, temperature::Real; kwargs...
) where {I <: InsulationAdmittanceFormulation, S <: SemiconAdmittanceFormulation}
    return constitutive(last(laws), material, frequency, temperature; kwargs...)
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

Evaluate the registered insulation and semicon constitutive relations for every physical
dielectric layer. The result excludes the radial geometry. `rho` holds the resistivity of
each layer at `temperature` \\[Ω·m\\], which the computation evaluates once with its
temperature law.

# Returns

- `admittivity`, overwritten in physical radial-layer order \\[S/m\\].
"""
function dielectric!(
        admittivity::AbstractVector{Complex{T}},
        cable::LocalCableData{T},
        rho::AbstractVector{T},
        formulas::NamedTuple,
        frequency::T,
        temperature::T; workspace = nothing
) where {T <: Real}
    laws = (formulas.insulation_admittance, formulas.semicon_admittance)
    @inbounds for (layer, material) in pairs(cable.dielectric_materials)
        admittivity[layer] = constitutive(laws, material, frequency, temperature;
            rho = rho[layer], workspace)
    end
    return admittivity
end

"""
$(TYPEDSIGNATURES)

Assemble the unreduced cable-internal potential-coefficient matrix from the
series sum of physical radial dielectric layers. The completed `shunt` blocks replace
the radial intervals that they cover.

# Returns

- `destination`, overwritten with the local potential coefficients \\[m/F\\].
"""
function cable_potential!(
        destination::AbstractMatrix{Complex{T}},
        cable::LocalCableData{T},
        shunt::NamedTuple,
        admittivity::AbstractVector{Complex{T}},
        s::Complex{T},
        layer_coefficients::AbstractVector{Complex{T}},
        coefficients::AbstractVector{Complex{T}},
        tails::AbstractVector{Complex{T}}
) where {T <: Real}
    fill!(destination, zero(Complex{T}))
    @. layer_coefficients = potential_coefficient(
        cable.r_layer_in, cable.r_layer_ext, admittivity, s)
    @inbounds for conductors in cable.assemblies
        count = length(conductors)
        for component in 1:count
            coefficients[component] = shunt.covered[conductors[component]] ? zero(s) :
                                      radial_coefficient(
                layer_coefficients,
                cable.dielectric_ranges[conductors[component]]
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
    return _shunt_potential!(destination, shunt.blocks)
end

"""
$(TYPEDSIGNATURES)

Assemble the earth-free nodal shunt-admittance matrix of independent
concentric assemblies.

For each conductor-owned radial interval, the physical layer coefficients
combine in series as ``p=\\sum_k p_k`` and contribute the branch admittance
``y=s/p``. An interval outside the last retained conductor terminates at the
external reference. The completed `shunt` blocks replace the intervals that they
cover.

# Returns

- `destination`, overwritten with the local shunt admittance [S/m].
"""
function cable_admittance!(
        destination::AbstractMatrix{Complex{T}},
        cable::LocalCableData{T},
        shunt::NamedTuple,
        admittivity::AbstractVector{Complex{T}},
        s::Complex{T},
        layer_coefficients::AbstractVector{Complex{T}}
) where {T <: Real}
    fill!(destination, zero(Complex{T}))
    @. layer_coefficients = potential_coefficient(
        cable.r_layer_in, cable.r_layer_ext, admittivity, s)
    @inbounds for conductors in cable.assemblies
        count = length(conductors)
        for position in 1:count
            index = conductors[position]
            shunt.covered[index] && continue
            coefficient = radial_coefficient(
                layer_coefficients,
                cable.dielectric_ranges[index]
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
    return _shunt_admittance!(destination, shunt.blocks, s)
end

"""
$(TYPEDSIGNATURES)

Add the already calculated exterior potential coefficients to the cable-local
primitive potential matrix \\[m/F\\] and retain its optional trace at `frequency`.
The caller then forms admittance from this matrix. No equation or
material law is evaluated here. Return the mutated `destination`.
"""
function admittance!(
        destination::AbstractMatrix{Complex{T}},
        workspace::LineParametersWorkspace{T},
        frequency::Int
) where {T <: Real}
    input = workspace.input
    indices = workspace.plan.cable_indices
    earth_matrix = workspace.buffers.Pearth
    trace = workspace.trace

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
                destination[column, row] += earth_matrix[right, left]
            end
        end
    end
    _stash!(trace, :P, frequency, destination)
    return destination
end
