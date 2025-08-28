import ..LineCableModels: get_description

"""
$(TYPEDEF)

Abstract type representing different frequency-dependent earth models (FDEM). Used in the multi-dispatch implementation of [`_calc_earth_properties`](@ref).

# Currently available formulations

- [`CPEarth`](@ref): Constant properties (CP) model.
"""
abstract type AbstractFDEMFormulation end

"""
$(TYPEDEF)

Represents an earth model with constant properties (CP), i.e. frequency-invariant electromagnetic properties.
"""
struct CPEarth <: AbstractFDEMFormulation end
get_description(::CPEarth) = "Constant properties (CP)"


"""
$(TYPEDSIGNATURES)

Functor implementation for `CPEarth`.

Computes frequency-dependent earth properties using the [`CPEarth`](@ref) formulation, which assumes frequency-invariant values for resistivity, permittivity, and permeability.

# Arguments

- `frequencies`: Vector of frequency values \\[Hz\\].
- `base_rho_g`: Base (DC) electrical resistivity of the soil \\[Ω·m\\].
- `base_epsr_g`: Base (DC) relative permittivity of the soil \\[dimensionless\\].
- `base_mur_g`: Base (DC) relative permeability of the soil \\[dimensionless\\].
- `formulation`: Instance of a subtype of [`AbstractFDEMFormulation`](@ref) defining the computation method.

# Returns

- `rho`: Vector of resistivity values \\[Ω·m\\] at the given frequencies.
- `epsilon`: Vector of permittivity values \\[F/m\\] at the given frequencies.
- `mu`: Vector of permeability values \\[H/m\\] at the given frequencies.

# Examples

```julia
frequencies = [1e3, 1e4, 1e5]

# Using the CP model
rho, epsilon, mu = $(FUNCTIONNAME)(frequencies, 100, 10, 1, CPEarth())
println(rho)     # Output: [100, 100, 100]
println(epsilon) # Output: [8.854e-11, 8.854e-11, 8.854e-11]
println(mu)      # Output: [1.2566e-6, 1.2566e-6, 1.2566e-6]
```

# See also

- [`EarthLayer`](@ref)
"""
function (f::CPEarth)(frequencies::Vector{T}, base_rho_g::T, base_epsr_g::T,
    base_mur_g::T) where {T<:REALSCALAR}

    # Preallocate for performance
    n_freq = length(frequencies)
    rho = Vector{T}(undef, n_freq)
    epsilon = Vector{typeof(ε₀ * base_epsr_g)}(undef, n_freq)
    mu = Vector{typeof(μ₀ * base_mur_g)}(undef, n_freq)

    # Vectorized assignment
    fill!(rho, base_rho_g)
    fill!(epsilon, ε₀ * base_epsr_g)
    fill!(mu, μ₀ * base_mur_g)

    return rho, epsilon, mu
end

function (f::CPEarth)(frequencies::AbstractVector, base_rho_g, base_epsr_g, base_mur_g)
    T = resolve_T(frequencies, base_rho_g, base_epsr_g, base_mur_g)
    return f(
        coerce_to_T(frequencies, T),
        coerce_to_T(base_rho_g, T),
        coerce_to_T(base_epsr_g, T),
        coerce_to_T(base_mur_g, T),
    )
end
