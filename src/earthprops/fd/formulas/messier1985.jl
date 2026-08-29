"Return the high-frequency permittivity assumption of the Messier relation."
assumptions(::Val{:Messier1985}) = (epsilon_infinity = 8.0,)

description(::Formula{:Messier1985}) = "Messier 1985 (square-root soil dispersion)"

#=
Evaluate the Messier square-root soil-dispersion relation.

Input and returned resistivity use Ω·m; conductivity is evaluated internally in
S/m and relative permittivity is dimensionless.

# Reference

M. Messier, *Another Soil Conductivity Model*, internal report, JAYCOR,
Santa Barbara, 1985.
=#
function (::Functor{:Messier1985})(
        material::EarthMaterial{T},
        frequency::T,
        values::NamedTuple
) where {T <: Real}
    conductivity_reference = inv(material.rho)
    epsilon_infinity = convert(T, values.epsilon_infinity)
    epsilon0 = vacuum_permittivity(frequency)
    pi_typed = one(frequency) * π
    relative_permittivity = epsilon_infinity + sqrt(
        conductivity_reference * epsilon_infinity /
        (pi_typed * frequency * epsilon0)
    )
    conductivity = conductivity_reference + sqrt(
        4 * pi_typed * frequency * conductivity_reference *
        epsilon0 * epsilon_infinity
    )
    return EarthMaterial{T}(inv(conductivity), relative_permittivity, material.mu_r)
end

# Return the stable discovery identifier.
:Messier1985
