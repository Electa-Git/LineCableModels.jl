const LONGMIRE_SMITH_COEFFICIENTS = (
    3.4e6, 2.74e5, 2.58e4, 3.38e3, 5.26e2, 1.33e2, 2.72e1,
    1.25e1, 4.8, 2.17, 0.98, 0.392, 0.173
)

"Return the high-frequency and relaxation assumptions of Longmire–Smith."
assumptions(::Val{:LongmireSmith1975}) = (
    epsilon_infinity = 5.0,
    corner_scale = 125.0,
    corner_exponent = 0.8312
)

description(::Formula{:LongmireSmith1975}) =
    "Longmire–Smith 1975 (13-term dielectric relaxation)"

#=
Evaluate the 13-term Longmire–Smith soil-dispersion relation.

The tabulated relaxation coefficients are embedded as a typed tuple rather than
loaded from a runtime text file. Input and returned resistivity use Ω·m;
relative permittivity is dimensionless.

# Reference

C. L. Longmire and K. S. Smith, *A Universal Impedance for Soils*, Defense
Nuclear Agency, 1975.
=#
function (::Functor{:LongmireSmith1975})(
        material::EarthMaterial{T},
        frequency::T,
        values::NamedTuple
) where {T <: Real}
    conductivity_reference = inv(material.rho)
    corner = (
        convert(T, values.corner_scale) * conductivity_reference
    )^convert(T, values.corner_exponent)
    permittivity_sum = zero(frequency)
    conductivity_sum = zero(frequency)
    decade = one(frequency)
    ten = convert(T, 10)
    @inbounds for coefficient in LONGMIRE_SMITH_COEFFICIENTS
        corner_frequency = corner * decade
        ratio = frequency / corner_frequency
        denominator = one(frequency) + ratio^2
        typed_coefficient = convert(T, coefficient)
        permittivity_sum += typed_coefficient / denominator
        conductivity_sum += typed_coefficient * ratio / denominator
        decade *= ten
    end
    relative_permittivity = convert(T, values.epsilon_infinity) + permittivity_sum
    conductivity = conductivity_reference +
                   2 * (one(frequency) * π) * frequency *
                   vacuum_permittivity(frequency) * conductivity_sum
    return EarthMaterial{T}(inv(conductivity), relative_permittivity, material.mu_r)
end

# Return the stable discovery identifier.
:LongmireSmith1975
