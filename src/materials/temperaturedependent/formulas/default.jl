"""
$(TYPEDSIGNATURES)

**Identification.** Linear electrical-resistivity temperature dependence.

**Expression.** ``\\rho(T)=\\rho_0[1+\\alpha(T-T_0)]``, with resistivity in Ω·m,
temperature in °C and ``\\alpha`` in K⁻¹. Each reference material supplies its
own calibration. Infinite passive resistivity remains infinite.

**Applicability.** The existing model policy requires ``|T-T_0|<150`` K and a
finite positive correction factor. This is a restriction of this approximation,
not a thermal-rating or material operating-temperature limit.

**Reference.** Package default linear resistivity approximation.
"""
description(::Formula{:default}) = "Linear electrical-resistivity temperature dependence"

function temperature_resistivity(::Val{:default}, material::Material, temperature::Real,
        parameters::NamedTuple, options::NamedTuple, workspace)
    difference = temperature - material.T0
    abs(difference) < oftype(difference, 150) || throw(DomainError(temperature,
        "temperature is outside the linear resistivity model range relative to $(material.T0) °C"))
    factor = one(difference) + material.alpha * difference
    isfinite(factor) && factor > zero(factor) || throw(DomainError(factor,
        "linear resistivity correction factor must be positive and finite"))
    return isinf(material.rho) ? material.rho : material.rho * factor
end

computation_options(::FormulaMethod{:default, typeof(temperature_resistivity)}) = (;)

:default
