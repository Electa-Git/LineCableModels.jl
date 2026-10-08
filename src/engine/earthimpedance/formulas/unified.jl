"""
$(TYPEDSIGNATURES)

**Identification.** Unified circumferentially averaged earth impedance for overhead, buried
and mixed conductor systems in two homogeneous half-spaces, with complete enclosed-current
normalization.

**Expression.** `EarthAdmittance.source_coefficients` computes the axial-field coefficient of
each conductor pair, together with its source-potential coefficient, for the Unified formula
of either earth family.

**Validity.** The formula holds within two ranges.
- A nonzero prescribed Γ lies between the propagation constants of air and earth. It is in
  range at a frequency where Im γ_air ≤ Im Γ ≤ Im γ_earth and Re Γ ≤ Re γ_earth. Γ = 0, the
  default, always lies in range.
- Each receiving conductor is represented by the mean field on its exterior circumference,
  which holds while |κ_m r_p| ≤ 0.1. Here κ_m is the outgoing root of γ_m² − Γ² in the
  medium of conductor p, and r_p is its exterior radius, the jacket for an insulated cable.
  Above the range, the voltage that a thick conductor receives departs from the full-field
  value, with an error that grows as (κ_m r_p)². Emission from a thick conductor is
  accurate.
"""
function description(::Type{<:Formula{:unified}}; compact::Bool = false)
    compact ? "Unified" :
    "Unified circumferential earth impedance with complete enclosed-current normalization"
end

function description(::Type{<:Formula{:unified}}, ::Val{:Γ}, value::Number; compact::Bool = false)
    "Γ="*string(value)*" m⁻¹"
end
function description(::Type{<:Formula{:unified}}, ::Val{:Γ}, value::AbstractVector; compact::Bool = false)
    "Γ=["*join(value, ", ")*"] m⁻¹ (frequency order)"
end

:unified
