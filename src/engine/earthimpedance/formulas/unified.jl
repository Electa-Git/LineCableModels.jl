"""
$(TYPEDSIGNATURES)

**Identification.** Unified circumferentially averaged earth impedance for overhead, buried
and mixed conductor systems in two homogeneous half-spaces, with complete enclosed-current
normalization.

**Expression.** `EarthAdmittance.source_coefficients` computes the axial-field coefficient of
each conductor pair, together with its source-potential coefficient, for the Unified formula
of either earth family.
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
