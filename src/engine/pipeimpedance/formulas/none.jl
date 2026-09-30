"""
$(TYPEDSIGNATURES)

**Identification.** No additional pipe impedance for coaxial geometry.

**Expression.** Ordinary coaxial assemblies do not require an additional pipe term.
An eccentric or multicore conductive enclosure requires a pipe formulation:
the coaxial formulation supplies none.

**Applicability.** This selector does not add an author-labeled pipe equation or numerical approximation.
"""
description(::Type{<:Formula{:none}}; compact::Bool=false) = compact ? "No additional pipe term" : "No additional pipe impedance for coaxial geometry"

function Formulation(::LineCableModelsCoaxial, ::Formula{:none}, ::Val{:coaxial})
    nothing
end

function Formulation(::LineCableModelsCoaxial, ::Formula{:none}, ::Val{:pipe})
    throw(ArgumentError(
        "Pipe-type cable formulation is not yet implemented for the coaxial backend. No pipe formulation is available."))
end

:none
