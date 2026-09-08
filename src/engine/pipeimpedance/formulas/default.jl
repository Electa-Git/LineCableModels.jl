"""
$(TYPEDSIGNATURES)

**Identification.** Default analytical pipe-type treatment.

**Expression.** Ordinary coaxial assemblies require no additional pipe term.
An eccentric or multicore conductive enclosure requires a pipe formulation:
the coaxial backend has none yet.

**Reference.** Backend applicability policy; no author-labelled pipe equation
or numerical approximation is introduced by this selector.
"""
description(::Formula{:default}) = "Default analytical pipe-type treatment"

function Formulation(::LineCableModelsCoaxial, ::Val{:default}, ::Formula{:default}, ::Val{:coaxial})
    nothing
end

function Formulation(::LineCableModelsCoaxial, ::Val{:default}, ::Formula{:default}, ::Val{:pipe})
    throw(ArgumentError(
        "Pipe-type cable formulation is not yet implemented for the coaxial backend. No default formulation is available."))
end

:default
