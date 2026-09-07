assumptions(::Val{:default}) = (;)

"""
$(TYPEDSIGNATURES)

**Identification.** Backend-dependent default pipe-type treatment.

**Expression.** Ordinary coaxial assemblies require no additional pipe term.
An eccentric or multicore conductive enclosure requires a pipe formulation:
the coaxial backend has none yet. FEM resolves supported enclosures directly
with its field equations, without an additional analytical pipe correction.

**Reference.** Backend applicability policy; no author-labelled pipe equation
or numerical approximation is introduced by this selector.
"""
description(::Formula{:default}) = "Default backend pipe-type treatment"

Formulation(::LineCableModelsCoaxial, ::Val{:default}, ::Formula{:default}, ::Val{:coaxial}) = nothing

function Formulation(::LineCableModelsCoaxial, ::Val{:default}, ::Formula{:default}, ::Val{:pipe})
    throw(ArgumentError(
        "Pipe-type cable formulation is not yet implemented for the coaxial backend. No default formulation is available."))
end

# FEM owns the enclosure surfaces and their material and terminal domains.
Formulation(::LineCableModelsFEM, ::Formula{:default}, ::CableDesign) = nothing

:default
