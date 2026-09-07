routes(::Val{:DirectNumericalIntegration}) = (;)
assumptions(::Val{:DirectNumericalIntegration}) = (;)
propagation(::Val{:DirectNumericalIntegration}) = Val(:backend)

"""
$(TYPEDSIGNATURES)

**Identification.** Direct numerical integration selected in an external
earth-return backend.

**Expression.** The backend evaluates its earth-return integral numerically.
The selected placement determines the integral and native solver setting.
No analytical kernel is supplied by this catalogue entry.

**Reference.** Backend numerical-integration mode, not an author-labelled
analytical approximation. PSCAD exposes this setting separately for overhead
and underground lines; the adapter records the selected native setting.
"""
description(::Formula{:DirectNumericalIntegration}) = "Direct numerical earth-return integration"

:DirectNumericalIntegration
