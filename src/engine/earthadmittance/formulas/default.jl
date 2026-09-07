function Formula(::Val{:default}; kwargs...)
    overrides = (; kwargs...)
    return Formula{:default, typeof(overrides), @NamedTuple{}}(overrides, (;))
end

function Formula(selected::Formula{:default}, ::Val{:overhead})
    resolved = Formula(Val(:Wise1948); selected.routes...)
    return Formula(Val(:Wise1948), resolved.routes,
        merge(resolved.assumptions, selected.assumptions))
end

function Formula(selected::Formula{:default}, ::Val{:underground})
    resolved = Formula(Val(:Xue2018); selected.routes...)
    return Formula(Val(:Xue2018), resolved.routes,
        merge(resolved.assumptions, selected.assumptions))
end

function Formula(::Formula{:default}, ::Val{Placement}) where {Placement}
    throw(ArgumentError(
        "earth-admittance :default is not yet implemented for $Placement placement; " *
        "select an explicit applicable formulation"))
end

"""
$(TYPEDSIGNATURES)

**Identification.** Problem-dependent default earth-admittance selection.

**Expression.** Select `:Wise1948` for overhead placement and `:Xue2018`
for underground placement. Resolve the selection before numerical evaluation;
retain explicit routes and the existing FD/EHEM policy. Mixed placement has no
default implementation.

**Reference.** The selected Wise or Xue formula owns its equations and
literature reference. This selector introduces no additional numerical formula.
"""
description(::Formula{:default}) = "Select earth admittance from the resolved problem context"

:default
