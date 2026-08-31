"""
$(TYPEDEF)

Represent a layered earth by one selected physical soil layer. This is an
explicit reconstruction policy, not a literature formula.

$(TYPEDFIELDS)
"""
struct Layer <: AbstractRule
    "Earth-model layer index; `-1` selects the bottommost layer."
    layer::Int

    function Layer(layer::Int)
        layer == -1 || layer >= 2 || throw(DomainError(
            layer,
            "EHEM layer must be -1 or an earth-model index of at least 2"
        ))
        return new(layer)
    end
end

description(rule::Layer) = rule.layer == -1 ?
                           "Bottommost earth layer" :
                           "Earth layer $(rule.layer)"

@inline function (rule::Layer)(
        ::Val,
        rho::AbstractVector,
        eps_r::AbstractVector,
        mu_r::AbstractVector,
        model::EarthModel,
        pair,
        frequency::Real
)
    _check(rho, eps_r, mu_r, model)
    return _material(rho, eps_r, mu_r, model, rule.layer)
end
