"""
$(TYPEDEF)

Reduce evaluated layered-earth properties to air plus one selected layer.

$(TYPEDFIELDS)
"""
struct EnforceLayer <: AbstractEHEMFormulation
    "Earth-layer index. `-1` selects the bottom layer."
    layer::Int

    function EnforceLayer(; layer::Int = -1)
        layer == -1 || layer >= 2 ||
            throw(DomainError(
                layer, "layer must be -1 or an earth-layer index of at least 2"
            ))
        return new(layer)
    end
end

function description(formulation::EnforceLayer)
    formulation.layer == -1 ?
    "Assume bottom layer" :
    formulation.layer == 2 ?
    "Assume top earth layer" : "Assume layer $(formulation.layer)"
end

function (formulation::EnforceLayer)(evaluated::NamedTuple, model::EarthModel)
    selected = formulation.layer == -1 ? length(model.layers) : formulation.layer
    selected in 2:length(model.layers) || throw(BoundsError(model.layers, selected))
    rows = [1, selected]
    return (
        rho = evaluated.rho[rows, :],
        epsilon = evaluated.epsilon[rows, :],
        mu = evaluated.mu[rows, :]
    )
end
