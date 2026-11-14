"""
$(TYPEDEF)

A layer of radius `r`.
"""
struct Layer
    r::Float64
    @doc """
    $(TYPEDSIGNATURES)

    Build a layer of radius `r`.
    """ function Layer(r::Real)
        return new(r)
    end
end
"""
$(SIGNATURES)

Build a layer of radius `r + t`.
"""
Layer(r, t) = Layer(r + t)
"""
    measure(layer)

Return the radius of `layer`.
"""
function measure end
"""
    measure(layer)

Return the radius of `layer`.
"""
measure(layer::Layer) = layer.r
