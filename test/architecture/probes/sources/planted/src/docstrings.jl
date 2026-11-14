struct Layer
    r::Float64
    @doc """
        Layer(r)

    Build a layer of radius `r`.
    """ function Layer(r::Real)
        return new(r)
    end
end
"""
    Layer(r, t)

Build a layer of radius `r + t`.
"""
Layer(r, t) = Layer(r + t)
"""
Build a layer of unit radius. $(TYPEDSIGNATURES)
"""
Layer() = Layer(1.0)
struct Wire{T}
    r::T
end
"""
    Wire{T}(r, n)

Build a wire of radius `r / n`.
"""
function Wire{T}(r, n) where {T}
    return Wire{T}(r / n)
end
