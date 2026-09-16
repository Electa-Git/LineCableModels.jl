TextDisplay.@showfields EarthMaterial "EarthMaterial" material -> (
    ρ = TextDisplay.engineering(material.rho, :ohm_meter),
    εᵣ = TextDisplay.value(material.eps_r),
    μᵣ = TextDisplay.value(material.mu_r)
)

TextDisplay.@showfields EarthLayer "EarthLayer" layer -> (
    ρ = TextDisplay.engineering(layer.rho, :ohm_meter),
    εᵣ = TextDisplay.value(layer.eps_r),
    μᵣ = TextDisplay.value(layer.mu_r),
    h = isinf(layer.thickness) ? "half-space" :
        TextDisplay.engineering(layer.thickness, :meter)
)

TextDisplay.name(::Type{<:EarthModel}) = "EarthModel"

function _earth_description(model::EarthModel)
    earth_count = length(model.layers) - 1
    earth_count == 1 && return "homogeneous"
    orientation = model.vertical_layers ? "vertical " : ""
    return "$earth_count $(orientation)earth layers"
end

function Base.summary(io::IO, model::EarthModel)
    print(io, "EarthModel · ", _earth_description(model))
end

function Base.show(io::IO, model::EarthModel)
    print(io, "EarthModel(", _earth_description(model), ")")
end

function _earth_layer_name(index::Int, count::Int)
    index == 1 && return "air"
    count == 2 && return "earth"
    index == count && return "basement"
    return "layer $(index - 1)"
end

function _earth_resistivity(layer::EarthLayer)
    return isinf(layer.rho) ? "∞" : TextDisplay.engineering(layer.rho, :ohm_meter)
end

function _earth_layer_text(
        name::AbstractString,
        layer::EarthLayer;
        show_thickness::Bool
)
    thickness = if !show_thickness
        ""
    elseif isinf(layer.thickness)
        "half-space"
    else
        "h=$(TextDisplay.engineering(layer.thickness, :meter))"
    end
    separator = show_thickness ? string("  ", rpad(thickness, 12)) : ""
    return string(
        rpad(name, 9), separator,
        "  ρ=", _earth_resistivity(layer),
        "  εᵣ=", TextDisplay.value(layer.eps_r),
        "  μᵣ=", TextDisplay.value(layer.mu_r)
    )
end

function Base.show(io::IO, ::MIME"text/plain", model::EarthModel)
    get(io, :compact, false) && return show(io, model)
    show_thickness = length(model.layers) > 2
    children = [
        _earth_layer_text(
            _earth_layer_name(index, length(model.layers)),
            layer;
            show_thickness
        )
        for (index, layer) in enumerate(model.layers)
    ]
    return TextDisplay.tree(
        io,
        "EarthModel · $(_earth_description(model))",
        children;
        noun = "layers"
    )
end
