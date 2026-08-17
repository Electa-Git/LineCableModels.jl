function Base.show(io::IO, ::MIME"text/plain", model::EarthModel)
    earth_count = length(model.layers) - 1
    orientation = model.vertical_layers ? "vertical" : "horizontal"
    kind = earth_count == 1 ? "homogeneous" : "multilayer"
    println(io,
        "EarthModel with $earth_count $orientation earth " *
        "$(earth_count == 1 ? "layer" : "layers") ($kind)")
    for (index, layer) in enumerate(model.layers)
        prefix = index == length(model.layers) ? "└─" : "├─"
        name = index == 1 ? "Layer 1 (air)" : "Layer $index"
        thickness = isinf(layer.thickness) ? "Inf" :
                    string(round(layer.thickness; sigdigits = 4))
        println(
            io,
            "$prefix $name: [rho=$(round(layer.rho; sigdigits=4)), " *
            "eps_r=$(round(layer.eps_r; sigdigits=4)), " *
            "mu_r=$(round(layer.mu_r; sigdigits=4)), thickness=$thickness]"
        )
    end
end
