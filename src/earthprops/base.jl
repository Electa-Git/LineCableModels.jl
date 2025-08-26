import Base: show, eltype, convert

eltype(::EarthModel{T}) where {T} = T

function convert(::Type{EarthModel{T}}, model::EarthModel) where {T}
    # If the model is already the target type, return it without modification.
    model isa EarthModel{T} && return model

    # Delegate the actual conversion logic to the existing coerce_to_T function.
    return coerce_to_T(model, T)
end

function convert(::Type{EarthLayer{T}}, layer::EarthLayer) where {T}
    # Avoid unnecessary work if the layer is already the correct type.
    layer isa EarthLayer{T} && return layer

    # Delegate the conversion logic to the specialized coerce_to_T function.
    return coerce_to_T(layer, T)
end

"""
$(TYPEDSIGNATURES)

Defines the display representation of a [`EarthModel`](@ref) object for REPL or text output.

# Arguments

- `io`: The output stream to write the representation to \\[IO\\].
- `mime`: The MIME type for plain text output \\[MIME"text/plain"\\].
- `model`: The [`EarthModel`](@ref) instance to be displayed.


# Returns

- Nothing. Modifies `io` to format the output.
"""
function show(io::IO, ::MIME"text/plain", model::EarthModel)
    # Determine model type based on num_layers and vertical_layers flag
    num_layers = length(model.layers)
    model_type = num_layers == 2 ? "homogeneous" : "multilayer"
    orientation = model.vertical_layers ? "vertical" : "horizontal"
    layer_word = (num_layers - 1) == 1 ? "layer" : "layers"

    # Count frequency samples from the first layer's property arrays
    num_freq_samples = length(model.layers[1].rho_g)
    freq_word = (num_freq_samples) == 1 ? "sample" : "samples"

    # Print header with key information
    println(
        io,
        "EarthModel with $(num_layers-1) $(orientation) earth $(layer_word) ($(model_type)) and $(num_freq_samples) frequency $(freq_word)",
    )

    # Print layers in treeview style
    for i in 1:num_layers
        layer = model.layers[i]
        # Determine prefix based on whether it's the last layer
        prefix = i == num_layers ? "└─" : "├─"

        # Format thickness value
        thickness_str = isinf(layer.t) ? "Inf" : "$(round(layer.t, sigdigits=4))"

        # Format layer name
        layer_name = i == 1 ? "Layer $i (air)" : "Layer $i"

        # Print layer properties with proper formatting
        println(
            io,
            "$prefix $layer_name: [rho_g=$(round(layer.base_rho_g, sigdigits=4)), " *
            "epsr_g=$(round(layer.base_epsr_g, sigdigits=4)), " *
            "mur_g=$(round(layer.base_mur_g, sigdigits=4)), " *
            "t=$thickness_str]",
        )
    end

    # Add formulation information as child nodes
    if !isnothing(model.freq_dependence)
        formulation_tag = get_description(model.freq_dependence)
        println(io, "   Frequency-dependent model: $(formulation_tag)")
    end

end