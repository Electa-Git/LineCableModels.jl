import ..LineCableModels: preview

"""
$(TYPEDSIGNATURES)

Generates a color representation for a [`Material`](@ref) based on its physical properties.

# Arguments

- `material_props`: Dictionary containing material properties:
  - `rho`: Electrical resistivity \\[Ω·m\\].
  - `eps_r`: Relative permittivity \\[dimensionless\\].
  - `mu_r`: Relative permeability \\[dimensionless\\].
- `rho_weight`: Weight assigned to resistivity in color blending (default: 1.0) \\[dimensionless\\].
- `epsr_weight`: Weight assigned to permittivity in color blending (default: 0.1) \\[dimensionless\\].
- `mur_weight`: Weight assigned to permeability in color blending (default: 0.1) \\[dimensionless\\].

# Returns

- An `RGBA` object representing the combined color based on the material's properties.

# Notes

Colors are normalized and weighted using property-specific gradients:
- Conductors (ρ ≤ 5ρ₀): White → Dark gray
- Poor conductors (5ρ₀ < ρ ≤ 10⁴): Bronze → Greenish-brown
- Insulators (ρ > 10⁴): Greenish-brown → Black
- Permittivity: Gray → Orange
- Permeability: Silver → Purple
- The overlay function combines colors with their respective alpha/weight values.

# Examples

```julia
material_props = Dict(
	:rho => 1.7241e-8,
	:eps_r => 2.3,
	:mu_r => 1.0
)
color = $(FUNCTIONNAME)(material_props)
println(color) # Expected output: RGBA(0.9, 0.9, 0.9, 1.0)
```
"""
function _get_material_color(
    material_props;
    rho_weight=1.0, #0.8,
    epsr_weight=0.1,
    mur_weight=0.1,
)

    # Auxiliar function to combine colors
    function _overlay_colors(colors::Vector{<:RGBA})
        # Handle edge cases
        if length(colors) == 0
            return RGBA(0, 0, 0, 0)
        elseif length(colors) == 1
            return colors[1]
        end

        # Initialize with the first color
        r, g, b, a = red(colors[1]), green(colors[1]), blue(colors[1]), alpha(colors[1])

        # Single-pass overlay for the remaining colors
        for i in 2:length(colors)
            r2, g2, b2, a2 =
                red(colors[i]), green(colors[i]), blue(colors[i]), alpha(colors[i])
            a_new = a2 + a * (1 - a2)

            if a_new == 0
                r, g, b, a = 0, 0, 0, 0
            else
                r = (r2 * a2 + r * a * (1 - a2)) / a_new
                g = (g2 * a2 + g * a * (1 - a2)) / a_new
                b = (b2 * a2 + b * a * (1 - a2)) / a_new
                a = a_new
            end
        end

        return RGBA(r, g, b, a)
    end

    # Fixed normalization bounds
    epsr_min, epsr_max = 1.0, 1000.0  # Adjusted permittivity range for semiconductors
    mur_min, mur_max = 1.0, 300.0  # Relative permeability range
    rho_base = 1.72e-8

    # Extract nominal values for uncertain measurements
    rho = to_nominal(material_props.rho)
    epsr_r = to_nominal(material_props.eps_r)
    mu_r = to_nominal(material_props.mu_r)

    # Handle air/void
    if isinf(rho)
        return RGBA(1.0, 1.0, 1.0, 1.0)  # Transparent white
    end

    # Normalize epsr and mur
    epsr_norm = (epsr_r - epsr_min) / (epsr_max - epsr_min)
    mur_norm = (mu_r - mur_min) / (mur_max - mur_min)

    # Define color gradients based on resistivity
    if rho <= 5 * rho_base
        # Conductors: Bright metallic white → Darker metallic gray (logarithmic scaling)
        rho_norm = log10(rho / rho_base) / log10(5)  # Normalize based on `5 * rho_base`

        rho_color = get(cgrad([
                RGB(0.9, 0.9, 0.9),  # Almost white
                RGB(0.6, 0.6, 0.6),  # Light gray
                RGB(0.4, 0.4, 0.4)  # Dark gray
            ]), clamp(rho_norm, 0.0, 1.0))

    elseif rho <= 10000
        # Poor conductors/semiconductors: Bronze → Gold → Reddish-brown → Dark orange → Greenish-brown
        rho_norm = (rho - 10e-8) / (10000 - 10e-8)
        rho_color = get(
            cgrad([
                RGB(0.8, 0.5, 0.2),  # Metallic bronze
                RGB(1.0, 0.85, 0.4),  # Metallic gold
                RGB(0.8, 0.4, 0.2),  # Reddish-brown
                RGB(0.8, 0.3, 0.1),  # Dark orange
                RGB(0.6, 0.4, 0.3),   # Greenish-brown
            ]), rho_norm)

    else
        # Insulators: Greenish-brown → Black
        rho_norm = (rho - 10000) / (1e5 - 10000)
        rho_color = get(cgrad([RGB(0.6, 0.4, 0.3), :black]), clamp(rho_norm, 0.0, 1.0))
    end

    # Normalize epsr and mur values to [0, 1]
    epsr_norm = clamp(epsr_norm, 0.0, 1.0)
    mur_norm = clamp(mur_norm, 0.0, 1.0)

    # Create color gradients for epsr and mur
    epsr_color = get(cgrad([:gray, RGB(1.0, 0.9, 0.7), :orange]), epsr_norm)  # Custom amber
    mur_color = get(
        cgrad([:silver, :gray, RGB(0.9, 0.8, 1.0), :purple, RGB(0.3, 0.1, 0.6)]),
        mur_norm,
    )  # Custom purple

    # Apply weights to each property
    rho_color_w = Colors.RGBA(rho_color.r, rho_color.g, rho_color.b, rho_weight)
    epsr_color_w = Colors.RGBA(epsr_color.r, epsr_color.g, epsr_color.b, epsr_weight)
    mur_color_w = Colors.RGBA(mur_color.r, mur_color.g, mur_color.b, mur_weight)

    # Combine weighted colors
    final_color = _overlay_colors([rho_color_w, epsr_color_w, mur_color_w])

    return final_color
end

"""
$(TYPEDSIGNATURES)

Displays the cross-section of a cable design.

# Arguments

- `design`: A [`CableDesign`](@ref) object representing the cable structure.
- `x_offset`: Horizontal offset for the plot \\[m\\].
- `y_offset`: Vertical offset for the plot \\[m\\].
- `plt`: An optional `Plots.Plot` object to use for plotting.
- `display_plot`: Boolean flag to display the plot after rendering.
- `display_legend`: Boolean flag to display the legend in the plot.
- `backend`: Optional plotting backend to use. If not specified, the function will choose a suitable backend based on the environment (e.g., GR for headless, PlotlyJS for interactive).
- `sz`: Optional plot dimensions (width, height). Default: (800, 600).

# Returns

- A `Plots.Plot` object representing the visualized cable design.

# Examples

```julia
conductor_group = ConductorGroup(central_conductor)
insulator_group = InsulatorGroup(main_insulation)
component = CableComponent("core", conductor_group, insulator_group)
design = CableDesign("example", component)
cable_plot = $(FUNCTIONNAME)(design)  # Cable cross-section is displayed
```

# See also

- [`CableDesign`](@ref)
- [`ConductorGroup`](@ref)
- [`InsulatorGroup`](@ref)
- [`WireArray`](@ref)
- [`Tubular`](@ref)
- [`Strip`](@ref)
- [`Semicon`](@ref)
"""
function preview(
    design::CableDesign;
    x_offset=0.0,
    y_offset=0.0,
    plt=nothing,
    display_plot=true,
    display_legend=true,
    backend=nothing,
    sz=(800, 600),
)
    if isnothing(plt)
        # Choose appropriate backend based on environment
        resolve_backend(backend)
        plt = plot(size=sz,
            aspect_ratio=:equal,
            legend=(0.875, 1.0),
            title="Cable design preview",
            xlabel="y [m]",
            ylabel="z [m]")
    end

    # Helper function to plot a layer
    function _plot_layer!(layer, label; x0=0.0, y0=0.0)
        if layer isa WireArray
            radius_wire = to_nominal(layer.radius_wire)
            num_wires = layer.num_wires

            lay_radius = num_wires == 1 ? 0.0 : to_nominal(layer.radius_in)
            material_props = layer.material_props
            color = _get_material_color(material_props)

            # Use the existing calc_wirearray_coords function to get wire centers
            wire_coords = calc_wirearray_coords(
                num_wires,
                radius_wire,
                to_nominal(lay_radius),
                C=(x0, y0),
            )

            # Plot each wire in the layer
            for (i, (x, y)) in enumerate(wire_coords)
                plot!(
                    plt,
                    Shape(
                        x .+ radius_wire * cos.(0:0.01:2π),
                        y .+ radius_wire * sin.(0:0.01:2π),
                    ),
                    linecolor=:black,
                    color=color,
                    label=(i == 1 && display_legend) ? label : "",  # Only add label for first wire
                )
            end

        elseif layer isa Strip || layer isa Tubular || layer isa Semicon ||
               layer isa Insulator
            radius_in = to_nominal(layer.radius_in)
            radius_ext = to_nominal(layer.radius_ext)
            material_props = layer.material_props
            color = _get_material_color(material_props)

            arcshape(θ1, θ2, rin, rext, x0=0.0, y0=0.0, N=100) = begin
                # Generate angles for the arc
                θ = range(θ1, θ2, length=N)

                # Outer circle coordinates
                x_outer = x0 .+ rext .* cos.(θ)
                y_outer = y0 .+ rext .* sin.(θ)

                # Inner circle coordinates (reversed to close the shape properly)
                x_inner = x0 .+ rin .* cos.(reverse(θ))
                y_inner = y0 .+ rin .* sin.(reverse(θ))

                # Concatenate and explicitly close the shape by repeating the first point
                x_coords = [x_outer; x_inner; x_outer[1]]
                y_coords = [y_outer; y_inner; y_outer[1]]

                Shape(x_coords, y_coords)
            end
            shape = arcshape(0, 2π, radius_in, radius_ext, x0, y0)
            plot!(
                plt,
                shape,
                linecolor=color,
                color=color,
                label=display_legend ? label : "",
            )
        end
    end

    # Iterate over all CableComponents in the design
    for component in design.components
        # Process conductor group layers
        for layer in component.conductor_group.layers
            # Check if layer is a compound structure
            if layer isa ConductorGroup
                # Special handling for nested conductor groups
                first_layer = true
                for sublayer in layer.layers
                    _plot_layer!(
                        sublayer,
                        first_layer ? lowercase(string(typeof(layer))) : "",
                        x0=x_offset,
                        y0=y_offset,
                    )
                    first_layer = false
                end
            else
                # Plot standard conductor layer
                _plot_layer!(
                    layer,
                    lowercase(string(typeof(layer))),
                    x0=x_offset,
                    y0=y_offset,
                )
            end
        end

        # Process insulator group layers
        for layer in component.insulator_group.layers
            _plot_layer!(
                layer,
                lowercase(string(typeof(layer))),
                x0=x_offset,
                y0=y_offset,
            )
        end
    end

    if display_plot
        if !is_in_testset()
            if is_headless()
                DisplayAs.Text(DisplayAs.PNG(plt))
            else
                display(plt)
            end
        end
    end

    return plt
end

"""
$(TYPEDSIGNATURES)

Displays the cross-section of a cable system.

# Arguments

- `system`: A [`LineCableSystem`](@ref) object containing the cable arrangement.
- `earth_model`: Optional [`EarthModel`](@ref) to display earth layers.
- `zoom_factor`: A scaling factor for adjusting the x-axis limits \\[dimensionless\\].
- `backend`: Optional plotting backend to use.
- `sz`: Optional plot dimensions (width, height). Default: (800, 600).

# Returns

- Nothing. Displays a plot of the cable system cross-section with cables, earth layers (if applicable), and the air/earth interface.

# Examples

```julia
system = LineCableSystem("test_system", 1000.0, cable_pos)
earth_params = EarthModel(f, 100.0, 10.0, 1.0)
$(FUNCTIONNAME)(system, earth_model=earth_params, zoom_factor=0.5)
```

# See also

- [`LineCableSystem`](@ref)
- [`EarthModel`](@ref)
- [`CablePosition`](@ref)
"""
function preview(
    system::LineCableSystem;
    earth_model=nothing,
    zoom_factor=0.25,
    backend=nothing,
    sz=(800, 600),
    display_plot=true,
)
    resolve_backend(backend)
    plt = plot(size=sz,
        aspect_ratio=:equal,
        legend=(0.8, 0.9),
        title="Cable system cross-section",
        xlabel="y [m]",
        ylabel="z [m]")

    # Plot the air/earth interface at y=0
    hline!(
        plt,
        [0],
        linestyle=:solid,
        linecolor=:black,
        linewidth=1.25,
        label="Air/earth interface",
    )

    # Determine explicit wide horizontal range for earth layer plotting
    x_positions = [to_nominal(cable.horz) for cable in system.cables]
    max_span = maximum(abs, x_positions) + 5  # extend 5 m beyond farthest cable position
    x_limits = [-max_span, max_span]

    # Plot earth layers if provided and vertical_layers == false
    if !isnothing(earth_model) && !earth_model.vertical_layers
        layer_colors = [:burlywood, :sienna, :peru, :tan, :goldenrod, :chocolate]
        cumulative_depth = 0.0
        for (i, layer) in enumerate(earth_model.layers[2:end])
            # Skip bottommost infinite layer
            if isinf(layer.t)
                break
            end

            # Compute the depth of the current interface
            cumulative_depth -= layer.t
            hline!(
                plt,
                [cumulative_depth],
                linestyle=:solid,
                linecolor=layer_colors[mod1(i, length(layer_colors))],
                linewidth=1.25,
                label="Earth layer $i",
            )

            # Fill the area for current earth layer
            y_coords = [cumulative_depth + layer.t, cumulative_depth]
            plot!(plt, [x_limits[1], x_limits[2], x_limits[2], x_limits[1]],
                [y_coords[1], y_coords[1], y_coords[2], y_coords[2]],
                seriestype=:shape, color=layer_colors[mod1(i, length(layer_colors))],
                alpha=0.25, linecolor=:transparent,
                label="")
        end
    end

    for cable_position in system.cables
        x_offset = to_nominal(cable_position.horz)
        y_offset = to_nominal(cable_position.vert)
        preview(
            cable_position.design_data;  # Changed from cable_position.design_data
            x_offset,
            y_offset,
            plt,
            display_plot=false,
            display_legend=false,
        )
    end

    plot!(plt, xlim=(x_limits[1], x_limits[2]) .* zoom_factor)

    if display_plot
        if !is_in_testset()
            if is_headless()
                DisplayAs.Text(DisplayAs.PNG(plt))
            else
                display(plt)
            end
        end
    end

    return plt
end