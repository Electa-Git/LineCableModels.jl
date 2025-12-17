module plot_runner

using Makie
import Makie: plot
using Colors
using StatsBase: fit, Histogram, normalize
using Statistics: quantile
using Base: basename
using Dates
using Printf: @sprintf

using LineCableModels.Engine: LP_FIG_SIZE, _ICON_FN
import LineCableModels.BackendHandler: BackendHandler
using LineCableModels.PlotUIComponents: ICON_TTF, _make_window,	with_plot_theme
using LineCableModels.UQ: LineParametersMCSummary, _hist_specs

export plot_frequency_ridgeline, 
       plot_frequency_fanchart, 
       plot_spectral_comparison, 
       plot_multi_voltage_ridgeline_grid,
       plot_multi_voltage_fanchart_grid

# --- Constants & Helpers ---
const SEQUENCE_MAP = Dict(
    :zero     => ((1, 1), "Zero Sequence"),
    :positive => ((2, 2), "Positive Sequence"),
    :negative => ((3, 3), "Negative Sequence")
)

"""
    CablePlotItem
Simple wrapper to carry geometry metadata alongside the simulation result.
"""
struct CablePlotItem
    obj::LineParametersMCSummary
    conductor::Int
    screen::Int
end

function _get_seq_info(seq::Symbol)
    if !haskey(SEQUENCE_MAP, seq)
        throw(ArgumentError("Invalid sequence: $seq. Must be :zero, :positive, or :negative"))
    end
    return SEQUENCE_MAP[seq]
end

"""
    _compute_ridgeline_data(obj, quantity, indices, seq_indices; ...)

Internal helper to compute histogram densities, normalize, and trim tails for ridgeline plots.
Returns a tuple: `(densities, global_max_dens, x_label)`
"""
function _compute_ridgeline_data(
    obj::LineParametersMCSummary,
    quantity::Symbol,
    indices::Vector{Int},
    (i_seq, j_seq)::Tuple{Int, Int};
    length_unit::Symbol = :kilo,
    per_length::Bool = true,
    quantity_units = nothing,
    nbins::Int = 80,
    scale_factor::Float64 = 1.0,
    step_style::Bool = false
)

    densities = Vector{Tuple{Int, Any, Vector{Float64}, Vector{Float64}}}()
    global_max_dens = 0.0
    xlabel = ""

    for k in indices
        # Fetch histogram spec
        spec = first(_hist_specs(
            obj, quantity; 
            ijk = (i_seq, j_seq, k), 
            length_unit = length_unit,
            per_length = per_length, 
            quantity_units = quantity_units, 
            nbins = nbins, 
            mode = :hist
        ))
        
        # Save label from first iteration
        if isempty(xlabel)
            xlabel = spec.xlabel
        end

        # Scale values
        vals = spec.values .* scale_factor
        
        # Dynamic binning (recalc based on scaled values)
        min_v, max_v = extrema(vals)
        margin = max((max_v - min_v) * 0.1, 1e-6)
        bins_use = range(min_v - margin, max_v + margin, length=nbins+1)

        h = fit(Histogram, vals, bins_use; closed=:left)
        h_norm = normalize(h, mode = :pdf)
        
        edges = h_norm.edges[1]
        weights = h_norm.weights

        # Trim Tails (Visual cleanup)
        cutoff = maximum(weights) * 0.005
        first_nz = findfirst(>(cutoff), weights)
        last_nz  = findlast(>(cutoff), weights)
        
        if !isnothing(first_nz) && !isnothing(last_nz)
            start_idx = max(1, first_nz - 1)
            stop_idx  = min(length(weights), last_nz + 1)
            weights = weights[start_idx:stop_idx]
            edges   = edges[start_idx:stop_idx+1]
        end

        # Format for Plotting (Step vs Smooth) - not used atm
        if step_style
            x_plot = Float64[]
            y_plot = Float64[]
            sizehint!(x_plot, 2*length(weights) + 2)
            sizehint!(y_plot, 2*length(weights) + 2)
            for j in 1:length(weights)
                push!(x_plot, edges[j], edges[j+1])
                push!(y_plot, weights[j], weights[j])
            end
        else
            x_plot = (edges[1:end-1] .+ edges[2:end]) ./ 2
            y_plot = weights
        end
        
        if !isempty(y_plot)
            current_max = maximum(y_plot)
            global_max_dens = max(global_max_dens, current_max)
            push!(densities, (k, spec, x_plot, y_plot))
        end
    end

    return densities, global_max_dens, xlabel
end

"""
    _compute_quantile_data(obj, quantity, sequence)

Internal helper to extract 5%, 25%, 50%, 75%, 95% quantiles across all frequencies.
Returns: (freqs, q_matrix) where q_matrix is 5 x N_freq.
"""
function _compute_quantile_data(
    obj::LineParametersMCSummary, 
    quantity::Symbol, 
    sequence::Symbol
)
    (ij, _) = _get_seq_info(sequence)
    
    # Scale settings
    scale_factor = (quantity == :L) ? 1000.0 : 1.0 # H -> mH
    
    # Data extraction
    raw_samples = getfield(obj.samples, quantity)
    freqs = obj.f
    n_freqs = length(freqs)
    
    # 5 quantiles x N frequencies
    qs = zeros(Float64, 5, n_freqs) 

    for k in 1:n_freqs
        vals = raw_samples[ij[1], ij[2], k, :] .* scale_factor
        qs[:, k] .= quantile(vals, [0.05, 0.25, 0.50, 0.75, 0.95])
    end
    
    return freqs, qs
end

# --- Main Plotting Functions ---

"""
    plot_frequency_ridgeline(obj, quantity, indices; ...)
    
Generates a standalone ridgeline plot showing the evolution of uncertainty over frequency.
"""
function plot_frequency_ridgeline(
    obj::LineParametersMCSummary,
    quantity::Symbol,
    indices::Vector{Int};
    sequence::Symbol = :positive,
    length_unit::Symbol = :kilo,
    per_length::Bool = true,
    fig_size::Tuple{Int, Int} = (800, 600),
    quantity_units = nothing,
    nbins::Int = 80,
    overlap::Float64 = 0.5,
    alpha::Float64 = 0.9,
    step_style::Bool = false,
    colormap::Symbol = :Blues,
    backend = nothing
)
    # Setup
    (ij, seq_title) = _get_seq_info(sequence)
    full_title = "Evolution of $(seq_title) $(quantity) over Frequency"
    
    # Compute Data
    densities, global_max_dens, x_lbl = _compute_ridgeline_data(
        obj, quantity, indices, ij;
        length_unit=length_unit, per_length=per_length, quantity_units=quantity_units,
        nbins=nbins, step_style=step_style
    )

    # Figure Setup
    backend_ctx = _make_window(BackendHandler, backend; title = full_title, icons = _ICON_FN, icons_font = ICON_TTF)
    fig = Makie.Figure(size = fig_size)
    ax = Makie.Axis(fig[1, 1], 
        xlabel = x_lbl, 
        ylabel = "Frequency [Hz] (Stacked)",
        title = full_title
    )

    # Plotting Loop
    y_step = global_max_dens * (1.0 - overlap)
    cmap = cgrad(colormap, length(indices), categorical = true)

    # Reverse order to draw top-down (back-to-front painter's algorithm)
    for (i, item) in enumerate(reverse(densities))
        (real_k, spec, x, y) = item
        
        # "i" here is the index in the reversed list, we need strict vertical stacking order
        stack_idx = length(densities) - i 
        base_y = stack_idx * y_step
        top_y = base_y .+ y
        
        freq_val = obj.f[real_k]
        color_idx = length(densities) - i + 1 # Align color with original frequency order

        band!(ax, x, fill(base_y, length(x)), top_y; color = (cmap[color_idx], alpha))
        lines!(ax, x, top_y; color = :black, linewidth = 1.0)
        
        # Label each ridge
        text!(ax, minimum(x), base_y + (y_step * 0.2); 
            text = @sprintf("%.1f Hz", freq_val), 
            align = (:left, :bottom), fontsize = 12
        )
    end

    hideydecorations!(ax, label = false, ticklabels = true, ticks = false, grid = false)

    if backend_ctx.interactive && !isnothing(backend_ctx.window)
        display(backend_ctx.window, fig)
    else
        BackendHandler.renderfig(fig)
    end
    
    return fig
end

"""
    _draw_ridgeline_on_axis!(axis, obj, quantity, indices; ...)
    
Internal function to draw ridgelines on a pre-existing axis (used for grids).
"""
function _draw_ridgeline_on_axis!(
    axis,
    obj::LineParametersMCSummary,
    quantity::Symbol,
    indices::Vector{Int};
    sequence::Symbol = :pos,
    scale_factor::Float64 = 1.0,
    length_unit::Symbol = :kilo,
    per_length::Bool = true,
    nbins::Int = 80,
    overlap::Float64 = 0.5,
    alpha::Float64 = 0.9,
    base_color::Colorant = colorant"blue",
)
    (ij, _) = _get_seq_info(sequence)

    # Compute Data
    densities, global_max_dens, _ = _compute_ridgeline_data(
        obj, quantity, indices, ij;
        length_unit=length_unit, per_length=per_length, nbins=nbins, scale_factor=scale_factor
    )

    if isempty(densities) return end

    # Plotting
    y_step = global_max_dens * (1.0 - overlap)
    n_plots = length(densities)

    # Create gradient from white -> base_color
    c_light = weighted_color_mean(0.7, colorant"white", base_color)
    cmap = cgrad([c_light, base_color], n_plots, categorical = true)

    for (i, (real_k, spec, x, y)) in enumerate(densities)
        base_y = (i - 1) * y_step
        top_y = base_y .+ y
        
        band!(axis, x, fill(base_y, length(x)), top_y; color = (cmap[i], alpha))
        lines!(axis, x, top_y; color = :black, linewidth = 0.5)
    end
end

"""
    plot_multi_voltage_ridgeline_grid(...)
    
Creates a giant grid of ridgeline plots (Row=Voltage, Col=Size).
"""
function plot_multi_voltage_ridgeline_grid(
    cable_lookup::AbstractDict, 
    sorted_volts::Vector{String}, 
    sorted_sizes::Vector{String}, 
    quantity::Symbol,
    indices::Vector{Int};
    sequence::Symbol = :positive,
    scale_factor::Float64 = 1.0,
    y_lbl::String = "Distribution Density",
    main_title::String = "",
)
    n_rows, n_cols = length(sorted_volts), length(sorted_sizes)
    f = Figure(size = (2200, 900))
    Label(f[0, 1:n_cols], main_title, fontsize = 24, font = :bold)

    col_colors = distinguishable_colors(n_cols, [colorant"white", colorant"black"], dropseed=true, lchoices=30:70)
    axes_grid = Matrix{Any}(nothing, n_rows, n_cols)

    for (r, volt) in enumerate(sorted_volts)
        Label(f[r, n_cols + 1], volt, rotation = -pi/2, font=:bold, fontsize=18)

        for (c, size_str) in enumerate(sorted_sizes)
            obj = get(cable_lookup, (volt, size_str), nothing)
            if isnothing(obj) continue end

            ax = Axis(f[r, c],
                title = r == 1 ? "$size_str mm²" : "",
                titlesize = 14,
                xticklabelsvisible = (r == n_rows),
                yticklabelsvisible = false,
                xgridvisible = false, ygridvisible = false
            )
            axes_grid[r, c] = ax

            _draw_ridgeline_on_axis!(
                ax, obj, quantity, indices;
                sequence = sequence,
                scale_factor = scale_factor,
                base_color = col_colors[c],
                overlap = 0.5, alpha = 0.85
            )

            if r == n_rows; ax.xlabel = y_lbl; end
        end
    end

    # Link X-axes within columns
    for c in 1:n_cols
        col_axes = filter(!isnothing, axes_grid[:, c])
        length(col_axes) > 1 && linkxaxes!(col_axes...)
    end

    return f
end

"""
    plot_frequency_fanchart(obj, quantity; ...)

Plots the trend of a parameter over frequency with shaded confidence intervals (90% and 50% bands).
"""
function plot_frequency_fanchart(
    obj::LineParametersMCSummary,
    quantity::Symbol;
    sequence::Symbol = :positive,
    fig_size::Tuple{Int, Int} = (800, 500),
    colormap::Symbol = :Blues,
    backend = nothing
)
    (_, seq_title) = _get_seq_info(sequence)
    
    # 1. Compute Data using Helper
    freqs, qs = _compute_quantile_data(obj, quantity, sequence)
    
    # 2. Setup Figure
    unit_label = (quantity == :L) ? "mH/km" : "Ω/km"
    full_title = "Frequency Dependence of $(seq_title) $(quantity)"
    
    backend_ctx = _make_window(BackendHandler, backend; title = full_title, icons = Engine._ICON_FN, icons_font = ICON_TTF)
    fig = Makie.Figure(size = fig_size)
    
    ax = Makie.Axis(fig[1, 1], 
        xlabel = "Frequency [Hz]", 
        ylabel = "$(quantity) [$(unit_label)]",
        title = full_title,
        xscale = log10
    )
    
    # 3. Draw Bands
    c_scheme = cgrad(colormap)
    # 90% CI (5% to 95%)
    band!(ax, freqs, qs[1,:], qs[5,:], color = (c_scheme[0.3], 0.5), label = "90% CI")
    # 50% CI (25% to 75%)
    band!(ax, freqs, qs[2,:], qs[4,:], color = (c_scheme[0.6], 0.6), label = "50% CI")
    # Median
    lines!(ax, freqs, qs[3,:], color = c_scheme[0.9], linewidth = 3, label = "Median")

    axislegend(ax, position = :rt)
    
    if backend_ctx.interactive && !isnothing(backend_ctx.window)
        display(backend_ctx.window, fig)
    else
        BackendHandler.renderfig(fig)
    end
    
    return fig
end

"""
    plot_multi_voltage_fanchart_grid(cable_lookup, sorted_volts, sorted_sizes, quantity; ...)

Creates a giant grid of fancharts (Row=Voltage, Col=Size).
"""
function plot_multi_voltage_fanchart_grid(
    cable_lookup::AbstractDict,
    sorted_volts::Vector{String},
    sorted_sizes::Vector{String}, 
    quantity::Symbol;
    sequence::Symbol = :positive,
    main_title::String = "",
)
    n_rows, n_cols = length(sorted_volts), length(sorted_sizes)
    f = Figure(size = (2000, 900))
    Label(f[0, 1:n_cols], main_title, fontsize = 24, font = :bold)

    # Unit label determination
    unit_label = (quantity == :L) ? "mH/km" : "Ω/km"

    # Distinct colors for columns so it's easy to track size vertically
    col_colors = distinguishable_colors(n_cols, [colorant"white", colorant"black"], dropseed=true, lchoices=30:70)
    
    axes_list = Axis[]

    for (r, volt) in enumerate(sorted_volts)
        # Row Label (Voltage)
        Label(f[r, n_cols + 1], volt, rotation = -pi/2, font=:bold, fontsize=18)

        for (c, size_str) in enumerate(sorted_sizes)
            # Fetch Object
            obj = get(cable_lookup, (volt, size_str), nothing)
            if isnothing(obj) continue end

            # Create Axis
            ax = Axis(f[r, c],
                title = r == 1 ? "$size_str mm²" : "",
                titlesize = 14,
                xscale = log10,
                xminorticksvisible = true, yminorticksvisible = true,
                xminorgridvisible = true, yminorgridvisible = true,
                # Labels only on edges
                xlabel = (r == n_rows ? "Frequency (Hz)" : ""),
                ylabel = (c == 1 ? "$(quantity) [$(unit_label)]" : "")
            )
            push!(axes_list, ax)

            # Compute Data
            freqs, qs = _compute_quantile_data(obj, quantity, sequence)

            # Draw Fanchart
            # We derive a specific color palette based on the column color
            base_col = col_colors[c]
            c_light  = weighted_color_mean(0.3, base_col, colorant"white")
            c_med    = weighted_color_mean(0.6, base_col, colorant"white")
            
            # 90% CI
            band!(ax, freqs, qs[1,:], qs[5,:], color = (c_light, 0.4))
            # 50% CI
            band!(ax, freqs, qs[2,:], qs[4,:], color = (c_med, 0.5))
            # Median
            lines!(ax, freqs, qs[3,:], color = base_col, linewidth = 2)

            # Hide tick labels if inside grid
            if r < n_rows; hidexdecorations!(ax, grid=false, minorgrid=false); end
            if c > 1; hideydecorations!(ax, grid=false, minorgrid=false); end
        end
    end

    # Link all axes for consistent zooming
    if !isempty(axes_list); linkaxes!(axes_list...); end

    return f
end

"""
    plot_spectral_comparison(cables, quantity; ...)

Plots the median value of `quantity` vs Frequency for multiple cables on one plot.
Useful for comparing different sizes (e.g. 150mm vs 240mm) directly.
"""
function plot_spectral_comparison(
    items::Vector{CablePlotItem},
    quantity::Symbol;
    sequence::Symbol = :positive,
    title_str::String = "Spectral Comparison",
    fig_size::Tuple{Int, Int} = (900, 600),
    backend = nothing
)
    # Sort by conductor size
    sort!(items, by = item -> item.conductor)
    
    unit_label = (quantity == :L) ? "mH/km" : "Ω/km"
    
    backend_ctx = _make_window(BackendHandler, backend; title = title_str, icons = _ICON_FN, icons_font = ICON_TTF)
    f = Figure(size = fig_size)
    
    ax = Axis(f[1, 1], 
        title = title_str, 
        xlabel = "Frequency (Hz)", 
        ylabel = "$(quantity) [$(unit_label)]",
        xscale = log10,
        xminorticksvisible = true, yminorticksvisible = true,
        xminorgridvisible = true, yminorgridvisible = true
    )
    
    colors_list = distinguishable_colors(length(items), [colorant"white", colorant"black"], dropseed=true, lchoices=30:70)
    
    for (i, item) in enumerate(items)
        # Pass the inner object (.obj) to the data extractor
        freqs, qs = _compute_quantile_data(item.obj, quantity, sequence)
        median_curve = qs[3, :]
        
        lines!(ax, freqs, median_curve, 
            label = "$(item.conductor)/$(item.screen) mm²", # <--- Access metadata from wrapper
            color = colors_list[i], 
            linewidth = 3
        )
    end
    
    axislegend(ax, position = :rt, framevisible = false)
    
    if backend_ctx.interactive && !isnothing(backend_ctx.window)
        display(backend_ctx.window, f)
    else
        BackendHandler.renderfig(f)
    end
    
    return f
end

end # module