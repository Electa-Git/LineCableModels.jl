using Revise
using LineCableModels
using LineCableModels.UQ: LineParametersMCSummary
using JSON
using Statistics
using Dates
using JLD2
using CairoMakie 
using Colors

# Include the plotting functions
include("plot_runner.jl") 
using .plot_runner 

# --- CONFIG ---
input_folder = "System_Matrices"
output_folder = "Figures"
grid_output = joinpath(output_folder, "NA2XS2Y Series")
if !isdir(output_folder); mkdir(output_folder); end
if !isdir(grid_output); mkdir(grid_output); end

# --- LOAD DATA ---
grouped_data = Dict{String, Vector{plot_runner.CablePlotItem}}()

full_obj_lookup = Dict{Tuple{String, String}, Any}()

files = filter(f -> endswith(f, ".jld2"), readdir(input_folder))
filename_pattern = r"1x(\d+)-(\d+)\s+(\d+-\d+kV)"

@info "Processing $(length(files)) files..."

# Parse each file and group by voltage level
for filename in files
    m = match(filename_pattern, filename)
    if !isnothing(m)
        cond_size = parse(Int, m.captures[1])
        scr_size  = parse(Int, m.captures[2])
        voltage   = m.captures[3]
        
        d = load(joinpath(input_folder, filename))
        
        if haskey(d, "res_sys")
            res = d["res_sys"]
            
            # Label for lookup
            lbl = "$(cond_size)/$(scr_size)"
            full_obj_lookup[(voltage, lbl)] = res
            
            if !haskey(grouped_data, voltage)
                grouped_data[voltage] = plot_runner.CablePlotItem[]
            end
            
            item = plot_runner.CablePlotItem(res, cond_size, scr_size)
            push!(grouped_data[voltage], item)
        end
    end
end

# --- SPECTRAL COMPARISONS ---
for (voltage, cables) in grouped_data
    safe_v = replace(voltage, "/" => "-")
    @info "Plotting group: $voltage"

    # 1. R Positive
    f1 = plot_runner.plot_spectral_comparison(
        cables, :R; 
        sequence=:positive, 
        title_str="Pos. Seq. Resistance ($voltage)"
    )
    Makie.save(joinpath(output_folder, "R_pos_NA2XS2Y_$(safe_v)_comparison.svg"), f1)
    @info "Saved R_pos comparison for $voltage"
    
    # 2. L Positive
    f2 = plot_runner.plot_spectral_comparison(
        cables, :L; 
        sequence=:positive, 
        title_str="Pos. Seq. Inductance ($voltage)"
    )
    Makie.save(joinpath(output_folder, "L_pos_NA2XS2Y_$(safe_v)_comparison.svg"), f2)
    @info "Saved L_pos comparison for $voltage"
    
    # 3. R Zero
    f3 = plot_runner.plot_spectral_comparison(
        cables, :R; 
        sequence=:zero, 
        title_str="Zero Seq. Resistance ($voltage)"
    )
    Makie.save(joinpath(output_folder, "R_zero_NA2XS2Y_$(safe_v)_comparison.svg"), f3)
    @info "Saved R_zero comparison for $voltage"
    
    # 4. L Zero
    f4 = plot_runner.plot_spectral_comparison(
        cables, :L; 
        sequence=:zero, 
        title_str="Zero Seq. Inductance ($voltage)"
    )
    Makie.save(joinpath(output_folder, "L_zero_NA2XS2Y_$(safe_v)_comparison.svg"), f4)
    @info "Saved L_zero comparison for $voltage"
end

# --- PREPARE GRIDS ---
# A. Sort Voltages (Rows)
sorted_volts = sort(collect(keys(grouped_data)))

# B. Sort Cross-Sections (Columns)
all_cables_flat = vcat(values(grouped_data)...)
unique_cables = unique(c -> (c.conductor, c.screen), all_cables_flat)
sort!(unique_cables, by = c -> c.conductor)
sorted_size_labels = ["$(c.conductor)/$(c.screen)" for c in unique_cables]

# --- GENERATE GRIDS ---
@info "Generating NA2XS2Y Series Grids..."
selected_indices = [1, 3, 7, 13, 20] # Select some frequencies for ridgeline plots

# === 1. R Positive ===
# Fanchart Grid
f5 = plot_runner.plot_multi_voltage_fanchart_grid(
    full_obj_lookup, sorted_volts, sorted_size_labels, :R;
    sequence = :positive,
    main_title = "Pos. Seq. Resistance Comparison"
)
Makie.save(joinpath(grid_output, "R_pos_NA2XS2Y_series_fan.svg"), f5)
@info "Saved R_pos Fanchart"

# Ridgeline Grid
f_ridge_1 = plot_runner.plot_multi_voltage_ridgeline_grid(
    full_obj_lookup, sorted_volts, sorted_size_labels, :R, selected_indices;    
    sequence = :positive,
    scale_factor = 1.0,
    y_lbl = "R₁ (Ω/km)",
    main_title = "Evolution of Positive Sequence Resistance Distribution",
)
Makie.save(joinpath(grid_output, "R_pos_NA2XS2Y_series_ridge.svg"), f_ridge_1)
@info "Saved R_pos Ridgeline"

# === 2. L Positive ===
# Fanchart Grid (Auto-scales to mH)
f6 = plot_runner.plot_multi_voltage_fanchart_grid(
    full_obj_lookup, sorted_volts, sorted_size_labels, :L;
    sequence = :positive,
    main_title = "Pos. Seq. Inductance Comparison"
)
Makie.save(joinpath(grid_output, "L_pos_NA2XS2Y_series_fan.svg"), f6)
@info "Saved L_pos Fanchart"

# Ridgeline Grid (Must manually scale to mH for consistency)
f_ridge_2 = plot_runner.plot_multi_voltage_ridgeline_grid(
    full_obj_lookup, sorted_volts, sorted_size_labels, :L, selected_indices;    
    sequence = :positive,
    scale_factor = 1000.0, # H -> mH
    y_lbl = "L₁ (mH/km)",
    main_title = "Evolution of Positive Sequence Inductance Distribution",
)
Makie.save(joinpath(grid_output, "L_pos_NA2XS2Y_series_ridge.svg"), f_ridge_2)
@info "Saved L_pos Ridgeline"

# === 3. R Zero ===
# Fanchart Grid
f7 = plot_runner.plot_multi_voltage_fanchart_grid(
    full_obj_lookup, sorted_volts, sorted_size_labels, :R;
    sequence = :zero,
    main_title = "Zero Seq. Resistance Comparison"
)
Makie.save(joinpath(grid_output, "R_zero_NA2XS2Y_series_fan.svg"), f7)
@info "Saved R_zero Fanchart"

# Ridgeline Grid
f_ridge_3 = plot_runner.plot_multi_voltage_ridgeline_grid(
    full_obj_lookup, sorted_volts, sorted_size_labels, :R, selected_indices;    
    sequence = :zero,
    scale_factor = 1.0,
    y_lbl = "R₀ (Ω/km)",
    main_title = "Evolution of Zero Sequence Resistance Distribution",
)
Makie.save(joinpath(grid_output, "R_zero_NA2XS2Y_series_ridge.svg"), f_ridge_3)
@info "Saved R_zero Ridgeline"

# === 4. L Zero ===
# Fanchart Grid (Auto-scales to mH)
f8 = plot_runner.plot_multi_voltage_fanchart_grid(
    full_obj_lookup, sorted_volts, sorted_size_labels, :L;
    sequence = :zero,
    main_title = "Zero Seq. Inductance Comparison"
)
Makie.save(joinpath(grid_output, "L_zero_NA2XS2Y_series_fan.svg"), f8)
@info "Saved L_zero Fanchart"

# Ridgeline Grid (Must manually scale to mH)
f_ridge_4 = plot_runner.plot_multi_voltage_ridgeline_grid(
    full_obj_lookup, sorted_volts, sorted_size_labels, :L, selected_indices;    
    sequence = :zero,
    scale_factor = 1000.0, # H -> mH
    y_lbl = "L₀ (mH/km)",
    main_title = "Evolution of Zero Sequence Inductance Distribution",
)
Makie.save(joinpath(grid_output, "L_zero_NA2XS2Y_series_ridge.svg"), f_ridge_4)
@info "Saved L_zero Ridgeline"

# display(f5), display(f_ridge_1), display(f6), display(f_ridge_2), display(f7), display(f_ridge_3), display(f8), display(f_ridge_4)