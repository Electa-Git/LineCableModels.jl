mutable struct UIContext
    backend::Symbol
    interactive::Bool
    window::Any
    figure::Any
    shell::Any
    canvas::Any
    status::Observable{String}
    panels::Vector{Any}
    widgets::Dict{Symbol, Any}
    legend::Any
    colorbars::Vector{Any}
    responsive_legend::Any
    legend_slot_grid::Any
    observers::Vector{Any}
    plot_reference::Base.RefValue{Any}
end

struct UIPanel
    metadata::Any
    axis::Any
    plots::Vector{Any}
    groups::Dict{Symbol, Vector{Any}}
    group_labels::Dict{Symbol, String}
    group_order::Vector{Symbol}
end
