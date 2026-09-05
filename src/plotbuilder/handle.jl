"""
$(TYPEDEF)

Return the small handle assembled around a native Makie figure by a
LineCableModels plotting call.

`figure`, `title`, every element of `axes`, every value in `controls`, `legend`,
every value in `panel_legends`, and every element of `colorbars` are the live
Makie objects used for display. They remain owned by the caller: mutating them
or adding native Makie plots changes the displayed figure and the state saved
by [`export_svg`](@ref).

The remaining fields are private addon state and export defaults.

$(TYPEDFIELDS)
"""
mutable struct UIPlot{F, A, C, B}
    "Native Makie figure owned by the caller."
    figure::F
    "Figure-wide native Makie title label, or `nothing`."
    title::Any
    "Native Makie axes in display order."
    axes::A
    "Native Makie controls keyed by their public purpose."
    controls::C
    "Native Makie legend, or `nothing`."
    legend::Any
    "Panel-scoped native Makie legends keyed by logical panel identity."
    panel_legends::Dict{Any, Any}
    "Native Makie colorbars in display order."
    colorbars::B
    "Private addon registry consumed by the optional Makie extension."
    addon_state::Any
    "Default base filename used by [`export_svg`](@ref)."
    export_name::String
    "Default SVG theme."
    export_theme::Symbol
    "Whether toolbar exports should be opened after writing."
    open_export::Bool
end

function UIPlot(
        figure,
        axes;
        title = nothing,
        controls = Dict{Symbol, Any}(),
        legend = nothing,
        panel_legends = Dict{Any, Any}(),
        colorbars = (),
        addon_state = nothing,
        export_name::AbstractString = "linecablemodels_plot",
        export_theme::Symbol = :default,
        open_export::Bool = true
)
    export_theme in (:default, :publication) || throw(ArgumentError(
        "export_theme must be :default or :publication",
    ))
    isempty(strip(export_name)) && throw(ArgumentError("export_name cannot be empty"))
    panel_legends isa AbstractDict || throw(ArgumentError(
        "panel_legends must be a dictionary keyed by logical panel identity",
    ))
    return UIPlot(
        figure,
        title,
        axes,
        controls,
        legend,
        Dict{Any, Any}(panel_legends),
        colorbars,
        addon_state,
        String(export_name),
        export_theme,
        open_export
    )
end

Base.summary(io::IO, ::UIPlot) = print(io, "Native Makie plot")
