# -------------------------
# UI Logic Specs (Recipes)
# -------------------------

struct UIContainerSpec
	name::Symbol
	parent::Union{Nothing, Symbol}
	at::Tuple # Relaxed from strict Union types to generic Tuple
	layout::NamedTuple
end

UIContainerSpec(name, parent, at; layout = (;)) =
	UIContainerSpec(name, parent, at, layout)

struct UISlotSpec
	name::Symbol
	parent::Symbol
	at::Tuple # Relaxed from strict Union types
	layout::NamedTuple # Grid properties (e.g. alignmode, height)
	attrs::NamedTuple  # Content properties (e.g. Axis background, Legend align)
end

# Robust helper constructor
UISlotSpec(name, parent, at; layout = (;), attrs = (;)) =
	UISlotSpec(name, parent, at, layout, attrs)

struct UILayoutSpec
	name::Symbol
	containers::Vector{UIContainerSpec}
	slots::Vector{UISlotSpec}
	rowsizes::Dict{Symbol, Vector{Any}}
	colsizes::Dict{Symbol, Vector{Any}}
end

abstract type UIWidgetSpec end

struct UIButtonSpec <: UIWidgetSpec
	label::String
	icon::Union{Nothing, String}
	action::Function # (ctx, uifig, button) -> nothing
	attrs::NamedTuple # Passed to Makie.Button
end

# Constructor
UIButtonSpec(label, icon, action; attrs = (;)) =
	UIButtonSpec(label, icon, action, attrs)

struct UIToggleSpec <: UIWidgetSpec
	label::String
	active::Bool
	action::Function # (ctx, uifig, active::Bool) -> nothing
	attrs::NamedTuple # Passed to Makie.Toggle
end

# Constructor
UIToggleSpec(label, active, action; attrs = (;)) =
	UIToggleSpec(label, active, action, attrs)

# -------------------------
# UI Instances (Objects)
# -------------------------

mutable struct UIContext
	backend::Symbol
	interactive::Bool
	use_latex_fonts::Bool
	window::Union{Nothing, Any}
	screen::Union{Nothing, Any}
	status::Union{Nothing, Makie.Observable{String}}
	theme::Makie.Theme
end

struct UIFigure
	figure::Makie.Figure
	layoutspec::UILayoutSpec
	containers::Dict{Symbol, Makie.GridLayout}
	slots::Dict{Symbol, Any}
	cursor::Base.RefValue{Int}
	panelshape::Tuple{Int, Int}
end

struct UIPanel
	view::ViewSpec
	axis::Any
	plots::Vector{Any}
end

struct UIPlot
	spec::DataType
	ctx::UIContext
	page::PageSpec
	uifig::UIFigure
	panels::Vector{UIPanel}
	widgets::Dict{Symbol, Any}
end