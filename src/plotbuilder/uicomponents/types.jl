# -------------------------
# Declarative layout specs
# -------------------------

struct UIContainerSpec
	name::Symbol
	parent::Union{Nothing, Symbol}
	at::Tuple{Union{Int, UnitRange{Int}}, Union{Int, UnitRange{Int}}}
	shape::Tuple{Int, Int}
	props::NamedTuple

	function UIContainerSpec(
		name::Symbol,
		parent::Union{Nothing, Symbol},
		at::Tuple{Union{Int, UnitRange{Int}}, Union{Int, UnitRange{Int}}},
		shape::Tuple{Int, Int};
		props::NamedTuple = NamedTuple(),
	)
		return new(name, parent, at, shape, props)
	end
end

struct UISlotSpec
	name::Symbol
	parent::Symbol
	at::Tuple{Union{Int, UnitRange{Int}}, Union{Int, UnitRange{Int}}}
	props::NamedTuple

	function UISlotSpec(
		name::Symbol,
		parent::Symbol,
		at::Tuple{Union{Int, UnitRange{Int}}, Union{Int, UnitRange{Int}}};
		props::NamedTuple = NamedTuple(),
	)
		return new(name, parent, at, props)
	end
end

struct UILayoutSpec
	name::Symbol
	containers::Vector{UIContainerSpec}
	slots::Vector{UISlotSpec}
end

# -------------------------
# Materialized UI objects
# -------------------------

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
	meta::NamedTuple
end

mutable struct UIContext
	backend::Symbol
	interactive::Bool
	window::Union{Nothing, Any}
	screen::Union{Nothing, Any}
	status::Union{Nothing, Makie.Observable{String}}
	theme::Union{Nothing, Makie.Theme}
end

struct PlotAssembly
	spec::DataType
	ctx::UIContext
	pbfig::PageSpec
	uifig::UIFigure
	panels::Vector{UIPanel}
	widgets::Dict{Symbol, Any}
	meta::NamedTuple
end
