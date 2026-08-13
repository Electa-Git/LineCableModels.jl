Base.@kwdef struct FEMOptions <: AbstractFormulationOptions
	common::LineParamOptions = LineParamOptions()

	"Build mesh only and preview (no solving)"
	mesh_only::Bool = false
	"Force mesh regeneration even if file exists"
	force_remesh::Bool = false
	"Generate field visualization outputs"
	plot_field_maps::Bool = true
	"Archive temporary files after each frequency run"
	keep_run_files::Bool = false

	"Base path for output files"
	save_path::String = joinpath(".", "fem_output")
	"Path to GetDP executable"
	getdp_executable::Union{String, Nothing} = nothing
end

const _FEM_OWN = Tuple(s for s in fieldnames(FEMOptions) if s != :common)
@inline Base.hasproperty(::FEMOptions, s::Symbol) =
	(s in _FEM_OWN) || (s in _COMMON_SYMS) || s === :common

@inline function Base.getproperty(o::FEMOptions, s::Symbol)
	s === :common && return getfield(o, :common)
	(s in _FEM_OWN) && return getfield(o, s)          # FEM-specific
	(s in _COMMON_SYMS) && return getfield(o.common, s)   # forwarded common
	throw(ArgumentError("Unknown option $(s) for $(typeof(o))"))
end

Base.propertynames(::FEMOptions, ::Bool = false) = (_COMMON_SYMS..., _FEM_OWN..., :common)
Base.get(o::FEMOptions, s::Symbol, default) =
	hasproperty(o, s) ? getproperty(o, s) : default
asnamedtuple(o::FEMOptions) = (; (k=>getproperty(o, k) for k in propertynames(o))...)
# asnamedtuple(o::FEMOptions) = (; (k=>getproperty(o,k) for k in propertynames(o) if k != :common)...)
