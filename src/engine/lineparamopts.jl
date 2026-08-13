Base.@kwdef struct LineParamOptions
	"Skip user confirmation for overwriting results"
	force_overwrite::Bool = false
	"Reduce bundle conductors to equivalent single conductor"
	reduce_bundle::Bool = true
	"Eliminate grounded conductors from the system (Kron reduction)"
	kron_reduction::Bool = true
	"Enforce ideal transposition/snaking"
	ideal_transposition::Bool = true
	"Temperature correction"
	temperature_correction::Bool = true
	"Store primitive matrices"
	store_primitive_matrices::Bool = true
	"Verbosity level"
	verbosity::Int = 0
	"Log file path"
	logfile::Union{String, Nothing} = nothing
end


# --- Helpers to turn anything into a NamedTuple ----------------------------

_to_nt(nt::NamedTuple) = nt
_to_nt(p::Base.Pairs) = (; p...)
_to_nt(d::AbstractDict) = (; d...)
_to_nt(::Nothing) = (;)

# --- Generic key splitter + builder ---------------------------------------

const _COMMON_KEYS = Set(fieldnames(LineParamOptions))

_select_keys(nt::NamedTuple, allowed::Set{Symbol}) =
	(; (k => v for (k, v) in pairs(nt) if k in allowed)...)

function build_options(::Type{O}, opts;
	strict::Bool = true,
) where {O <: AbstractFormulationOptions}

	nt = _to_nt(opts)

	own_allowed = Set(filter(!=(:common), fieldnames(O)))
	common_nt   = _select_keys(nt, _COMMON_KEYS)
	own_nt      = _select_keys(nt, own_allowed)

	unknown = setdiff(Set(keys(nt)), union(_COMMON_KEYS, own_allowed))
	if strict && !isempty(unknown)
		throw(ArgumentError("Unknown option keys for $(O): $(collect(unknown))"))
	end

	return O(; common = LineParamOptions(; common_nt...), own_nt...)
end

# Convenience overloads (accept already-built things)
build_options(::Type{O}, o::O; kwargs...) where {O <: AbstractFormulationOptions} = o
build_options(
	::Type{O},
	c::LineParamOptions;
	kwargs...,
) where {O <: AbstractFormulationOptions} = O(; common = c)

# save_path stays solver-specific (different sensible defaults).
Base.@kwdef struct EMTOptions <: AbstractFormulationOptions
	common::LineParamOptions = LineParamOptions()
	"Save path for output files"
	save_path::String = joinpath(".", "lineparams_output")
end

const _COMMON_SYMS = Tuple(fieldnames(LineParamOptions))
const _EMT_OWN = Tuple(s for s in fieldnames(EMTOptions) if s != :common)
@inline Base.hasproperty(::EMTOptions, s::Symbol) =
	(s in _EMT_OWN) || (s in _COMMON_SYMS) || s === :common

@inline function Base.getproperty(o::EMTOptions, s::Symbol)
	s === :common && return getfield(o, :common)
	(s in _EMT_OWN) && return getfield(o, s)          # EMT-specific
	(s in _COMMON_SYMS) && return getfield(o.common, s)   # forwarded common
	throw(ArgumentError("Unknown option $(s) for $(typeof(o))"))
end

Base.propertynames(::EMTOptions, ::Bool = false) = (_COMMON_SYMS..., _EMT_OWN..., :common)
Base.get(o::EMTOptions, s::Symbol, default) =
	hasproperty(o, s) ? getproperty(o, s) : default
asnamedtuple(o::EMTOptions) = (; (k=>getproperty(o, k) for k in propertynames(o))...)





