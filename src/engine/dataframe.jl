import DataFrames: DataFrame, metadata!

const _LP_FREQ_COL = :frequency

const _SERIES_DUMMY = SeriesImpedance(zeros(ComplexF64, 1, 1, 1))
const _SHUNT_DUMMY = ShuntAdmittance(zeros(ComplexF64, 1, 1, 1))

_freq_units_label(unit::Symbol) = unit_text(unit, "Hz")

_length_unit(per::Symbol) = per

function _column_name(meta::ComponentMetadata)
	component = meta.component
	if component in (:resistance, :inductance, :conductance, :capacitance)
		return Symbol(meta.symbol)
	else
		return Symbol(component)
	end
end

function _normalize_quantity_units(units)
	return normalize_quantity_units(units)
end

function _frequency_vector(obj, freqs)
	if freqs === nothing
		return float.(collect(axes(obj, 3)))
	else
		f = collect(freqs)
		length(f) == size(obj, 3) ||
			Base.error("Frequency vector length does not match object samples")
		return float.(f)
	end
end

function _frequency_vector(slice::AbstractVector, freqs::AbstractVector)
	f = collect(freqs)
	length(f) == length(slice) ||
		Base.error("Frequency vector length must match slice length")
	return float.(f)
end

function _build_dataframe(
	slice,
	freq_raw::Vector{<:Real},
	comps::Vector{ComponentMetadata},
	units::Dict{Symbol, Symbol},
	length_unit::Symbol,
	freq_unit::Symbol,
	tol::Real,
)
	freq_scale = frequency_scale(freq_unit)
	freq_values = freq_raw .* freq_scale
	unit_map = Dict{Symbol, String}(
		_LP_FREQ_COL => _freq_units_label(freq_unit),
	)
	df = DataFrame(_LP_FREQ_COL => freq_values)
	for meta in comps
		q_prefix = resolve_quantity_prefix(meta.quantity, units)
		scale = quantity_scale(q_prefix)
		l_scale = meta.unit.per_length ? length_scale(length_unit) : 1.0
		raw_vals = component_values(meta.component, slice, freq_raw)
		col_data = map(raw_vals) do x
			_clip_field(x * (scale * l_scale), tol)
		end
		col_name = _column_name(meta)
		df[!, col_name] = col_data
		unit_map[col_name] =
			composite_unit(q_prefix, meta.unit.symbol, meta.unit.per_length, length_unit)
	end
	metadata!(df, "units", unit_map, style = :note)
	return df
end

function _matrix_dataframes(
	obj,
	freq_raw::Vector{<:Real},
	comps::Vector{ComponentMetadata},
	units::Dict{Symbol, Symbol},
	length_unit::Symbol,
	freq_unit::Symbol,
	tol::Real,
)
	nx, ny, _ = size(obj.values)
	result = Matrix{DataFrame}(undef, nx, ny)
	for i in 1:nx, j in 1:ny
		slice = @view obj.values[i, j, :]
		result[i, j] = _build_dataframe(
			slice,
			freq_raw,
			comps,
			units,
			length_unit,
			freq_unit,
			tol,
		)
	end
	return result
end

function _slice_dataframe(
	slice::AbstractVector,
	kind::Symbol,
	freq_raw::Vector{<:Real},
	mode::Symbol,
	coord::Symbol,
	units::Dict{Symbol, Symbol},
	length_unit::Symbol,
	freq_unit::Symbol,
	tol::Real,
)
	resolved_kind = _resolve_kind(slice, kind, tol)
	comps =
		resolved_kind == :series_impedance ?
		components_for(_SERIES_DUMMY, mode, coord) :
		components_for(_SHUNT_DUMMY, mode, coord)
	return _build_dataframe(
		slice,
		freq_raw,
		comps,
		units,
		length_unit,
		freq_unit,
		tol,
	)
end

"""
	DataFrame(Z::SeriesImpedance; freqs=nothing, mode=:RLCG, coord=:cart,
			  freq_unit=:base, length_unit=:kilo, quantity_units=nothing,
			  tol=sqrt(eps(Float64)))

Convert the entries of a `SeriesImpedance` object into per-element `DataFrame`s
indexed by frequency. Returns an `n×n` matrix of `DataFrame`s whose rows
correspond to conductor indices.

- `freqs`: explicit frequency vector in Hz. Defaults to `1:length(freq axis)`.
- `mode`: `:RLCG` (default) or `:ZY`. For `:ZY`, `coord` may be `:cart` or `:polar`.
- `length_unit`: metric prefix for per-length units (e.g. `:kilo` ⇒ per km).
- `quantity_units`: optional overrides for the quantity metric prefixes used in each column.
- `tol`: absolute tolerance used to zero-out tiny numerical noise.
"""
function DataFrame(
	Z::SeriesImpedance;
	freqs = nothing,
	mode::Symbol = :RLCG,
	coord::Symbol = :cart,
	freq_unit::Symbol = :base,
	length_unit::Symbol = :kilo,
	quantity_units = nothing,
	tol::Real = sqrt(eps(Float64)),
)
	freq_raw = _frequency_vector(Z, freqs)
	units = _normalize_quantity_units(quantity_units)
	comps = components_for(Z, mode, coord)
	return _matrix_dataframes(
		Z,
		freq_raw,
		comps,
		units,
		length_unit,
		freq_unit,
		float(tol),
	)
end

"""
	DataFrame(Y::ShuntAdmittance; freqs=nothing, mode=:RLCG, coord=:cart,
			  freq_unit=:base, length_unit=:kilo, quantity_units=nothing,
			  tol=sqrt(eps(Float64)))

Convert the entries of a `ShuntAdmittance` object into per-element `DataFrame`s
indexed by frequency. Returns an `n×n` matrix of `DataFrame`s.

Keyword arguments mirror those of `DataFrame(::SeriesImpedance)`.
"""
function DataFrame(
	Y::ShuntAdmittance;
	freqs = nothing,
	mode::Symbol = :RLCG,
	coord::Symbol = :cart,
	freq_unit::Symbol = :base,
	length_unit::Symbol = :kilo,
	quantity_units = nothing,
	tol::Real = sqrt(eps(Float64)),
)
	freq_raw = _frequency_vector(Y, freqs)
	units = _normalize_quantity_units(quantity_units)
	comps = components_for(Y, mode, coord)
	return _matrix_dataframes(
		Y,
		freq_raw,
		comps,
		units,
		length_unit,
		freq_unit,
		float(tol),
	)
end

function _clip_field(x::Real, tol)
	isfinite(x) || return x
	return _clip(x, tol)
end

function _clip_field(m::Measurements.Measurement, tol)
	v = _clip(value(m), tol)
	u = _clip(uncertainty(m), tol)
	return Measurements.measurement(v, u)
end

_clip_field(x, _) = x

function _resolve_kind(slice, kind::Symbol, tol::Real)
	kind != :auto && return kind
	max_real = 0.0
	max_imag = 0.0
	for z in slice
		r = real(z)
		i = imag(z)
		val_r = _scalar_abs(r)
		val_i = _scalar_abs(i)
		isfinite(val_r) && val_r > max_real && (max_real = val_r)
		isfinite(val_i) && val_i > max_imag && (max_imag = val_i)
	end
	if max_real <= tol && max_imag > tol
		return :shunt_admittance
	else
		return :series_impedance
	end
end

_scalar_abs(x::Real) = abs(x)
_scalar_abs(m::Measurements.Measurement) = abs(value(m))

"""
	DataFrame(LP::LineParameters; mode=:RLCG, coord=:cart,
			  freq_unit=:base, length_unit=:kilo, quantity_units=nothing,
			  tol=sqrt(eps(Float64)))
Convert `LP.Z` and `LP.Y` to per-element, frequency-indexed `DataFrame`s
using `LP.f` as the authoritative frequency vector. Returns `(df_z, df_y)`,
each an `n×n` `Matrix{DataFrame}`.
"""
function DataFrame(
	LP::LineParameters;
	mode::Symbol = :RLCG,
	coord::Symbol = :cart,
	freq_unit::Symbol = :base,
	length_unit::Symbol = :kilo,
	quantity_units = nothing,
	tol::Real = sqrt(eps(Float64)),
)
	# --- validations: LP is the source of truth for frequency samples ----
	@assert eltype(LP.f) <: Real "LP.f must be real-valued frequencies."
	nzx, nzy, nfZ = size(LP.Z.values)
	nyx, nyy, nfY = size(LP.Y.values)
	nfZ == nfY ||
		Base.error("Z and Y have different number of frequency samples: $nfZ ≠ $nfY.")
	length(LP.f) == nfZ || Base.error(
		"Length of LP.f ($(length(LP.f))) does not match samples in Z/Y ($nfZ).",
	)

	# --- delegate with LP.f explicitly (no guessing, no manual input) ----
	df_z = DataFrame(
		LP.Z;
		freqs = LP.f,
		mode = mode,
		coord = coord,
		freq_unit = freq_unit,
		length_unit = length_unit,
		quantity_units = quantity_units,
		tol = tol,
	)

	df_y = DataFrame(
		LP.Y;
		freqs = LP.f,
		mode = mode,
		coord = coord,
		freq_unit = freq_unit,
		length_unit = length_unit,
		quantity_units = quantity_units,
		tol = tol,
	)

	return df_z, df_y
end
