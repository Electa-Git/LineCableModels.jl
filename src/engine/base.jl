
Base.eltype(::LineParametersProblem{T}) where {T} = T
Base.eltype(::Type{LineParametersProblem{T}}) where {T} = T

Base.eltype(::LineParameters{T}) where {T} = T
Base.eltype(::Type{LineParameters{T}}) where {T} = T

Base.eltype(::EMTWorkspace{T}) where {T} = T
Base.eltype(::Type{EMTWorkspace{T}}) where {T} = T

abstract type UnitLen end
struct PerMeter <: UnitLen end
struct PerKilometer <: UnitLen end
_len_scale(::PerMeter) = 1.0
_len_scale(::PerKilometer) = 1_000.0
_len_label(::PerMeter) = "m"
_len_label(::PerKilometer) = "km"

abstract type DisplayMode end
struct AsZY <: DisplayMode end
struct AsRLCG <: DisplayMode end

using Printf
using Measurements: value, uncertainty

"""
ResultsView: pretty, non-mutating renderer with zero-clipping and units.

- `mode = AsZY()` prints Z [Ω/len], Y [S/len].
- `mode = AsRLCG()` prints R [Ω/len], L [mH/len], G [S/len], C [µF/len];
  frequency is taken from the `LineParameters.f` vector of the view.
- `tol` clips tiny magnitudes to 0.0 in display only (value & uncertainty).
"""
struct ResultsView{LP <: LineParameters, U <: UnitLen, M <: DisplayMode}
	lp::LP
	unit::U
	mode::M
	tol::Float64
end

# Builders
resultsview(
	lp::LineParameters;
	per::Symbol = :km,
	mode::Symbol = :ZY,
	tol::Real = sqrt(eps(Float64)),
) = ResultsView(
	lp,
	per === :km ? PerKilometer() : PerMeter(),
	mode === :ZY ? AsZY() : AsRLCG(),
	float(tol),
)

# --- Scalar formatting with zero-clipping -------------------------------------

# Zero-clip helpers (display only)
_clip(x::Real, tol) = (abs(x) < tol ? 0.0 : x)

_format_real(io, x::Real, tol) = @printf(io, "%.6g", _clip(x, tol))

_format_meas(io, m, tol) = begin
	v = _clip(value(m), tol)
	u = _clip(uncertainty(m), tol)
	@printf(io, "%.6g±%.6g", v, u)
end

_format_complex(io, z, tol) = begin
	# z may be Complex{<:Real} or Complex{<:Measurement}
	print(io, "")
	if z.re isa Real
		_format_real(io, real(z), tol)
	else
		_format_meas(io, real(z), tol)
	end
	print(io, "+")
	if z.im isa Real
		_format_real(io, imag(z), tol)
	else
		_format_meas(io, imag(z), tol)
	end
	print(io, "im")
end

_format_any(io, x, tol) =
	x isa Complex ? _format_complex(io, x, tol) :
	x isa Measurements.Measurement ? _format_meas(io, x, tol) :
	_format_real(io, x, tol)

# String versions (for aligned, copy-pastable matrix literals)
_repr_real(x::Real, tol) = @sprintf("%.6g", _clip(x, tol))
_repr_meas(m, tol) = begin
	v = _clip(value(m), tol)
	u = _clip(uncertainty(m), tol)
	@sprintf("%.6g±%.6g", v, u)
end
function _repr_complex(z, tol)
	if z.re isa Real
		rs = _repr_real(real(z), tol)
	else
		rs = _repr_meas(real(z), tol)
	end
	if z.im isa Real
		is = _repr_real(imag(z), tol)
	else
		is = _repr_meas(imag(z), tol)
	end
	return string(rs, "+", is, "im")
end
_repr_any(x, tol) =
	x isa Complex ? _repr_complex(x, tol) :
	x isa Measurements.Measurement ? _repr_meas(x, tol) : _repr_real(x, tol)

# Detect if any element would be clipped by tolerance after mapping
_would_clip(x::Real, tol) = (x != 0 && abs(x) < tol)
_would_clip_meas(m, tol) = _would_clip(value(m), tol) || _would_clip(uncertainty(m), tol)
function _would_clip_complex(z, tol)
	(z.re isa Real ? _would_clip(real(z), tol) : _would_clip_meas(real(z), tol)) ||
		(z.im isa Real ? _would_clip(imag(z), tol) : _would_clip_meas(imag(z), tol))
end
_would_clip_any(x, tol) =
	x isa Complex ? _would_clip_complex(x, tol) :
	x isa Measurements.Measurement ? _would_clip_meas(x, tol) : _would_clip(x, tol)

function _any_clipped(A::AbstractMatrix; tol::Float64, map::Function = identity)
	n1, n2 = size(A, 1), size(A, 2)
	@inbounds for i in 1:n1, j in 1:n2
		x = map(A[i, j])
		_would_clip_any(x, tol) && return true
	end
	return false
end



# --- Show methods --------------------------------------------------------------

function _show_matrix(io::IO, A::AbstractArray; tol::Float64, map::Function = identity)
	n1, n2 = size(A, 1), size(A, 2)
	for i in 1:n1
		for j in 1:n2
			j > 1 && print(io, "  ")
			_format_any(io, map(A[i, j]), tol)
		end
		i < n1 && print(io, '\n')
	end
end

# Copy-pastable Julia matrix literal with column alignment
function _show_matrix_literal(
	io::IO,
	A::AbstractMatrix;
	tol::Float64,
	map::Function = identity,
)
	n1, n2 = size(A, 1), size(A, 2)
	# Build string table
	S = [_repr_any(map(A[i, j]), tol) for i in 1:n1, j in 1:n2]
	# Column widths
	widths = [maximum(length(S[i, j]) for i in 1:n1) for j in 1:n2]
	# Print rows
	for i in 1:n1
		print(io, i == 1 ? "[" : " ")
		for j in 1:n2
			s = S[i, j]
			pad = widths[j] - length(s)
			# right align
			print(io, " "^pad, s)
			if j < n2
				print(io, " ")
			end
		end
		if i < n1
			print(io, ";\n")
		else
			print(io, "]")
		end
	end
end

function Base.show(io::IO, ::MIME"text/plain", rv::ResultsView)
	lp = rv.lp
	unit = rv.unit
	tol = rv.tol
	scale = _len_scale(unit)
	ulabel = _len_label(unit)
	_, _, nf = size(lp.Z)

	# Determine if any value would be clipped across displayed content
	any_clipped = false
	if rv.mode isa AsZY
		@inbounds for k in 1:nf
			Zk = lp.Z.values[:, :, k]
			Yk = lp.Y.values[:, :, k]
			any_clipped |= _any_clipped(Zk; tol = tol, map = x -> scale * x)
			any_clipped && break
			any_clipped |= _any_clipped(Yk; tol = tol, map = x -> scale * x)
			any_clipped && break
		end
	else
		@inbounds for k in 1:nf
			Zk = lp.Z.values[:, :, k]
			Yk = lp.Y.values[:, :, k]
			fk = lp.f[k]
			ω = 2 * pi * float(fk)
			any_clipped |=
				_any_clipped(Zk; tol = tol, map = x -> scale * real(x)) ||
				_any_clipped(Zk; tol = tol, map = x -> (scale * 1e3 / ω) * imag(x)) ||
				_any_clipped(Yk; tol = tol, map = x -> scale * real(x)) ||
				_any_clipped(Yk; tol = tol, map = x -> (scale * 1e6 / ω) * imag(x))
			any_clipped && break
		end
	end

	# Styled header similar to DataFrame-like formatting
	n, _, _ = size(lp.Z)
	mode_label = rv.mode isa AsZY ? "ZY" : "RLCG"
	tol_str = @sprintf("%.1e", tol)
	header_plain =
		@sprintf("%dx%dx%d LineParameters |  mode = %s  |  units per %s  |  tol = %s%s",
			n, n, nf, mode_label, ulabel, tol_str, any_clipped ? " (!)" : "")
	printstyled(io, @sprintf("%dx%dx%d LineParameters", n, n, nf); bold = true)
	print(io, " |  mode = ")
	printstyled(io, mode_label; bold = true, color = :cyan)
	print(io, "  |  units per ", ulabel, "  |  tol = ", tol_str)
	if any_clipped
		print(io, " ")
		printstyled(io, "(!)"; bold = true, color = :yellow)
	end
	print(io, "\n")
	print(io, repeat("─", length(header_plain)))
	print(io, "\n\n")

	@views for k in 1:nf
		# Slice header with frequency of the slice
		fk = lp.f[k]
		print(io, "\n[:, :, ", k, "]  @ f=")
		print(io, @sprintf("%.6g", float(fk)))
		print(io, " Hz\n")
		Zk = lp.Z.values[:, :, k]
		Yk = lp.Y.values[:, :, k]

		if rv.mode isa AsZY
			println(io, "Z [Ω/", ulabel, "] =")
			_show_matrix_literal(io, Zk; tol = tol, map = x -> scale * x)

			print(io, "\n\nY [S/", ulabel, "] =\n")
			_show_matrix_literal(io, Yk; tol = tol, map = x -> scale * x)
		else
			# derive ω from frequency vector for this slice
			ω = 2 * pi * float(fk)

			println(io, "R [Ω/", ulabel, "] =")
			_show_matrix_literal(io, Zk; tol = tol, map = x -> scale * real(x))

			print(io, "\n\nL [mH/", ulabel, "] =\n")
			_show_matrix_literal(io, Zk; tol = tol, map = x -> (scale * 1e3 / ω) * imag(x))

			print(io, "\n\nG [S/", ulabel, "] =\n")
			_show_matrix_literal(io, Yk; tol = tol, map = x -> scale * real(x))

			print(io, "\n\nC [µF/", ulabel, "] =\n")
			_show_matrix_literal(io, Yk; tol = tol, map = x -> (scale * 1e6 / ω) * imag(x))
		end

		k < nf && print(io, "\n", "---"^10, "\n")
	end
end

function Base.show(io::IO, ::MIME"text/plain", Z::SeriesImpedance)
	n, _, nf = size(Z.values)
	header_plain = @sprintf("%dx%dx%d SeriesImpedance [Ω/m]", n, n, nf)
	printstyled(io, header_plain; bold = true)
	print(io, "\n")
	print(io, repeat("─", length(header_plain)))
	print(io, "\n")
	@views _show_matrix(io, Z.values[:, :, 1]; tol = sqrt(eps(Float64)))
	size(Z, 3) > 1 && print(
		io,
		"\n… (",
		size(Z, 3) - 1,
		" more slice",
		size(Z, 3) - 1 == 1 ? "" : "s",
		")",
	)
end

function Base.show(io::IO, ::MIME"text/plain", Y::ShuntAdmittance)
	n, _, nf = size(Y.values)
	header_plain = @sprintf("%dx%dx%d ShuntAdmittance [S/m]", n, n, nf)
	printstyled(io, header_plain; bold = true)
	print(io, "\n")
	print(io, repeat("─", length(header_plain)))
	print(io, "\n")
	@views _show_matrix(io, Y.values[:, :, 1]; tol = sqrt(eps(Float64)))
	size(Y, 3) > 1 && print(
		io,
		"\n… (",
		size(Y, 3) - 1,
		" more slice",
		size(Y, 3) - 1 == 1 ? "" : "s",
		")",
	)
end



# ---- SeriesImpedance array-ish interface ----
Base.size(Z::SeriesImpedance) = size(Z.values)
Base.size(Z::SeriesImpedance, d::Int) = size(Z.values, d)
Base.axes(Z::SeriesImpedance) = axes(Z.values)
Base.ndims(::Type{SeriesImpedance{T}}) where {T} = 3
Base.eltype(::Type{SeriesImpedance{T}}) where {T} = T
Base.getindex(Z::SeriesImpedance, I...) = @inbounds Z.values[I...]

# ---- ShuntAdmittance array-ish interface ----
Base.size(Y::ShuntAdmittance) = size(Y.values)
Base.size(Y::ShuntAdmittance, d::Int) = size(Y.values, d)
Base.axes(Y::ShuntAdmittance) = axes(Y.values)
Base.ndims(::Type{ShuntAdmittance{T}}) where {T} = 3
Base.eltype(::Type{ShuntAdmittance{T}}) where {T} = T
Base.getindex(Y::ShuntAdmittance, I...) = @inbounds Y.values[I...]

# --- Frequency-slice sugar ----------------------------------------------------
@inline Base.getindex(lp::LineParameters, k::Integer) = LineParameters(
	SeriesImpedance(@view lp.Z.values[:, :, k:k]),
	ShuntAdmittance(@view lp.Y.values[:, :, k:k]),
	lp.f[k:k],
)

# --- One-argument k, derive ω from freq (or accept ω directly) ---------------
function per_km(lp::LineParameters, k::Integer = 1;
	mode::Symbol = :ZY,
	tol::Real = sqrt(eps(Float64)))
	lpk = lp[k]
	return resultsview(lpk; per = :km, mode = mode, tol = tol)
end

function per_m(lp::LineParameters, k::Integer = 1;
	mode::Symbol = :ZY,
	tol::Real = sqrt(eps(Float64)))
	lpk = lp[k]
	return resultsview(lpk; per = :m, mode = mode, tol = tol)
end

# Helper: detect uncertainties in element type
_has_uncertainty_type(::Type{Complex{S}}) where {S} = S <: Measurement
_has_uncertainty_type(::Type) = false

# Terse summary (used inside collections)
function Base.show(io::IO, lp::LineParameters)
	n, _, nf = size(lp.Z)
	T = eltype(lp.Z)
	print(io, "LineParameters{$(T)} ", n, "×", n, "×", nf, "  [Z:Ω/m, Y:S/m]")
	_has_uncertainty_type(T) && print(io, " (±)")
end

function Base.show(io::IO, ::MIME"text/plain", lp::LineParameters)
	n, _, nf = size(lp.Z)
	T = eltype(lp.Z)
	tol = sqrt(eps(Float64))
	scale = 1_000.0   # per km preview
	ulabel = "km"

	# Styled header similar to ResultsView
	header_plain = string(
		n, "x", n, "x", nf, " LineParameters  |  eltype = ", T,
		_has_uncertainty_type(T) ? "  |  uncertainties: yes" : "",
	)

	printstyled(io, string(n, "x", n, "x", nf, " LineParameters"); bold = true)
	print(io, "  |  eltype = ", T)
	_has_uncertainty_type(T) && print(io, "  |  uncertainties: yes")
	print(io, "\n")
	print(io, repeat("─", length(header_plain)))
	print(io, "\n\n")

	# Preview: slice 1, per km, Z then Y
	@views begin
		Z1 = view(lp.Z.values,:,:,1)
		Y1 = view(lp.Y.values,:,:,1)

		println(io, "Preview (slice 1/", nf, ")  per ", ulabel)
		println(io, "Z [Ω/", ulabel, "] =")
		_show_matrix(io, Z1; tol = tol, map = x -> scale * x)

		print(io, "\n\nY [S/", ulabel, "] =\n")
		_show_matrix(io, Y1; tol = tol, map = x -> scale * x)
	end

	if nf > 1
		print(io, "\n\n… (", nf - 1, " more frequency slice", nf - 1 == 1 ? "" : "s", ")")
	end
end

