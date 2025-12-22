
"""
	build_axes(::Type{S}, nt::NamedTuple) where {S<:AbstractPlotSpec}

Build axes for spec `S` using quantity tags stored in `nt` as
fields `x_quantity`, `y_quantity`, `z_quantity` (when applicable).

Returns:
	(xaxis = PBAxis or nothing,
	 yaxis = PBAxis or nothing,
	 zaxis = PBAxis or nothing)
"""
function build_axes(::Type{S}, nt::NamedTuple) where {S <: AbstractPlotSpec}
	dims = geom_axes(S)

	xaxis = nothing
	yaxis = nothing
	zaxis = nothing

	for dim in dims
		qfield = Symbol(dim, :_quantity)   # :x_quantity, :y_quantity, :z_quantity
		q = getproperty(nt, qfield)        # expected to be a QuantityTag

		u  = axis_unit(S, q, dim)
		ax = build_axis(S, dim, q, u)

		if dim === :x
			xaxis = ax
		elseif dim === :y
			yaxis = ax
		elseif dim === :z
			zaxis = ax
		else
			Base.error("Unsupported axis dim $(dim) in geom_axes for $(S)")
		end
	end

	return (xaxis = xaxis, yaxis = yaxis, zaxis = zaxis)
end

"""
Build the final axis label from quantity + units.

Default pattern: "Label [symbol]" if units are non-empty.
Uses `get_label(q::QuantityTag)` and `get_label(u::Units)` from UnitHandler.
"""
function axis_label(q::QuantityTag, u::Units)
	base = get_label(q)        # human-readable quantity name
	usym = get_label(u)        # unit symbol string, e.g. "Ω/km"
	return isempty(usym) ? base : string(base, " [", usym, "]")
end

"""
Build a PBAxis for spec `S`, axis dim `dim`, quantity `q`, units `u`.
"""
function build_axis(::Type{S},
	dim::Symbol,
	q::QuantityTag,
	u::Units) where {S <: AbstractPlotSpec}

	lab = axis_label(q, u)
	dims = enable_logscale(S)
	sc = dim in dims ? :log10 : :linear
	return PBAxis(dim, q, u, lab, sc)
end

# High-level builder: dim only (static semantics)
function build_axis(::Type{S}, ::Val{dim}) where {S <: AbstractPlotSpec, dim}
	q = axis_quantity(S, Val(dim))
	u = axis_unit(S, q, dim)
	return build_axis(S, dim, q, u)
end

# High-level builder: dim + semantic code (e.g. :R/:L/:C/:G)
function build_axis(
	::Type{S},
	::Val{dim},
	::Val{qty},
) where {S <: AbstractPlotSpec, dim, qty}
	q = axis_quantity(S, Val(dim), Val(qty))
	u = axis_unit(S, q, dim)
	return build_axis(S, dim, q, u)
end
