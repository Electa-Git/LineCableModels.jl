Base.IteratorEltype(::Type{CableBuilderSpec}) = Base.HasEltype()
Base.eltype(::Type{CableBuilderSpec}) = DataModel.CableDesign
Base.IteratorSize(::Type{CableBuilderSpec}) = Base.SizeUnknown()

function Base.iterate(cbs::CableBuilderSpec)
	ch = iterate_designs(cbs)
	try
		d = take!(ch)
		return (d, ch)
	catch
		return nothing
	end
end

function Base.iterate(::CableBuilderSpec, ch::Channel)
	try
		d = take!(ch)
		return (d, ch)
	catch
		return nothing
	end
end


# how many choices are in a "range-like" thing
_choice_count(x) =
	x === nothing                   ? 1 :
	(x isa Tuple && length(x) == 2) ? _choice_count(x[1]) * _choice_count(x[2]) :
	(x isa AbstractVector)          ? length(x) :
	(x isa Tuple && length(x) == 3) ? last(x) : 1

# count choices for a MaterialSpec (rho/eps/mu/T/α product)
_choice_count(ms::MaterialSpec) = length(_make_range(ms))

# args: each entry can be scalar | vector | (lo,hi,n) | (value_spec, pct_spec)
_arg_choice_count(a) =
	(a isa Tuple && length(a) == 2) ? (_choice_count(a[1]) * _choice_count(a[2])) :
	_choice_count(a)

_args_choice_count(args::Tuple) =
	isempty(args) ? 1 : prod(_arg_choice_count(a) for a in args)

function cardinality(cbs::CableBuilderSpec)
	comp_names = unique(p.component for p in cbs.parts)
	by_comp = Dict{Symbol, Vector{PartSpec}}()
	for p in cbs.parts
		get!(by_comp, p.component, PartSpec[]) |> v -> push!(v, p)
	end

	total = 1
	for cname in comp_names
		ps = by_comp[cname]
		cond = [p for p in ps if p.part_type <: DataModel.AbstractConductorPart]
		insu = [p for p in ps if p.part_type <: DataModel.AbstractInsulatorPart]
		isempty(cond) && Base.error("component '$cname' has no conductors")
		isempty(insu) && Base.error("component '$cname' has no insulators")

		# first conductor axes
		p1c    = cond[1]
		c_dim  = _choice_count(p1c.dim[1]) * _choice_count(p1c.dim[2])
		c_args = _args_choice_count(p1c.args)
		c_mat  = _choice_count(p1c.material)

		# uncoupled extras from later conductors (couple when tuples compare equal)
		for pc in cond[2:end]
			pc_dim_same  = (pc.dim == p1c.dim)
			pc_args_same = (pc.args == p1c.args)
			pc_mat_same  = (pc.material == p1c.material)

			c_dim  *= pc_dim_same ? 1 : (_choice_count(pc.dim[1]) * _choice_count(pc.dim[2]))
			c_args *= pc_args_same ? 1 : _args_choice_count(pc.args)
			c_mat  *= pc_mat_same ? 1 : _choice_count(pc.material)
		end
		cond_factor = c_dim * c_args * c_mat

		# first insulator axes
		p1i    = insu[1]
		i_dim  = _choice_count(p1i.dim[1]) * _choice_count(p1i.dim[2])
		i_args = _args_choice_count(p1i.args)
		i_mat  = _choice_count(p1i.material)

		for pi in insu[2:end]
			pi_dim_same  = (pi.dim == p1i.dim)
			pi_args_same = (pi.args == p1i.args)
			pi_mat_same  = (pi.material == p1i.material)

			i_dim  *= pi_dim_same ? 1 : (_choice_count(pi.dim[1]) * _choice_count(pi.dim[2]))
			i_args *= pi_args_same ? 1 : _args_choice_count(pi.args)
			i_mat  *= pi_mat_same ? 1 : _choice_count(pi.material)
		end
		insu_factor = i_dim * i_args * i_mat

		total *= cond_factor * insu_factor
	end
	return total
end

Base.length(cbs::CableBuilderSpec) = cardinality(cbs)

function Base.show(io::IO, ::MIME"text/plain", cbs::CableBuilderSpec)
	comp_names = unique(p.component for p in cbs.parts)
	by_comp = Dict{Symbol, Vector{PartSpec}}()
	for p in cbs.parts
		get!(by_comp, p.component, PartSpec[]) |> v -> push!(v, p)
	end

	println(io, "CableBuilderSpec(\"", cbs.cable_id, "\")")
	println(io, "  components: ", join(string.(comp_names), ", "))

	total = 1
	for cname in comp_names
		ps = by_comp[cname]
		cond = [p for p in ps if p.part_type <: DataModel.AbstractConductorPart]
		insu = [p for p in ps if p.part_type <: DataModel.AbstractInsulatorPart]
		isempty(cond) && Base.error("component '$cname' has no conductors")
		isempty(insu) && Base.error("component '$cname' has no insulators")

		# conductors: couple to first when tuples compare equal
		p1c    = cond[1]
		c_dim  = _choice_count(p1c.dim[1]) * _choice_count(p1c.dim[2])
		c_args = _args_choice_count(p1c.args)
		c_mat  = _choice_count(p1c.material)
		for pc in cond[2:end]
			c_dim  *= (pc.dim == p1c.dim) ? 1 : (_choice_count(pc.dim[1]) * _choice_count(pc.dim[2]))
			c_args *= (pc.args == p1c.args) ? 1 : _args_choice_count(pc.args)
			c_mat  *= (pc.material == p1c.material) ? 1 : _choice_count(pc.material)
		end
		cond_factor = c_dim * c_args * c_mat

		# insulators: same coupling rule vs first insulator
		p1i    = insu[1]
		i_dim  = _choice_count(p1i.dim[1]) * _choice_count(p1i.dim[2])
		i_args = _args_choice_count(p1i.args)
		i_mat  = _choice_count(p1i.material)
		for pi in insu[2:end]
			i_dim  *= (pi.dim == p1i.dim) ? 1 : (_choice_count(pi.dim[1]) * _choice_count(pi.dim[2]))
			i_args *= (pi.args == p1i.args) ? 1 : _args_choice_count(pi.args)
			i_mat  *= (pi.material == p1i.material) ? 1 : _choice_count(pi.material)
		end
		insu_factor = i_dim * i_args * i_mat

		fac = cond_factor * insu_factor
		total *= fac
		print(io, "  • ", cname, ": ")
		print(io, "cond(dim=", c_dim, ", args=", c_args, ", mat=", c_mat, "); ")
		println(io, "insu(dim=", i_dim, ", args=", i_args, ", mat=", i_mat, ")  ⇒ ×", fac)
	end

	println(io, "  cardinality: ", total)
	if cbs.nominal !== nothing
		println(io, "  nominal: ", typeof(cbs.nominal))
	end
end


function show_trace(tr::DesignTrace)
	println("Design: ", tr.cable_id)
	for comp in tr.components
		println("  Component: ", comp.name)
		for c in comp.choices
			mat = c.mat
			println("    [", c.role, "] ", c.T,
				" layers=", c.layers,
				" dim=", c.dim,
				" args=", c.args,
				" ρ=", mat.rho, " εr=", mat.eps_r, " μr=", mat.mu_r)
		end
	end
end

Base.show(io::IO, ::MIME"text/plain", tr::DesignTrace) = show_trace(tr)

# == Choice-count helpers (reusing ParametricBuilder counting) ==

_position_choice_count(p::PositionSpec) = _choice_count(p.dx) * _choice_count(p.dy)

_earth_choice_count(e::EarthSpec) =
	_choice_count(e.rho) * _choice_count(e.eps_r) * _choice_count(e.mu_r) *
	_choice_count(e.t)

# == Public cardinality API ==
function cardinality(spec::SystemBuilderSpec)
	pos_prod =
		isempty(spec.positions) ? 0 :
		prod(_position_choice_count(p) for p in spec.positions)
	pos_prod == 0 && return 0
	return cardinality(spec.builder) *
		   _choice_count(spec.length) *
		   pos_prod *
		   _choice_count(spec.temperature) *
		   _earth_choice_count(spec.earth)
end

Base.length(spec::SystemBuilderSpec) = cardinality(spec)

# == Iterator over fully formed LineParametersProblem (skips overlaps silently) ==
function Base.iterate(spec::SystemBuilderSpec)
	ch = iterate_problems(spec)
	try
		x = take!(ch);
		return (x, ch)
	catch e
		@error "SystemBuilderSpec iteration failed before first yield" exception=(
			e,
			catch_backtrace(),
		)

		rethrow()
	end
end

function Base.iterate(::SystemBuilderSpec, ch::Channel{LineParametersProblem})
	try
		x = take!(ch);
		return (x, ch)
	catch
		return nothing
	end
end

Base.IteratorEltype(::Type{SystemBuilderSpec}) = Base.HasEltype()
Base.eltype(::Type{SystemBuilderSpec}) = LineParametersProblem
Base.IteratorSize(::Type{SystemBuilderSpec}) = Base.SizeUnknown()

# == Terse pretty printer (because why the hell not?) ==
# show at most `limit` values: "v1, v2, ..., vN (N=total)"
_fmt_vals(vals; limit = 8) = begin
	v = collect(vals)
	n = length(v)
	if n == 0
		"∅"
	elseif n <= limit
		string(join(v, ", "))
	else
		string(join(v[1:limit], ", "), ", …  (N=", n, ")")
	end
end

# expand one knob (your (valuespec,pct) grammar) into concrete values
_vals_pair(p) = collect(_expand_pair(p))                      # :contentReference[oaicite:2]{index=2}
# axis around anchor (handles (nothing, pct) → uncertain anchor)
_vals_axis(anchor, dspec) = collect(_axis(anchor, dspec))     # :contentReference[oaicite:3]{index=3}

# deterministic freq summary: list if tiny, else min..max (N)
_fmt_freqs(f::AbstractVector) =
	length(f) ≤ 8 ? join(f, ", ") :
	string(first(f), " … ", last(f), "  (N=", length(f), ")")

# stable, human order for phases: core,sheath,jacket first if present, then alphabetical
function _fmt_map(conn::Dict{String, Int})
	prio = Dict("core"=>1, "sheath"=>2, "jacket"=>3)
	ks = collect(keys(conn))
	sort!(ks, by = k -> (get(prio, k, 1000), k))
	# FIX: use getindex, not get
	vs = getindex.(Ref(conn), ks)          # or: map(k -> conn[k], ks)
	return join(string.(ks, "=>", vs), ", ")
end

function Base.show(io::IO, ::MIME"text/plain", spec::SystemBuilderSpec)
	println(io, "SystemBuilder(\"", spec.system_id, "\")")
	println(io, "  designs × = ", cardinality(spec.builder))  # we don’t explode details here; PB has its own show. :contentReference[oaicite:4]{index=4}

	# positions block
	println(io, "  positions = ", length(spec.positions))
	for (i, p) in enumerate(spec.positions)
		dxvals = _vals_axis(p.x0, p.dx)
		dyvals = _vals_axis(p.y0, p.dy)
		println(io, "    • p", i, "  x: ", _fmt_vals(dxvals), ", y: ", _fmt_vals(dyvals),
			", phases: {", _fmt_map(p.conn), "}")
	end

	# system scalars
	println(io, "  length  = ", _fmt_vals(_vals_pair(spec.length)))
	println(io, "  temp    = ", _fmt_vals(_vals_pair(spec.temperature)))

	# earth knobs (each axis separately)
	println(io, "  earth:")
	println(io, "    ρ     = ", _fmt_vals(_vals_pair(spec.earth.rho)))
	println(io, "    εr    = ", _fmt_vals(_vals_pair(spec.earth.eps_r)))
	println(io, "    μr    = ", _fmt_vals(_vals_pair(spec.earth.mu_r)))
	println(io, "    t     = ", _fmt_vals(_vals_pair(spec.earth.t)))

	# frequencies (deterministic vector coming from the user/spec)
	if hasfield(SystemBuilderSpec, :frequencies) && !isempty(getfield(spec, :frequencies))
		f = getfield(spec, :frequencies)
		println(io, "  f       = ", _fmt_freqs(f))
	end

	println(io, "  cardinality (upper bound): ", cardinality(spec))
end
