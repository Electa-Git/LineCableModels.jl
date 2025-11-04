# spec: (value_spec, pct_spec) — pct_spec can be `nothing | number | vector | (lo,hi,n)`
"""
PartSpec:
- component::Symbol           # e.g. :core, :sheath, :jacket
- part_type::Type             # WireArray, Tubular, Strip, Insulator, Semicon, …
- n_layers::Int               # how many stacked layers of this part_type
- dim::Tuple                  # diameter OR thickness OR radius (spec, pct)
- args::Tuple                 # specialized ctor positional args
- material::MaterialSpec
"""
struct PartSpec
	component::Symbol
	part_type::Type
	n_layers::Int
	dim::Tuple # (spec, pct)
	args::Tuple # positional args; each entry is either a number or (spec, pct)
	material::MaterialSpec
end

PartSpec(component::Symbol, part_type::Type, n_layers::Int;
	dim, args = (), material::MaterialSpec) =
	PartSpec(component, part_type, n_layers, dim, args, material)

"""
CableBuilderSpec:
- cable_id::String
- parts::Vector{PartSpec}  # may interleave conductor/insulator arbitrarily
- nominal::Union{Nothing,DataModel.NominalData}
"""
struct CableBuilderSpec
	cable_id::String
	parts::Vector{PartSpec}
	nominal::Union{Nothing, DataModel.NominalData}
end
CableBuilder(id::AbstractString, parts::Vector{PartSpec}; nominal = nothing) =
	CableBuilderSpec(String(id), parts, nominal)

# --- minimal flattening helpers (accept PartSpec or collections of them) -----
function _collect_parts!(acc::Vector{PartSpec}, x)
	if x isa PartSpec
		push!(acc, x)
	elseif x isa AbstractVector
		@inbounds for y in x
			_collect_parts!(acc, y)
		end
	elseif isnothing(x)
		@warn "Ignoring `nothing` in parts collection."
	else
		Base.error("Expected PartSpec or a collection of PartSpec; got $(typeof(x))")
	end
	return acc
end

# ctor that accepts a vector with possible nested vectors (no splat needed)
function CableBuilder(id::AbstractString, parts_any::AbstractVector; nominal = nothing)
	acc = PartSpec[]
	_collect_parts!(acc, parts_any)
	return CableBuilderSpec(String(id), acc, nominal)  # calls your primary ctor
end

# ctor that accepts varargs (mixed PartSpec and vectors), plus nominal kw
function CableBuilder(id::AbstractString, parts...; nominal = nothing)
	acc = PartSpec[]
	@inbounds for p in parts
		_collect_parts!(acc, p)
	end
	return CableBuilderSpec(String(id), acc, nominal)
end

struct PartChoice
	idx::Int                  # index in ps vector (1-based)
	role::Symbol              # :conductor or :insulator
	T::Type
	dim::Any                  # chosen scalar (Diameter/Thickness proxy input)
	args::Tuple               # chosen positional args (scalars)
	mat::Materials.Material   # concrete material used
	layers::Int               # n_layers replicated with that choice
end

struct ComponentTrace
	name::String
	choices::Vector{PartChoice}
end

struct DesignTrace
	cable_id::String
	components::Vector{ComponentTrace}
end



# ----- anchor: last physical layer, not the container -----
@inline _anchor(x::Real) = x
@inline _anchor(x::DataModel.AbstractConductorPart) = x
@inline _anchor(x::DataModel.AbstractInsulatorPart) = x

@inline function _anchor(g::DataModel.ConductorGroup)
	L = getfield(g, :layers)
	@assert !isempty(L) "ConductorGroup has no layers to anchor on."
	return L[end]
end
@inline function _anchor(g::DataModel.InsulatorGroup)
	L = getfield(g, :layers)
	@assert !isempty(L) "InsulatorGroup has no layers to anchor on."
	return L[end]
end

# ----- proxy for radius_ext by CONTRACT -----
@inline function _resolve_dim(T::Type, is_abs_first::Bool)
	return T <: DataModel.AbstractWireArray ? :diameter :
		   (is_abs_first && T === DataModel.Tubular ? :diameter : :thickness)
end

@inline _make_dim(::Val{:diameter}, d)  = DataModel.Diameter(d)
@inline _make_dim(::Val{:thickness}, d) = DataModel.Thickness(d)
@inline _make_dim(::Val{:radius}, r)    = r  # if direct radius
@inline _make_dim(sym::Symbol, d)       = _make_dim(Val(sym), d)


function _init_cg(T::Type, base, dim_val, args_pos::Tuple, mat; abs_first::Bool)

	r_in = _anchor(base)
	sym  = _resolve_dim(T, abs_first)
	_    = _make_dim(sym, dim_val)   # keeps intent (WireArray ignores this)

	if T <: DataModel.AbstractWireArray
		@assert length(args_pos) ≥ 1 "WireArray needs (n, [lay])."
		n   = args_pos[1]
		lay = length(args_pos) ≥ 2 ? args_pos[2] : 0.0
		return DataModel.ConductorGroup(
			DataModel.WireArray(r_in, DataModel.Diameter(dim_val), n, lay, mat),
		)
	else
		return DataModel.ConductorGroup(T(r_in, _make_dim(sym, dim_val), args_pos..., mat))
	end
end


function _add_conductor!(
	cg::DataModel.ConductorGroup,
	T::Type,
	dim_val,
	args_pos::Tuple,
	mat;
	layer::Int,
)
	if T <: DataModel.AbstractWireArray
		@assert length(args_pos) ≥ 1 "WireArray needs (n, [lay])."
		n   = args_pos[1]
		lay = length(args_pos) ≥ 2 ? args_pos[2] : 0.0
		add!(cg, DataModel.WireArray, DataModel.Diameter(dim_val), layer*n, lay, mat)
	else
		# thickness by contract for all non-wire additions
		add!(cg, T, _make_dim(:thickness, dim_val), args_pos..., mat)
	end
end


function _init_ig(T::Type, cg::DataModel.ConductorGroup, dim_val, args_pos::Tuple, mat)
	c_last = _anchor(cg)
	obj = T(c_last, DataModel.Thickness(dim_val), args_pos..., mat)  # insulators use THICKNESS
	return DataModel.InsulatorGroup(obj)
end


function _add_insulator!(
	ig::DataModel.InsulatorGroup,
	T::Type,
	dim_val,
	args_pos::Tuple,
	mat,
)
	add!(ig, T, DataModel.Thickness(dim_val), args_pos..., mat)
end

# Build all variants of ONE component, anchored at `base` (0.0 for the very first)
function _make_variants(ps::Vector{PartSpec}, base)
	cond = [p for p in ps if p.part_type <: DataModel.AbstractConductorPart]
	insu = [p for p in ps if p.part_type <: DataModel.AbstractInsulatorPart]
	isempty(cond) && error("component has no conductors")
	isempty(insu) && error("component has no insulators")

	variants = Tuple{DataModel.CableComponent, DataModel.InsulatorGroup, ComponentTrace}[]

	# ---------------- first conductor choice spaces ----------------
	p1c    = cond[1]
	mats1  = _make_range(p1c.material)
	dims1  = _make_range(p1c.dim[1]; pct = p1c.dim[2])
	args1s = collect(_expand_args(p1c.args))  # Vector{<:Tuple}

	# remaining conductors — spaces, with COUPLING flags to p1c
	# Tuple layout: (pc, mcs_or_nothing, dcs_or_nothing, acs_or_nothing)
	rest_cond_spaces = Tuple{PartSpec, Any, Union{Nothing, Any}, Union{Nothing, Any}}[]
	for pc in cond[2:end]
		same_mat  = (pc.material == p1c.material)
		same_dim  = (pc.dim == p1c.dim)
		same_args = (pc.args == p1c.args)

		mcs = same_mat ? nothing : _make_range(pc.material)
		dcs = same_dim ? nothing : _make_range(pc.dim[1]; pct = pc.dim[2])
		acs = same_args ? nothing : collect(_expand_args(pc.args))

		push!(rest_cond_spaces, (pc, mcs, dcs, acs))
	end

	# ---------------- first insulator choice spaces ----------------
	p1i    = insu[1]
	matsi  = _make_range(p1i.material)
	dimsi  = _make_range(p1i.dim[1]; pct = p1i.dim[2])
	args1i = collect(_expand_args(p1i.args))

	# remaining insulators — spaces, with COUPLING flags to p1i
	rest_ins_spaces = Tuple{PartSpec, Any, Union{Nothing, Any}, Union{Nothing, Any}}[]
	for pi in insu[2:end]
		same_mat  = (pi.material == p1i.material)
		same_dim  = (pi.dim == p1i.dim)
		same_args = (pi.args == p1i.args)

		m2 = same_mat ? nothing : _make_range(pi.material)
		d2 = same_dim ? nothing : _make_range(pi.dim[1]; pct = pi.dim[2])
		a2 = same_args ? nothing : collect(_expand_args(pi.args))

		push!(rest_ins_spaces, (pi, m2, d2, a2))
	end

	# ---------------- selection stacks (resolved tuples) -------------
	# chosen_c stores (pc, mc, dc, ac)
	# chosen_i stores (pi, m2i, d2i, a2i)
	chosen_c = Vector{NTuple{4, Any}}()
	chosen_i = Vector{NTuple{4, Any}}()

	# ---------------- build with current resolved choices ------------
	function build_with_current_selection(mat1, d1, a1, mi, di, ai)
		# 1) conductors
		cg = _init_cg(p1c.part_type, base, d1, a1, mat1; abs_first = base == 0.0)
		for k in 2:p1c.n_layers
			_add_conductor!(cg, p1c.part_type, d1, a1, mat1; layer = k)
		end
		for (pc, mc, dc, ac) in chosen_c
			for k in 1:pc.n_layers
				_add_conductor!(cg, pc.part_type, dc, ac, mc; layer = k)
			end
		end

		# 2) insulators
		ig = _init_ig(p1i.part_type, cg, di, ai, mi)
		for (pi, m2i, d2i, a2i) in chosen_i
			for k in 1:pi.n_layers
				_add_insulator!(ig, pi.part_type, d2i, a2i, m2i)
			end
		end

		# assemble trace
		choices = PartChoice[]
		# first conductor spec
		push!(choices, PartChoice(1, :conductor, p1c.part_type, d1, a1, mat1, p1c.n_layers))
		# remaining conductors
		for (j, (pc, mc, dc, ac)) in enumerate(chosen_c)
			push!(
				choices,
				PartChoice(1 + j, :conductor, pc.part_type, dc, ac, mc, pc.n_layers),
			)
		end
		# first insulator spec
		push!(
			choices,
			PartChoice(
				length(choices)+1,
				:insulator,
				p1i.part_type,
				di,
				ai,
				mi,
				p1i.n_layers,
			),
		)
		# remaining insulators
		for (pi, m2i, d2i, a2i) in chosen_i
			push!(
				choices,
				PartChoice(
					length(choices)+1,
					:insulator,
					pi.part_type,
					d2i,
					a2i,
					m2i,
					pi.n_layers,
				),
			)
		end
		ctrace = ComponentTrace(String(ps[1].component), choices)

		push!(
			variants,
			(DataModel.CableComponent(String(ps[1].component), cg, ig), ig, ctrace),
		)

		# push!(variants, (DataModel.CableComponent(String(ps[1].component), cg, ig), ig))
	end

	# ---------------- enumerate insulators with coupling -------------
	function choose_ins(idx::Int, mi, di, ai, mat1, d1, a1)
		if idx > length(rest_ins_spaces)
			build_with_current_selection(mat1, d1, a1, mi, di, ai)
			return
		end
		pi, m2, d2, a2 = rest_ins_spaces[idx]

		Ms = (m2 === nothing) ? (mi,) : m2
		Ds = (d2 === nothing) ? (di,) : d2
		As = (a2 === nothing) ? (ai,) : a2

		for m2i in Ms, d2i in Ds, a2i in As
			push!(chosen_i, (pi, m2i, d2i, a2i))
			choose_ins(idx + 1, mi, di, ai, mat1, d1, a1)
			pop!(chosen_i)
		end
	end

	# ---------------- enumerate conductors with coupling -------------
	function choose_cond(idx::Int, mat1, d1, a1, mi, di, ai)
		if idx > length(rest_cond_spaces)
			empty!(chosen_i)
			choose_ins(1, mi, di, ai, mat1, d1, a1)
			return
		end
		pc, mcs, dcs, acs = rest_cond_spaces[idx]

		Ms = (mcs === nothing) ? (mat1,) : mcs
		Ds = (dcs === nothing) ? (d1,) : dcs
		As = (acs === nothing) ? (a1,) : acs

		for mc in Ms, dc in Ds, ac in As
			push!(chosen_c, (pc, mc, dc, ac))
			choose_cond(idx + 1, mat1, d1, a1, mi, di, ai)
			pop!(chosen_c)
		end
	end

	# ---------------- top-level selection loops ----------------------
	for mat1 in mats1, d1 in dims1, a1 in args1s
		for mi in matsi, di in dimsi, ai in args1i
			empty!(chosen_c)
			empty!(chosen_i)
			choose_cond(1, mat1, d1, a1, mi, di, ai)
		end
	end

	return variants
end


function build(cbs::CableBuilderSpec; trace::Bool = false)
	comp_names = unique(p.component for p in cbs.parts)
	by_comp = Dict{Symbol, Vector{PartSpec}}()
	for p in cbs.parts
		get!(by_comp, p.component, PartSpec[]) |> v -> push!(v, p)
	end

	# partials: (built_components, last_ig_or_nothing)
	partials = Tuple{
		Vector{DataModel.CableComponent},
		Union{Nothing, DataModel.InsulatorGroup},
		Vector{ComponentTrace},
	}[(DataModel.CableComponent[], nothing, ComponentTrace[])]

	for cname in comp_names
		ps = by_comp[cname]
		new_partials = Tuple{
			Vector{DataModel.CableComponent},
			Union{Nothing, DataModel.InsulatorGroup},
			Vector{ComponentTrace},
		}[]
		for (built, last_ig, tr) in partials
			base = last_ig === nothing ? 0.0 : last_ig
			for (comp, ig, ctrace) in _make_variants(ps, base)
				push!(new_partials, (vcat(built, comp), ig, [tr...; ctrace]))
			end
		end
		partials = new_partials
	end



	if !trace
		designs = DataModel.CableDesign[]
		for (comps, _) in ((x[1], x[2]) for x in partials)
			des = DataModel.CableDesign(cbs.cable_id, comps[1]; nominal_data = cbs.nominal)
			for k in Iterators.drop(eachindex(comps), 1)
				add!(des, comps[k])
			end
			push!(designs, des)
		end
		return designs
	else
		designs = DataModel.CableDesign[]
		traces  = DesignTrace[]
		for (comps, _, ctraces) in partials
			des = DataModel.CableDesign(cbs.cable_id, comps[1]; nominal_data = cbs.nominal)
			for k in 2:length(comps)
				;
				add!(des, comps[k]);
			end
			push!(designs, des)
			push!(traces, DesignTrace(cbs.cable_id, ctraces))
		end
		return designs, traces
	end
end

"""
	iterate_designs(cbs) -> Channel{DataModel.CableDesign}

Lazy stream of `CableDesign`s built from `CableBuilderSpec` without allocating all of them.
Works with `for d in iterate_designs(cbs)`.
"""
function iterate_designs(cbs::CableBuilderSpec)
	# group by component
	comp_names = unique(p.component for p in cbs.parts)
	by_comp = Dict{Symbol, Vector{PartSpec}}()
	for p in cbs.parts
		get!(by_comp, p.component, PartSpec[]) |> v -> push!(v, p)
	end

	return Channel{DataModel.CableDesign}(32) do ch
		built  = DataModel.CableComponent[]
		lastig = Ref{Union{Nothing, DataModel.InsulatorGroup}}(nothing)

		function dfs(i::Int)
			if i > length(comp_names)
				des = DataModel.CableDesign(
					cbs.cable_id,
					built[1];
					nominal_data = cbs.nominal,
				)
				for k in 2:length(built)
					;
					add!(des, built[k]);
				end
				put!(ch, des)
				return
			end
			cname = comp_names[i]
			ps    = by_comp[cname]
			base  = (lastig[] === nothing) ? 0.0 : lastig[]

			for (comp, ig) in _make_variants(ps, base)
				push!(built, comp)
				prev = lastig[];
				lastig[] = ig
				dfs(i + 1)
				lastig[] = prev
				pop!(built)
			end
		end

		dfs(1)
	end
end

module Conductor

using ..ParametricBuilder: PartSpec, _spec
using ...DataModel: DataModel

# wire: args are (n, lay)
Wires(component::Symbol; layers::Int, d, n::Int, lay = 11.0, m) =
	PartSpec(component, DataModel.WireArray, layers;
		dim = _spec(d), args = (n, _spec(lay)), material = m)

# tube: no extra args
Tubular(component::Symbol; layers::Int, t, m) =
	PartSpec(component, DataModel.Tubular, layers;
		dim = _spec(t), args = (), material = m)

# strip: args are (width, lay)
Strip(component::Symbol; layers::Int, t, w, lay = 0.0, m) =
	PartSpec(component, DataModel.Strip, layers;
		dim = _spec(t), args = (_spec(w), _spec(lay)), material = m)

# solid: inherits inner radius = 0.0, builds from diameter
Solid(component::Symbol; d, m) =
	PartSpec(component, DataModel.Tubular, 1;
		dim = _spec(d), args = (), material = m)

# central + hex rings sugar
function Stranded(component::Symbol; layers::Int, d, n::Int, lay = 11.0, m)
	@assert layers >= 1 "stranded: layers must be ≥ 1 (includes the central wire)."
	specs = PartSpec[]
	dspec = _spec(d)

	# 1) central wire: 1 layer, n=1, lay=0.0
	push!(
		specs,
		PartSpec(component, DataModel.WireArray, 1;
			dim = dspec, args = (1, (0.0, nothing)), material = m),
	)

	# 2) rings: (layers-1) layers, base n, common lay
	if layers > 1
		push!(
			specs,
			PartSpec(component, DataModel.WireArray, layers - 1;
				dim = dspec, args = (n, _spec(lay)), material = m),
		)
	end

	return specs
end

end

module Insulator

using ..ParametricBuilder: PartSpec, _spec
using ...DataModel: DataModel

Tubular(component::Symbol; layers::Int, t, m) =
	PartSpec(component, DataModel.Insulator, layers;
		dim = _spec(t), args = (), material = m)

Semicon(component::Symbol; layers::Int, t, m) =
	PartSpec(component, DataModel.Semicon, layers;
		dim = _spec(t), args = (), material = m)
end
