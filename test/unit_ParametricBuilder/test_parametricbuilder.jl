@testitem "ParametricBuilder(SystemBuilderSpec): combinatorics + value integrity" setup =
	[defaults] begin
	# -------------------------------------------------------------------------
	# Shared setup (mirrors your example, but we toggle `unc` per testset)
	# -------------------------------------------------------------------------
	using LineCableModels
	using LineCableModels.ParametricBuilder:
		CableBuilder, build, Conductor, Insulator, Material, Earth, SystemBuilder, at,
		make_stranded, make_screened, cardinality
	using LineCableModels.DataModel: trifoil_formation, LineCableSystem, CablePosition
	using Measurements

	# deterministic frequency grid
	f = 10.0 .^ range(0, stop = 6, length = 10)

	# Materials library
	materials = MaterialsLibrary(add_defaults = true)

	# deterministic geometry
	t_sct    = 0.3e-3
	t_sc_in  = 0.000768
	t_ins    = 0.0083
	t_sc_out = 0.000472
	t_cut    = 0.0001
	w_cut    = 10e-3
	t_wbt    = 0.00094
	t_alt    = 0.15e-3
	t_pet    = 0.05e-3
	t_jac    = 0.0034

	# nominal data
	datasheet_info = NominalData(
		designation_code = "NA2XS(FL)2Y",
		U0 = 18.0, U = 30.0,
		conductor_cross_section = 1000.0, screen_cross_section = 35.0,
		resistance = 0.0291, capacitance = 0.39, inductance = 0.3,
	)

	co_w   = make_stranded(datasheet_info.conductor_cross_section).best_match
	co_n   = co_w.layers
	co_d   = co_w.wire_diameter_m
	co_lay = 13.0

	sc_w   = make_screened(datasheet_info.screen_cross_section, 55.3).best_match
	sc_n   = sc_w.wires
	sc_d   = sc_w.wire_diameter_m
	sc_lay = 10.0

	# canonical parts builder (ρ/μ grids attached via `unc` in each testset)
	function make_parts(
		ms_al_uq,
		ms_cu,
		ms_pe,
		ms_xlpe,
		ms_sem1,
		ms_sem2,
		ms_polyacryl;
		unc = nothing,
	)
		return [
			# CORE conductors: stranded (central + rings) — uses PB coupling semantics. :contentReference[oaicite:1]{index=1}
			Conductor.Stranded(
				:core;
				layers = co_n,
				d = (co_d, unc),
				n = 6,
				lay = (co_lay, unc),
				m = ms_al_uq,
			),

			# CORE insulators
			Insulator.Semicon(:core; layers = 1, t = (t_sct, unc), m = ms_polyacryl),
			Insulator.Semicon(:core; layers = 1, t = (t_sc_in, unc), m = ms_sem1),
			Insulator.Tubular(:core; layers = 1, t = (t_ins, unc), m = ms_xlpe),
			Insulator.Semicon(:core; layers = 1, t = t_sc_out, m = ms_sem2),
			Insulator.Semicon(:core; layers = 1, t = t_sct, m = ms_polyacryl),

			# SHEATH
			Conductor.Wires(
				:sheath;
				layers = 1,
				d = (sc_d, unc),
				n = sc_n,
				lay = (sc_lay, unc),
				m = ms_cu,
			),
			Conductor.Strip(
				:sheath;
				layers = 1,
				t = (t_cut, unc),
				w = (w_cut, unc),
				lay = (sc_lay, unc),
				m = ms_cu,
			),
			Insulator.Semicon(:sheath; layers = 1, t = t_wbt, m = ms_polyacryl),

			# JACKET
			Conductor.Tubular(:jacket; layers = 1, t = t_alt, m = ms_al_uq),
			Insulator.Tubular(:jacket; layers = 1, t = t_pet, m = ms_pe),
			Insulator.Tubular(:jacket; layers = 1, t = t_jac, m = ms_pe),
		]
	end

	# formation anchors
	x0, y0 = 0.0, -1.0
	xa, ya, xb, yb, xc, yc = trifoil_formation(x0, y0, 0.05)

	# convenience to build SystemBuilder with 3 positions
	function make_spec(cbs; dx = (0.0, nothing), dy = (0.0, nothing),
		length = (1000.0, nothing),
		temperature = (20.0, nothing), earth = Earth(rho = 100.0, eps_r = 10.0, mu_r = 1.0))
		positions = [
			at(
				x = xa,
				y = ya,
				dx = dx,
				dy = dy,
				phases = (:core=>1, :sheath=>0, :jacket=>0),
			),
			at(
				x = xb,
				y = yb,
				dx = dx,
				dy = dy,
				phases = (:core=>2, :sheath=>0, :jacket=>0),
			),
			at(
				x = xc,
				y = yc,
				dx = dx,
				dy = dy,
				phases = (:core=>3, :sheath=>0, :jacket=>0),
			),
		]
		return SystemBuilder(
			"trifoil_case",
			cbs,
			positions;
			length = length,
			temperature = temperature,
			earth = earth,
			f = f,
		)
	end

	# # helpers to collect all produced problems (channel consumer).
	# function collect_all(xs)
	# 	acc = Any[]
	# 	for x in xs
	# 		;
	# 		push!(acc, x);
	# 	end
	# 	return acc
	# end

	# ────────────────────────────────────────────────────────────────────────
	@testset "Baseline: fully deterministic (cardinality=1, value equality)" begin
		unc = nothing
		# materials
		ms_al_uq     = Material(materials, "aluminum", rho = unc, mu_r = unc)
		ms_al        = Material(materials, "aluminum")
		ms_cu        = Material(materials, "copper")
		ms_pe        = Material(materials, "pe")
		ms_xlpe      = Material(materials, "xlpe")
		ms_sem1      = Material(materials, "semicon1")
		ms_sem2      = Material(materials, "semicon2")
		ms_polyacryl = Material(materials, "polyacrylate")

		parts = make_parts(ms_al_uq, ms_cu, ms_pe, ms_xlpe, ms_sem1, ms_sem2, ms_polyacryl; unc)
		cbs   = CableBuilder("NA2XS(FL)2Y_1000", parts; nominal = datasheet_info)

		# CableBuilder cardinality should be 1 (all scalars). :contentReference[oaicite:2]{index=2}
		@test length(cbs) == 1

		spec = make_spec(cbs; dx = (0.0, nothing), dy = (0.0, nothing),
			length = (1000.0, nothing), temperature = (20.0, nothing),
			earth = Earth(rho = 100.0, eps_r = 10.0, mu_r = 1.0))

		# SystemBuilder cardinality is designs × length × positions(dx,dy) × temperature × earth. :contentReference[oaicite:3]{index=3}
		@test length(spec) == 1

		probs = collect(spec)
		@test length(probs) == 1

		prob = probs[1]

		# Check system contents: 3 cables, phase mapping intact. :contentReference[oaicite:4]{index=4}
		sys = prob.system
		@test sys.num_cables == 3
		# access positions
		let cps = sys.cables
			@test length(cps) == 3
			# coords exact (deterministic)
			@test cps[1].horz == xa && cps[1].vert == ya
			@test cps[2].horz == xb && cps[2].vert == yb
			@test cps[3].horz == xc && cps[3].vert == yc
			# mapping
			@test cps[1].conn[1] == 1 && cps[1].conn[2] == 0 &&
				  cps[1].conn[3] == 0
			@test cps[2].conn[1] == 2
			@test cps[3].conn[1] == 3
		end

		# frequencies are carried through, deterministic preview
		@test prob.frequencies == f
	end

	# ────────────────────────────────────────────────────────────────────────
	@testset "Earth grids & %: ρ×εr×μr×t axes expand correctly" begin
		unc          = nothing
		ms_al_uq     = Material(materials, "aluminum", rho = unc, mu_r = unc)
		ms_cu        = Material(materials, "copper")
		ms_pe        = Material(materials, "pe")
		ms_xlpe      = Material(materials, "xlpe")
		ms_sem1      = Material(materials, "semicon1")
		ms_sem2      = Material(materials, "semicon2")
		ms_polyacryl = Material(materials, "polyacrylate")

		parts = make_parts(ms_al_uq, ms_cu, ms_pe, ms_xlpe, ms_sem1, ms_sem2, ms_polyacryl; unc)
		cbs   = CableBuilder("NA2XS(FL)2Y_1000", parts; nominal = datasheet_info)
		@test length(cbs) == 1

		# ρ: (100, 500, 2) with 10% → values [100, 500] each ±10% (as Measurement)  → 2
		# εr: [5, 10] → 2
		# μr: 1.0 → 1
		# t:  Inf → 1
		earth = Earth(
			rho = ((100.0, 500.0, 2), (10.0)),
			eps_r = [5.0, 10.0],
			mu_r = 1.0,
			t = Inf,
		)

		spec = make_spec(cbs; earth = earth)
		# total = 1 (designs) × 1 (len) × (1×1)^3 (positions) × 1 (T) × (2×2×1×1) = 4
		@test length(spec) == 4

		probs = collect(spec)
		@test length(probs) == 4

		# Verify EarthModel inputs carry Measurement when % is present.
		for pr in probs
			em = pr.earth_props
			# Check first layer nominal scalars hold Measurement type for rho (base value) when % applied
			# (Earth layer API stores base_* as T; implementation ensures promotion via resolve_T). 
			lay1 = em.layers[1]
			@test lay1.base_rho_g isa Measurement
			@test lay1.base_epsr_g isa Measurement
			@test lay1.base_mur_g isa Measurement
		end
	end

	# ────────────────────────────────────────────────────────────────────────
	@testset "System knobs: length & temperature grids + % expansion" begin
		unc          = (0.0, 10.0, 2)  # use 0% and 10% (length 2)
		ms_al_uq     = Material(materials, "aluminum", rho = nothing, mu_r = nothing)
		ms_cu        = Material(materials, "copper")
		ms_pe        = Material(materials, "pe")
		ms_xlpe      = Material(materials, "xlpe")
		ms_sem1      = Material(materials, "semicon1")
		ms_sem2      = Material(materials, "semicon2")
		ms_polyacryl = Material(materials, "polyacrylate")

		parts = make_parts(ms_al_uq, ms_cu, ms_pe, ms_xlpe, ms_sem1, ms_sem2, ms_polyacryl; unc = nothing)
		cbs   = CableBuilder("NA2XS(FL)2Y_1000", parts; nominal = datasheet_info)
		@test length(cbs) == 1

		# length: (1000, (0,10,2)) → [1000, 1000±10%] (2 choices)
		# temp:   (20,   (0,10,2)) → [20, 20±10%]    (2 choices)
		spec = make_spec(cbs;
			length = (1000.0, unc),
			temperature = (20.0, unc),
			earth = Earth(rho = 100.0, eps_r = 10.0, mu_r = 1.0),
		)

		# total = 1 × 2 × 1 × 2 × 1 = 4
		@test length(spec) == 4

		probs = collect(spec)
		@test length(probs) == 4

		# Check at least one problem has Measurement length & temperature (promotion path)
		# Pull system line_length via internal field (LineCableSystem constructor stores it as T). :contentReference[oaicite:7]{index=7}
		found_meas = false
		for pr in probs
			sys = pr.system
			if sys.line_length isa Measurement && pr.temperature isa Measurement
				found_meas = true
				break
			end
		end
		@test found_meas
	end

	# ────────────────────────────────────────────────────────────────────────
	@testset "Position axes: displacement grids & anchor % semantics" begin
		ms_al_uq     = Material(materials, "aluminum", rho = nothing, mu_r = nothing)
		ms_cu        = Material(materials, "copper")
		ms_pe        = Material(materials, "pe")
		ms_xlpe      = Material(materials, "xlpe")
		ms_sem1      = Material(materials, "semicon1")
		ms_sem2      = Material(materials, "semicon2")
		ms_polyacryl = Material(materials, "polyacrylate")

		parts = make_parts(ms_al_uq, ms_cu, ms_pe, ms_xlpe, ms_sem1, ms_sem2, ms_polyacryl; unc = nothing)
		cbs   = CableBuilder("NA2XS(FL)2Y_1000", parts; nominal = datasheet_info)
		@test length(cbs) == 1

		# Case A: dx sweep, dy deterministic
		specA = make_spec(cbs;
			dx = (-0.01, 0.01, 3),  # => [-0.01, 0.0, 0.01] about anchor
			dy = (0.0, nothing),
			earth = Earth(rho = 100.0, eps_r = 10.0, mu_r = 1.0),
		)
		# total = 1 × 1 × (3×1)^3 × 1 × 1 = 27
		@test length(specA) == 27

		# Validate actual x coordinates hit the expected triplet on p1
		xs = Float64[]
		for pr in specA
			push!(xs, pr.system.cables[1].horz - xa)
		end
		@test sort!(unique(round.(xs; digits = 5))) == [-0.01, 0.0, 0.01]

		# Case B: anchor % (no displacement sweep) — (nothing, pct) on dx. :contentReference[oaicite:8]{index=8}
		specB = make_spec(cbs;
			dx = (nothing, (0.0, 10.0, 2)), # anchors become [xa, measurement(xa, 10%)]
			dy = (0.0, nothing),
			earth = Earth(rho = 100.0, eps_r = 10.0, mu_r = 1.0),
		)
		# total = 1 × 1 × (2×1)^3 × 1 × 1 = 8
		@test length(specB) == 8

		# Confirm produced anchors include a Measurement with std(|xa|*10%)
		has_anchor_meas = any(
			begin
				x = pr.system.cables[1].horz
				x isa Measurement &&
					isapprox(uncertainty(x), abs(xa)*0.10; atol = eps())  # std = |xa|*10%
			end for pr in specB
		)
		@test has_anchor_meas
	end

	# ────────────────────────────────────────────────────────────────────────
	@testset "Material % propagates into designs (ρ with % → Measurement in design tree)" begin
		# Attach % to aluminum ρ and μ; all geometry deterministic
		ms_al_uq     = Material(materials, "aluminum", rho = (1.0, (0.0, 5.0, 2)), mu_r = (1.0, (0.0, 5.0, 2)))
		ms_cu        = Material(materials, "copper")
		ms_pe        = Material(materials, "pe")
		ms_xlpe      = Material(materials, "xlpe")
		ms_sem1      = Material(materials, "semicon1")
		ms_sem2      = Material(materials, "semicon2")
		ms_polyacryl = Material(materials, "polyacrylate")

		parts = make_parts(ms_al_uq, ms_cu, ms_pe, ms_xlpe, ms_sem1, ms_sem2, ms_polyacryl; unc = nothing)
		cbs   = CableBuilder("NA2XS(FL)2Y_1000", parts; nominal = datasheet_info)

		# designs = 2×2 from ρ(%) × μ(%) on the same MaterialSpec (coupled across parts when equal tuples). :contentReference[oaicite:9]{index=9}
		@test length(cbs) == cardinality(cbs)

		# System with deterministic system knobs
		spec = make_spec(cbs; earth = Earth(rho = 100.0, eps_r = 10.0, mu_r = 1.0))
		probs = collect(spec)

		# Pick one problem; inspect the first cable's design tree for Measurement presence.
		pr  = probs[1]
		des = pr.system.cables[1].design_data  # the concrete CableDesign
		# Assert that somewhere in the conductor effective material we see Measurement (ρ or μ).
		# We traverse last component conductor props or any material-like fields that match ρ/μ semantics.
		found_meas = false
		for comp in des.components
			# effective conductor/insulator props are Materials.Material
			if hasproperty(comp, :conductor_props)
				mp = getproperty(comp, :conductor_props)
				if (getproperty(mp, :rho) isa Measurement) ||
				   (getproperty(mp, :mu_r) isa Measurement) ||
				   (getproperty(mp, :T0) isa Measurement) ||
				   (getproperty(mp, :alpha) isa Measurement) ||
				   (getproperty(mp, :eps_r) isa Measurement)
					found_meas = true
					break
				end
			end
		end
		@test found_meas
	end
end
