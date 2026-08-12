@testitem "Engine(ParallelRC): analytical layer and lossless limit" setup = [defaults] begin
	using LinearAlgebra

	formulation = InsulationAdmittance.ParallelRC()
	lossless = InsulationAdmittance.Lossless()

	r_in = 0.010
	r_ex = 0.018
	rho = 2.0e11
	eps_r = 2.4
	s = Complex(0.0, 2π * 50.0)

	log_ratio = log(r_ex / r_in)
	capacitance = 2π * ε₀ * eps_r / log_ratio
	conductance = 2π / (rho * log_ratio)
	p = formulation(r_in, r_ex, rho, eps_r, s)

	@test s / p ≈ conductance + s * capacitance
	@test real(s / p) > 0
	@test imag(s / p) > 0

	# Radially stacked layers are series admittances, hence their potential
	# coefficients add before the solver constructs the Maxwell matrix.
	r_mid = 0.013
	rho_2 = 8.0e10
	eps_r_2 = 3.6
	p_1 = formulation(r_in, r_mid, rho, eps_r, s)
	p_2 = formulation(r_mid, r_ex, rho_2, eps_r_2, s)
	y_1 = s / p_1
	y_2 = s / p_2
	y_series = inv(inv(y_1) + inv(y_2))
	@test s / (p_1 + p_2) ≈ y_series

	@test formulation(r_in, r_ex, Inf, eps_r, s) ≈
		lossless(r_in, r_ex, eps_r, s, 0.0)
end

@testitem "Engine(ParallelRC): physical layers reach the full EMT solve" setup = [defaults] begin
	using LinearAlgebra

	function two_terminal_problem(; uncertain = false)
		copper = Material(1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
		rho_1 = uncertain ? measurement(2.0e11, 2.0e10) : 2.0e11
		eps_1 = uncertain ? measurement(2.4, 0.12) : 2.4
		t_1 = uncertain ? measurement(3.0e-3, 1.5e-4) : 3.0e-3
		dielectric_1 = Material(rho_1, eps_1, 1.0, 20.0, 0.0)
		dielectric_2 = Material(8.0e10, 3.6, 1.0, 20.0, 0.0)
		outer_dielectric = Material(1.0e14, 2.3, 1.0, 20.0, 0.0)

		core_conductor = ConductorGroup(Tubular(0.0, Diameter(0.020), copper))
		core_insulation = InsulatorGroup(
			Insulator(core_conductor, Thickness(t_1), dielectric_1),
		)
		core_insulation = add!(
			core_insulation,
			Insulator,
			Thickness(4.0e-3),
			dielectric_2,
		)
		core = CableComponent("core", core_conductor, core_insulation)

		sheath_conductor = ConductorGroup(
			Tubular(core_insulation, Thickness(1.0e-3), copper),
		)
		sheath_insulation = InsulatorGroup(
			Insulator(sheath_conductor, Thickness(2.0e-3), outer_dielectric),
		)
		sheath = CableComponent("sheath", sheath_conductor, sheath_insulation)

		design = CableDesign("parallel-rc-test", core)
		design = add!(design, sheath)
		position = CablePosition(
			design,
			0.0,
			-1.0,
			Dict("core" => 1, "sheath" => 2),
		)
		system = LineCableSystem("parallel-rc-test", 1000.0, position)
		frequencies = [1.0e-3, 50.0, 1.0e6]
		earth = EarthModel(frequencies, 100.0, 10.0, 1.0)
		return LineParametersProblem(
			system;
			temperature = 20.0,
			earth_props = earth,
			frequencies = frequencies,
		)
	end

	formulation = FormulationSet(
		Val(:EMT);
		insulation_admittance = InsulationAdmittance.ParallelRC(),
		modal_transform = nothing,
		options = (
			reduce_bundle = false,
			kron_reduction = false,
			ideal_transposition = true,
			store_primitive_matrices = true,
		),
	)

	problem = two_terminal_problem()
	ws, line_parameters = compute!(problem, formulation)

	@test ws.insulator_layer_ranges == [1:2, 3:3]
	@test ws.r_ins_layer_in[2] ≈ ws.r_ins_layer_ext[1]
	@test size(line_parameters.Y) == (2, 2, 3)

	for k in eachindex(ws.freq)
		s = ws.jω[k]
		expected_p = sum(
			formulation.insulation_admittance(
				ws.r_ins_layer_in[layer_idx],
				ws.r_ins_layer_ext[layer_idx],
				ws.rho_ins_layer[layer_idx],
				ws.eps_ins_layer[layer_idx],
				s,
			)
			for layer_idx in ws.insulator_layer_ranges[1]
		)
		@test ws.Pin[1, 1, k] ≈ expected_p

		Y = line_parameters.Y.values[:, :, k]
		@test Y ≈ transpose(Y)
		@test all(isfinite, real.(Y))
		@test all(isfinite, imag.(Y))
		G = real.(Y)
		tolerance = 1.0e-9 * max(opnorm(G), eps())
		@test minimum(eigvals(Symmetric(G))) >= -tolerance
	end

	# Conductive leakage dominates the inner radial branch near DC, whereas
	# capacitive current dominates it at the top of the frequency sweep.
	y_near_dc = ws.jω[1] / ws.Pin[1, 1, 1]
	y_high = ws.jω[end] / ws.Pin[1, 1, end]
	@test real(y_near_dc) > imag(y_near_dc)
	@test imag(y_high) > real(y_high)

	# LEP must survive primitives -> components -> design -> system -> Z/Y.
	ws_uq, line_parameters_uq = compute!(two_terminal_problem(uncertain = true), formulation)
	y_gap_uq = ws_uq.jω[2] / ws_uq.Pin[1, 1, 2]
	@test uncertainty(real(y_gap_uq)) > 0
	@test uncertainty(imag(y_gap_uq)) > 0
	@test uncertainty(real(line_parameters_uq.Y[1, 1, 2])) > 0
	@test uncertainty(imag(line_parameters_uq.Y[1, 1, 2])) > 0
end

@testitem "UQ.mc: ParallelRC is sampled before line-parameter assembly" setup = [defaults] begin
	using Statistics
	using LineCableModels.ParametricBuilder:
		CableBuilder, Conductor, Insulator, Material, Earth, SystemBuilder, at
	import LineCableModels.UQ

	copper = Material(
		rho = 1.7241e-8,
		eps_r = 1.0,
		mu_r = 1.0,
		T0 = 20.0,
		alpha = 0.00393,
	)
	dielectric = Material(
		rho = (2.0e11, 10.0),
		eps_r = (2.4, 5.0),
		mu_r = 1.0,
		T0 = 20.0,
		alpha = 0.0,
	)
	outer_dielectric = Material(
		rho = 1.0e14,
		eps_r = 2.3,
		mu_r = 1.0,
		T0 = 20.0,
		alpha = 0.0,
	)

	parts = [
		Conductor.Solid(:core; d = 0.020, m = copper),
		Insulator.Tubular(:core; layers = 1, t = (7.0e-3, 5.0), m = dielectric),
		Conductor.Tubular(:sheath; layers = 1, t = 1.0e-3, m = copper),
		Insulator.Tubular(:sheath; layers = 1, t = 2.0e-3, m = outer_dielectric),
	]
	builder = CableBuilder("parallel-rc-mc", parts; nominal = NominalData())
	spec = SystemBuilder(
		"parallel-rc-mc",
		builder,
		at(
			x = 0.0,
			y = -1.0,
			phases = (:core => 1, :sheath => 2),
		);
		length = 1000.0,
		temperature = 20.0,
		earth = Earth(rho = 100.0, eps_r = 10.0, mu_r = 1.0),
		f = [50.0],
	)
	formulation = FormulationSet(
		Val(:EMT);
		insulation_admittance = InsulationAdmittance.ParallelRC(),
		modal_transform = nothing,
		options = (
			reduce_bundle = false,
			kron_reduction = false,
			ideal_transposition = true,
		),
	)

	result = UQ.mc(
		spec,
		formulation;
		trials = 24,
		distribution = :normal,
		seed = 20260812,
		print_step = 1000,
		return_samples = true,
	)

	@test size(result.measurements.Y) == (2, 2, 1)
	@test size(result.samples.G) == (2, 2, 1, 24)
	@test std(result.samples.G[1, 1, 1, :]) > 0
	@test std(result.samples.C[1, 1, 1, :]) > 0
	@test uncertainty(real(result.measurements.Y[1, 1, 1])) > 0
	@test uncertainty(imag(result.measurements.Y[1, 1, 1])) > 0
end
