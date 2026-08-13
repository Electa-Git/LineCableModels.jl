@testitem "Utils(coerce_to_T): Measurement covariance preservation" setup = [defaults] begin
	using Measurements
	import Measurements: derivative

	@testset "Exact-type coercion is an identity operation" begin
		x = measurement(2.0, 0.1)
		y = 3x

		@test coerce_to_T(x, typeof(x)) === x
		@test coerce_to_T(y, typeof(y)) === y
		@test derivative(coerce_to_T(y, typeof(y)), x) ≈ 3.0
	end

	@testset "Inner precision conversion retains the dependency graph" begin
		x = measurement(2.0, 0.1)
		y = 3x
		x32 = coerce_to_T(x, Measurement{Float32})
		y32 = coerce_to_T(y, Measurement{Float32})

		@test x32.tag == x.tag
		@test y32.tag == y.tag
		@test derivative(y32, x32) ≈ Float32(3.0)
	end

	@testset "Mixed BaseParams calls retain primitive sensitivities" begin
		r_ex = measurement(0.02, 0.001)
		rho = measurement(1.7241e-8, 1.0e-9)
		R = calc_tubular_resistance(0.0, r_ex, rho, 0.0, 20.0, 20.0)

		@test !iszero(derivative(R, r_ex))
		@test !iszero(derivative(R, rho))
	end
end

@testitem "DataModel: covariance survives the complete assembly path" setup = [defaults] begin
	using Measurements
	import Measurements: derivative
	copper_props = Material(1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
	insulator_props = Material(1.0e14, 2.3, 1.0, 20.0, 0.0)
	semicon_props = Material(1000.0, 1000.0, 1.0, 20.0, 0.0)

	diameter = measurement(0.02, 0.001)
	insulation_thickness = measurement(0.01, 0.001)
	semicon_thickness = measurement(0.005, 0.0005)

	core = Tubular(0.0, Diameter(diameter), copper_props)
	insulators = InsulatorGroup(
		Insulator(core, Thickness(insulation_thickness), insulator_props),
	)

	@test insulators.r_in === core.r_ex
	@test iszero(uncertainty(core.r_ex - insulators.r_in))

	insulators = add!(
		insulators,
		Semicon,
		Thickness(semicon_thickness),
		semicon_props,
	)

	expected_radius = diameter / 2 + insulation_thickness + semicon_thickness
	@test value(insulators.r_ex) ≈ value(expected_radius)
	@test uncertainty(insulators.r_ex) ≈ uncertainty(expected_radius)
	@test derivative(insulators.r_ex, diameter) ≈ 0.5
	@test derivative(insulators.r_ex, insulation_thickness) ≈ 1.0
	@test derivative(insulators.r_ex, semicon_thickness) ≈ 1.0

	conductors = ConductorGroup(core)
	component = CableComponent("core", conductors, insulators)
	design = CableDesign("uq-path", component)
	position = CablePosition(design, 0.0, -1.0, Dict("core" => 1))
	system = LineCableSystem("uq-path", 1000.0, position)
	assembled_radius =
		system.cables[1].design_data.components[1].insulator_group.r_ex

	@test system.cables[1].design_data === design
	@test derivative(assembled_radius, diameter) ≈ 0.5
	@test derivative(assembled_radius, insulation_thickness) ≈ 1.0
	@test derivative(assembled_radius, semicon_thickness) ≈ 1.0
end
