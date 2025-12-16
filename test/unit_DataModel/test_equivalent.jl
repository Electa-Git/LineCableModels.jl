@testsnippet simplify_fixtures begin
	# Aliases
	const LM = LineCableModels
	const DM = LM.DataModel
	const MAT = LM.Materials
	using Measurements: measurement

	# Basic materials
	copper_props = MAT.Material(1.7241e-8, 1.0, 1.0, 20.0, 0.00393, 401.0)
	xlpe_props = MAT.Material(1e10, 2.3, 1.0, 20.0, 0.0, 0.3)
	semi_props = MAT.Material(1e3, 2.6, 1.0, 20.0, 0.0, 148.0)

	# Geometry helpers
	d_wire = 3e-3
	rin0 = 0.0

	function make_conductor_group()
		core = DM.WireArray(rin0, DM.Diameter(d_wire), 1, 0.0, copper_props)
		g = DM.ConductorGroup(core)
		add!(g, DM.WireArray, DM.Diameter(d_wire), 6, 10.0, copper_props)
		add!(g, DM.Strip, DM.Thickness(0.5e-3), 0.02, 8.0, copper_props)
		add!(g, DM.Tubular, DM.Thickness(0.8e-3), copper_props)
		g
	end

	function make_insulator_group(conductor_group)
		ins1 = DM.Insulator(conductor_group.radius_ext, DM.Thickness(2.0e-3), xlpe_props)
		ig = DM.InsulatorGroup(ins1)
		add!(ig, DM.Semicon, DM.Thickness(0.8e-3), semi_props)
		add!(ig, DM.Insulator, DM.Thickness(2.0e-3), xlpe_props)
		ig
	end

	function make_component(id::AbstractString)
		g = make_conductor_group()
		ig = make_insulator_group(g)
		DM.CableComponent(String(id), g, ig)
	end

	function make_design(id::AbstractString; ncomponents::Int = 1)
		comps = [make_component(n == 1 ? "core" : "comp$(n)") for n in 1:ncomponents]
		des = DM.CableDesign(String(id), comps[1])
		for c in comps[2:end]
			add!(des, c.id, c.conductor_group, c.insulator_group)
		end
		des
	end
end

@testitem "simplify unit tests" setup =
	[defaults, deps_datamodel, defs_materials, simplify_fixtures] begin
	const DM = LineCableModels.DataModel

	@testset "Input Validation" begin
		des = make_design("CAB-V-0")

		# Missing required positional argument
		@test_throws MethodError DM.equivalent()

		# Invalid first argument type
		@test_throws MethodError DM.equivalent(42)

		# Invalid keyword type for new_id
		@test_throws TypeError DM.equivalent(des; new_id = 123)


	end

	@testset "Basic Functionality" begin
		des = make_design("CAB-BASIC"; ncomponents = 2)
		des_s = DM.equivalent(des)

		@test des_s isa DM.CableDesign
		@test des_s.cable_id == "CAB-BASIC_equivalent"
		@test length(des_s.components) == length(des.components)

		# Geometry continuity: outer radius preserved by equivalence
		for (orig, simp) in zip(des.components, des_s.components)
			@test simp.conductor_group.radius_in ≈ orig.conductor_group.radius_in atol =
				TEST_TOL
			@test simp.conductor_group.radius_ext ≈ orig.conductor_group.radius_ext atol =
				TEST_TOL
			@test simp.insulator_group.radius_ext ≈ orig.insulator_group.radius_ext atol =
				TEST_TOL
			@test simp.id == orig.id
		end

		# new_id override
		des_s2 = DM.equivalent(des; new_id = "CAB-SIMPLE")
		@test des_s2.cable_id == "CAB-SIMPLE"
	end

	@testset "Equivalence Preservation" begin
		# The simplified design must preserve the component-equivalent properties
		des = make_design("CAB-EQ"; ncomponents = 2)
		des_s = DM.equivalent(des)

		for (orig, simp) in zip(des.components, des_s.components)
			# Compare equivalent material properties (conductor)
			@test simp.conductor_props.rho ≈ orig.conductor_props.rho atol = TEST_TOL
			@test simp.conductor_props.eps_r ≈ orig.conductor_props.eps_r atol = TEST_TOL
			@test simp.conductor_props.mu_r ≈ orig.conductor_props.mu_r atol = TEST_TOL
			@test simp.conductor_props.T0 ≈ orig.conductor_props.T0 atol = TEST_TOL
			@test simp.conductor_props.alpha ≈ orig.conductor_props.alpha atol = TEST_TOL

			# Compare equivalent material properties (insulator)
			@test simp.insulator_props.rho ≈ orig.insulator_props.rho atol = TEST_TOL
			@test simp.insulator_props.eps_r ≈ orig.insulator_props.eps_r atol = TEST_TOL
			@test simp.insulator_props.mu_r ≈ orig.insulator_props.mu_r atol = TEST_TOL
			@test simp.insulator_props.T0 ≈ orig.insulator_props.T0 atol = TEST_TOL
			@test simp.insulator_props.alpha ≈ orig.insulator_props.alpha atol = TEST_TOL

			# Compare group lumped parameters (should be preserved by construction)
			@test simp.conductor_group.resistance ≈ orig.conductor_group.resistance atol =
				TEST_TOL
			@test simp.conductor_group.gmr ≈ orig.conductor_group.gmr atol = TEST_TOL
			@test simp.insulator_group.shunt_capacitance ≈
				  orig.insulator_group.shunt_capacitance atol = TEST_TOL
			@test simp.insulator_group.shunt_conductance ≈
				  orig.insulator_group.shunt_conductance atol = TEST_TOL
		end
	end

	@testset "Edge Cases" begin
		# Use Measurement geometry to ensure robustness with promoted numeric types
		des = make_design("CAB-EDGE")
		desM = DM.CableDesign(
			"CAB-EDGE-M",
			DM.CableComponent(
				des.components[1].id,
				DM.coerce_to_T(
					des.components[1].conductor_group,
					Measurements.Measurement{Float64},
				),
				DM.coerce_to_T(
					des.components[1].insulator_group,
					Measurements.Measurement{Float64},
				),
			),
		)

		des_sM = DM.equivalent(desM)
		@test typeof(des_sM.components[1].conductor_group.radius_in) <:
			  Measurements.Measurement
		@test typeof(des_sM.components[1].conductor_group.radius_ext) <:
			  Measurements.Measurement
		@test typeof(des_sM.components[1].insulator_group.radius_in) <:
			  Measurements.Measurement
		@test typeof(des_sM.components[1].insulator_group.radius_ext) <:
			  Measurements.Measurement
	end

	@testset "Physical Behavior" begin
		des = make_design("CAB-PHYS")
		des_s = DM.equivalent(des)
		for comp in des_s.components
			@test comp.conductor_props.rho > 0
			@test comp.conductor_group.gmr > 0
			@test comp.insulator_props.eps_r > 0
			@test comp.insulator_group.shunt_capacitance > 0
		end
	end

	@testset "Type Stability & Promotion" begin
		des = make_design("CAB-TYPES")
		cF = des.components[1]

		# Base: Float64 -> Float64
		desF = DM.CableDesign("CAB-F", cF)
		sF = DM.equivalent(desF)
		@test eltype([sF.components[1].conductor_group.radius_in]) == Float64

		# Fully promoted: Measurement -> Measurement
		gM = DM.coerce_to_T(cF.conductor_group, Measurements.Measurement{Float64})
		igM = DM.coerce_to_T(cF.insulator_group, Measurements.Measurement{Float64})
		desM = DM.CableDesign("CAB-M", DM.CableComponent("coreM", gM, igM))
		sM = DM.equivalent(desM)
		@test typeof(sM.components[1].conductor_group.radius_in) <: Measurements.Measurement
		@test typeof(sM.components[1].insulator_group.radius_ext) <:
			  Measurements.Measurement

		# Mixed cases
		desC = DM.CableDesign("CAB-C", DM.CableComponent("c1", gM, cF.insulator_group))
		sC = DM.equivalent(desC)
		@test typeof(sC.components[1].conductor_group.radius_in) <: Measurements.Measurement
		@test typeof(sC.components[1].insulator_group.radius_in) <: Measurements.Measurement

		desI = DM.CableDesign("CAB-I", DM.CableComponent("c2", cF.conductor_group, igM))
		sI = DM.equivalent(desI)
		@test typeof(sI.components[1].conductor_group.radius_in) <: Measurements.Measurement
		@test typeof(sI.components[1].insulator_group.radius_in) <: Measurements.Measurement
	end
end
