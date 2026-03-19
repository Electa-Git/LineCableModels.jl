@testitem "BaseParams: calc_shunt_capacitance unit tests" setup = [defaults] begin
	@testset "Basic Functionality" begin
		# Example from docstring
		r_in = 0.01
		r_ex = 0.02
		epsr = 2.3
		cap = calc_shunt_capacitance(r_in, r_ex, epsr)
		@test isapprox(cap, 1.241e-10, atol = TEST_TOL)
		# Vacuum (epsr = 1)
		cap_vac = calc_shunt_capacitance(0.01, 0.02, 1.0)
		@test cap_vac < cap
	end

	@testset "Edge Cases" begin
		# Collapsing geometry: r_in -> r_ex
		cap = calc_shunt_capacitance(0.02, 0.02, 2.3)
		@test isinf(cap) || isnan(cap)
		# Very large radii
		cap = calc_shunt_capacitance(1e2, 1e3, 2.3)
		@test isfinite(cap)
		# Inf/NaN input
		@test isnan(calc_shunt_capacitance(NaN, 0.02, 2.3))
		@test isnan(calc_shunt_capacitance(0.01, NaN, 2.3))
		@test isnan(calc_shunt_capacitance(0.01, 0.02, NaN))
	end

	@testset "Numerical Consistency" begin
		# Float32 vs Float64
		cap_f32 = calc_shunt_capacitance(Float32(0.01), Float32(0.02), Float32(2.3))
		cap_f64 = calc_shunt_capacitance(0.01, 0.02, 2.3)
		@test isapprox(cap_f32, cap_f64, atol = TEST_TOL)
	end

	@testset "Physical Behavior" begin
		# Capacitance increases with epsr
		c1 = calc_shunt_capacitance(0.01, 0.02, 2.3)
		c2 = calc_shunt_capacitance(0.01, 0.02, 3.0)
		@test c2 > c1
		# Capacitance decreases as radii get closer
		c3 = calc_shunt_capacitance(0.01, 0.011, 2.3)
		@test c3 > c1
	end

	@testset "Type Stability & Promotion" begin
		using Measurements
		r_in = 0.01
		r_ex = 0.02
		epsr = 2.3
		min = measurement(r_in, 1e-4)
		mex = measurement(r_ex, 1e-4)
		mepsr = measurement(epsr, 1e-2)
		# All Float64
		res1 = calc_shunt_capacitance(r_in, r_ex, epsr)
		@test typeof(res1) == Float64
		# All Measurement
		res2 = calc_shunt_capacitance(min, mex, mepsr)
		@test res2 isa Measurement{Float64}
		# Mixed: first argument Measurement
		res3 = calc_shunt_capacitance(min, r_ex, epsr)
		@test res3 isa Measurement{Float64}
		# Mixed: second argument Measurement
		res4 = calc_shunt_capacitance(r_in, mex, epsr)
		@test res4 isa Measurement{Float64}
		# Mixed: third argument Measurement
		res5 = calc_shunt_capacitance(r_in, r_ex, mepsr)
		@test res5 isa Measurement{Float64}
	end

	@testset "Uncertainty Quantification" begin
		using Measurements
		min = measurement(0.01, 1e-4)
		mex = measurement(0.02, 1e-4)
		mepsr = measurement(2.3, 1e-2)
		cap = calc_shunt_capacitance(min, mex, mepsr)
		@test cap isa Measurement{Float64}
		@test uncertainty(cap) > 0
	end
end
