@testitem "BaseParams: calc_tubular_gmr unit tests" setup = [defaults] begin
	using Measurements: measurement, value, uncertainty

	@testset "Basic Functionality" begin
		# Example from docstring
		r_ex = 0.02
		r_in = 0.01
		mu_r = 1.0
		gmr = calc_tubular_gmr(r_ex, r_in, mu_r)
		# Manual calculation for expected value
		term1 = (r_in^4 / (r_ex^2 - r_in^2)^2) * log(r_ex / r_in)
		term2 = (3 * r_in^2 - r_ex^2) / (4 * (r_ex^2 - r_in^2))
		Lin = (μ₀ * mu_r / (2 * π)) * (term1 - term2)
		expected = exp(log(r_ex) - (2 * π / μ₀) * Lin)
		@test isapprox(gmr, expected; atol = TEST_TOL)
		@test gmr > 0
		@test_throws ArgumentError calc_tubular_gmr(r_in, r_ex, mu_r)
		@test_throws ArgumentError calc_tubular_gmr(0.0, r_in, mu_r)
	end

	@testset "Edge Cases" begin
		# Thin shell: r_ex ≈ r_in
		r_ex = 0.01
		r_in = 0.01
		mu_r = 1.0
		gmr = calc_tubular_gmr(r_ex, r_in, mu_r)
		@test isapprox(gmr, r_ex; atol = TEST_TOL)

		# Infinitely thick tube: r_in ≫ 0, r_in / r_ex ≈ 0
		r_ex = 1.0
		r_in = 1e-12
		mu_r = 1.0
		gmr = calc_tubular_gmr(r_ex, r_in, mu_r)
		@test isapprox(gmr, 0.7788; atol = 1e-4)

		# r_in = 0 (solid cylinder)
		r_ex = 0.02
		r_in = 0.0
		mu_r = 1.0
		gmr = calc_tubular_gmr(r_ex, r_in, mu_r)
		@test isapprox(gmr, 0.7788 * r_ex; atol = 1e-4)

		# r_ex < r_in (should throw)
		r_ex = 0.01
		r_in = 0.02
		mu_r = 1.0
		@test_throws ArgumentError calc_tubular_gmr(r_ex, r_in, mu_r)
	end

	@testset "Numerical Consistency" begin
		# Float64
		gmr1 = calc_tubular_gmr(0.02, 0.01, 1.0)
		# Measurement{Float64}
		gmr2 = calc_tubular_gmr(
			measurement(0.02, 1e-4),
			measurement(0.01, 1e-4),
			measurement(1.0, 0.01),
		)
		@test isapprox(value(gmr2), gmr1; atol = TEST_TOL)
		@test uncertainty(gmr2) > 0
	end

	@testset "Physical Behavior" begin
		# GMR increases with r_ex
		gmr1 = calc_tubular_gmr(0.01, 0.005, 1.0)
		gmr2 = calc_tubular_gmr(0.02, 0.005, 1.0)
		@test gmr2 > gmr1
		# GMR decreases with mu_r
		gmr1 = calc_tubular_gmr(0.02, 0.01, 0.5)
		gmr2 = calc_tubular_gmr(0.02, 0.01, 2.0)
		@test gmr2 < gmr1
	end

	@testset "Type Stability & Promotion" begin
		# All Float64
		gmr = calc_tubular_gmr(0.02, 0.01, 1.0)
		@test typeof(gmr) == Float64
		# All Measurement
		gmr = calc_tubular_gmr(
			measurement(0.02, 1e-4),
			measurement(0.01, 1e-4),
			measurement(1.0, 0.01),
		)
		@test gmr isa Measurement{Float64}
		# Mixed: r_ex as Measurement
		gmr = calc_tubular_gmr(measurement(0.02, 1e-4), 0.01, 1.0)
		@test gmr isa Measurement{Float64}
		# Mixed: r_in as Measurement
		gmr = calc_tubular_gmr(0.02, measurement(0.01, 1e-4), 1.0)
		@test gmr isa Measurement{Float64}
		# Mixed: mu_r as Measurement
		gmr = calc_tubular_gmr(0.02, 0.01, measurement(1.0, 0.01))
		@test gmr isa Measurement{Float64}
	end

	@testset "Uncertainty Quantification" begin
		r_ex = measurement(0.02, 1e-4)
		r_in = measurement(0.01, 1e-4)
		mu_r = measurement(1.0, 0.01)
		gmr = calc_tubular_gmr(r_ex, r_in, mu_r)
		@test gmr isa Measurement{Float64}
		@test uncertainty(gmr) > 0
	end
end
