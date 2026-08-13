@testitem "UQ.mc: joint measurement surrogate retains empirical covariance" setup = [defaults] begin
	using LinearAlgebra
	using Statistics
	using Measurements
	using LineCableModels.Commons: PhaseDomain
	import LineCableModels.UQ

	ntrials = 5
	tensor_size = (2, 2, 2)
	q = Float64[-2, -1, 0, 1, 2]
	Rsamp = Array{Float64}(undef, tensor_size..., ntrials)
	Lsamp = similar(Rsamp)
	Gsamp = similar(Rsamp)
	Csamp = similar(Rsamp)

	for (coordinate, index) in enumerate(CartesianIndices(tensor_size)), trial in 1:ntrials
		u = q[trial]
		Rsamp[index, trial] = coordinate * (10 + u)
		Lsamp[index, trial] = coordinate * (20 + 2u)
		Gsamp[index, trial] = coordinate * (30 - 3u)
		Csamp[index, trial] = coordinate * (40 + 4u)
	end

	f = [50.0, 150.0]
	lp = UQ._joint_line_parameters(PhaseDomain, Rsamp, Lsamp, Gsamp, Csamp, f)
	ω = reshape(2π .* f, 1, 1, :)
	Rmeas = real.(lp.Z.values)
	Lmeas = imag.(lp.Z.values) ./ ω
	Gmeas = real.(lp.Y.values)
	Cmeas = imag.(lp.Y.values) ./ ω

	X = vcat(
		reshape(Rsamp, :, ntrials),
		reshape(Lsamp, :, ntrials),
		reshape(Gsamp, :, ntrials),
		reshape(Csamp, :, ntrials),
	)
	joint = vcat(vec(Rmeas), vec(Lmeas), vec(Gmeas), vec(Cmeas))
	μ = vec(mean(X; dims = 2))
	centered = X .- μ
	empirical_covariance = centered * transpose(centered) / (ntrials - 1)

	@test value.(joint) ≈ μ rtol = 4eps()
	@test uncertainty.(joint) ≈ vec(std(X; dims = 2)) rtol = 16eps()
	@test Measurements.cov(joint) ≈ empirical_covariance rtol = 64eps()

	# Negative covariance and every cross-boundary covariance use the same
	# primitive set: Z/Y, entries, complex components, and frequencies.
	@test Measurements.cov(Rmeas[1, 1, 1], Gmeas[1, 1, 1]) < 0
	@test !iszero(Measurements.cov(Rmeas[1, 1, 1], Rmeas[2, 1, 1]))
	@test !iszero(Measurements.cov(real(lp.Z[1, 1, 1]), imag(lp.Z[1, 1, 1])))
	@test !iszero(Measurements.cov(real(lp.Z[1, 1, 1]), real(lp.Z[1, 1, 2])))

	# R[2,1,1] is exactly twice R[1,1,1] in every trial.
	linear_residual = Rmeas[2, 1, 1] - 2Rmeas[1, 1, 1]
	@test iszero(value(linear_residual))
	@test iszero(uncertainty(linear_residual))

	# A single trial is represented explicitly with zero uncertainty.
	one_trial = UQ._joint_line_parameters(
		PhaseDomain,
		fill(1.0, 1, 1, 1, 1),
		fill(2.0, 1, 1, 1, 1),
		fill(3.0, 1, 1, 1, 1),
		fill(4.0, 1, 1, 1, 1),
		[50.0],
	)
	@test iszero(uncertainty(real(one_trial.Z[1, 1, 1])))
	@test iszero(uncertainty(imag(one_trial.Z[1, 1, 1])))
	@test iszero(uncertainty(real(one_trial.Y[1, 1, 1])))
	@test iszero(uncertainty(imag(one_trial.Y[1, 1, 1])))
end
