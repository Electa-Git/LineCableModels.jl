"""
	mc(cbs::CableBuilderSpec; trials=nothing, distribution=:normal, seed=nothing,
	   trial_sampler=nothing, conf=0.95, tol=0.02, print_step=1000,
	   return_samples=false, return_pdf=false, nbins=nothing)

Propagate cable-design uncertainty by Monte Carlo sampling.

# Arguments

- `cbs`: Cable-builder specification containing the uncertain primitives.
- `trials`: Number of trials, or `nothing` to size the run with the DKW bound.
- `distribution`: Primitive sampling law, either `:normal` or `:uniform`.
- `seed`: Random-number seed, or `nothing` to preserve the active RNG state.
- `trial_sampler`: Optional callback called as
  `trial_sampler(deterministic_spec, trial_index, distribution)`. It must return
  a sampled cable design accepted by `DataFrame(design, :baseparams)`.
- `conf`: Confidence level used for DKW sizing and mean confidence intervals.
- `tol`: DKW absolute CDF tolerance when `trials === nothing`.
- `print_step`: Number of trials between progress messages.
- `return_samples`: Retain the scalar R, L, and C trial arrays when `true`.
- `return_pdf`: Construct histogram-based PDFs when `true`.
- `nbins`: Histogram bin count, or `nothing` for automatic selection.

# Returns

- A `CableDesignMC` containing R, L, and C statistics and Measurements.jl values.

# Notes

`trial_sampler` permits callers to impose shared or correlated primitive draws.
The callback receives one-based trial indices and is invoked exactly once per
trial. Its result bypasses the package's default independent primitive sampler.
"""
function mc(cbs::CableBuilderSpec;
	trials::Union{Int, Nothing} = nothing,
	distribution::Symbol        = :normal,
	seed::Union{Int, Nothing}   = nothing,
	trial_sampler::Union{Nothing, Function} = nothing,
	conf::Float64               = 0.95,
	tol::Float64                = 0.02,     # used only if trials === nothing (DKW sizing)
	print_step::Int             = 1000,
	return_samples::Bool        = false,
	return_pdf::Bool            = false,
	nbins::Union{Int, Nothing}  = nothing,
)
	seed !== nothing && Random.seed!(seed)
	z = quantile(Distributions.Normal(), 0.5 + conf/2)

	# 3 scalar observables: R, L, C
	M = 3
	ntrials = if trials === nothing
		α = 1 - conf
		ceil(Int, log(2 * M / α) / (2 * tol^2))
	else
		trials
	end

	if trials === nothing
		@info "mc: estimate number of trials using DKW inequality" scalars = M conf = conf tol =
			tol trials = ntrials
	end

	@info "mc: starting draws" draws = ntrials conf = conf tol = tol distribution =
		(distribution === :uniform ? "Uniform(μ ± √3·σ)" : "Normal(μ, σ)")

	# Base float type — enforced upstream
	T = BASE_FLOAT

	μR = Vector{T}(undef, ntrials)
	μL = Vector{T}(undef, ntrials)
	μC = Vector{T}(undef, ntrials)
	cbs_det = determinize(cbs)

	@inline function _draw!(i::Int)
		des = if trial_sampler === nothing
			sample(cbs_det; distribution = distribution)
		else
			trial_sampler(cbs_det, i, distribution)
		end
		params = DataFrame(des, :baseparams).computed  # invariant ordering: R, L, C
		r      = params[1]
		l      = params[2]
		c      = params[3]
		@inbounds begin
			μR[i] = T(r) / 1e3 # ohm/km to ohm/m
			μL[i] = T(l) / 1e6 # mH/km to H/m
			μC[i] = T(c) / 1e9 # μF/km to F/m
		end
		return nothing
	end

	for i in 1:ntrials
		_draw!(i)
		(i % print_step == 0) && @info "mc: progress" done = i
	end
	@info "mc: done" total = ntrials

	# stats kernel (scalar real vector → NamedTuple)
	_stats = function (arr::AbstractVector{<:Real})
		m  = mean(arr)
		s  = std(arr)
		N  = length(arr)
		ci = z * s / sqrt(N)
		return (mean    = m, std     = s, min     = minimum(arr),
			q05     = quantile(arr, 0.05),
			q50     = quantile(arr, 0.50),
			q95     = quantile(arr, 0.95),
			max     = maximum(arr),
			n       = N,
			ci_half = ci,
			ci_rel  = ci / max(abs(m), eps()))
	end

	sR = _stats(μR)
	sL = _stats(μL)
	sC = _stats(μC)

	meas = Measurement{T}[
		measurement(sR.mean, sR.std),
		measurement(sL.mean, sL.std),
		measurement(sC.mean, sC.std),
	]

	# PDFs (optional)
	pdf_nt = nothing
	if return_pdf
		pdfR = _pdf_from_hist(μR; nbins = nbins)
		pdfL = _pdf_from_hist(μL; nbins = nbins)
		pdfC = _pdf_from_hist(μC; nbins = nbins)
		pdf_nt = (R = pdfR, L = pdfL, C = pdfC)
	end

	# Samples as NamedTuple of vectors (R,L,C) or nothing
	samples_nt = return_samples ? (R = μR, L = μL, C = μC) : nothing

	return CableDesignMC{T}(
		(R = sR, L = sL, C = sC),
		pdf_nt,
		samples_nt,
		meas,
	)
end


"""
	mc(sbs::SystemBuilderSpec, F::EMTFormulation; trials=nothing,
	   distribution=:normal, seed=nothing, trial_sampler=nothing, conf=0.95,
	   tol=0.02, print_step=1000, return_samples=false, return_pdf=false,
	   per_length=true, nbins=nothing)

Propagate system uncertainty through a frequency-domain EMT formulation.

# Arguments

- `sbs`: Cable-system specification containing the uncertain primitives and
  frequency vector \\[Hz\\].
- `F`: EMT formulation used to compute the series-impedance and shunt-admittance
  matrices.
- `trials`: Number of trials, or `nothing` to size the run with the DKW bound.
- `distribution`: Primitive sampling law, either `:normal` or `:uniform`.
- `seed`: Random-number seed, or `nothing` to preserve the active RNG state.
- `trial_sampler`: Optional callback called as
  `trial_sampler(deterministic_spec, trial_index, distribution)`. It must return
  a sampled system specification accepted by `compute!`.
- `conf`: Confidence level used for DKW sizing and mean confidence intervals.
- `tol`: DKW absolute CDF tolerance when `trials === nothing`.
- `print_step`: Number of trials between progress messages.
- `return_samples`: Retain R, L, C, and G trial tensors when `true`.
- `return_pdf`: Construct histogram-based PDFs when `true`.
- `per_length`: Return per-unit-length quantities when `true`; otherwise scale
  impedance and admittance by the physical line length \\[m\\].
- `nbins`: Histogram bin count, or `nothing` for automatic selection.

# Returns

- A `LineParametersMC` containing entrywise R, L, C, and G statistics, optional
  samples and PDFs, and a Measurements.jl-valued `LineParameters` object.

# Notes

The first sampled problem is solved before tensor allocation, so matrix size is
taken from the actual solver result after bundling and reduction. When supplied,
`trial_sampler` is invoked exactly once for every one-based trial index and its
result bypasses the package's default independent primitive sampler.
"""
function mc(
	sbs::SystemBuilderSpec,
	F::EMTFormulation;
	trials::Union{Int, Nothing} = nothing,
	distribution::Symbol        = :normal,      # :uniform => Uniform(μ ± √3·σ), :normal => Normal(μ,σ)
	seed::Union{Int, Nothing}   = nothing,
	trial_sampler::Union{Nothing, Function} = nothing,
	conf::Float64               = 0.95,
	tol::Float64                = 0.02,
	print_step::Int             = 1000,
	return_samples::Bool        = false,         # returns Vector{LineParameters} (one per trial)
	return_pdf::Bool            = false,         # hist-based LineParametersPDF per R/L/C/G & freq
	per_length::Bool            = true,         # scale results per length
	nbins::Union{Int, Nothing}  = nothing,
)

	seed !== nothing && Random.seed!(seed)
	z = quantile(Distributions.Normal(), 0.5 + conf/2)

	fvec  = sbs.frequencies
	nfreq = length(fvec)

	trials === nothing || trials > 0 ||
		throw(ArgumentError("mc: trials must be greater than zero"))

	# Materialize and solve the first draw before allocating result tensors. The
	# output dimension depends on the complete solver reduction policy (bundling,
	# Kron reduction, retained grounded terminals, and modal transformation), so
	# it cannot be inferred reliably from position dictionaries alone.
	sys_det = determinize(sbs)
	first_problem = if trial_sampler === nothing
		sample(sys_det; distribution = distribution)
	else
		trial_sampler(sys_det, 1, distribution)
	end
	first_ws, first_lp = compute!(first_problem, F)
	nph = size(first_lp.Z, 1)
	size(first_lp.Z) == size(first_lp.Y) ||
		throw(DimensionMismatch("mc: first-trial Z and Y dimensions differ"))
	size(first_lp.Z, 3) == nfreq ||
		throw(DimensionMismatch("mc: first-trial frequency dimension differs from the specification"))

	# Total scalar observables under DKW: Z & Y, Real & Imag, upper-triangular (incl. diag) per freq
	M = 2 * nph * (nph + 1) * nfreq

	ntrials = if trials === nothing
		α = 1 - conf
		ceil(Int, log(2 * M / α) / (2 * tol^2))
	else
		trials
	end

	if trials === nothing
		@info "mc: estimate number of trials using DKW inequality" scalars = M conf = conf tol =
			tol trials = ntrials
	end

	@info "mc[Z,Y]: starting" draws = ntrials conf = conf tol = tol distribution =
		(distribution === :uniform ? "Uniform(μ ± √3·σ)" : "Normal(μ, σ)")

	# Stats kernel on reals
	_stats = function (arr::AbstractVector{<:Real})
		m  = mean(arr)
		s  = std(arr)
		N  = length(arr)
		ci = z * s / sqrt(N)
		return (mean = m, std = s, min = minimum(arr),
			q05 = quantile(arr, 0.05), q50 = quantile(arr, 0.50),
			q95 = quantile(arr, 0.95),
			max = maximum(arr), n = N, conf = conf, z = z,
			ci_half = ci, ci_rel = ci / max(abs(m), eps()))
	end

	U = eltype(fvec)

	# Concrete vectors of RLCG samples
	Rsamp = Array{U, 4}(undef, nph, nph, nfreq, ntrials)
	Lsamp = Array{U, 4}(undef, nph, nph, nfreq, ntrials)
	Gsamp = Array{U, 4}(undef, nph, nph, nfreq, ntrials)
	Csamp = Array{U, 4}(undef, nph, nph, nfreq, ntrials)

	# ─────────────────────────────────────────────────────────────────────────
	# Monte Carlo over FULL frequency vector: one LineParameters per trial
	# ─────────────────────────────────────────────────────────────────────────
	Dlp = domain(first_lp)

	for i in 1:ntrials
		if i == 1
			ws, lp = first_ws, first_lp
		else
			prob = if trial_sampler === nothing
				sample(sys_det; distribution = distribution)
			else
				trial_sampler(sys_det, i, distribution)
			end
			ws, lp = compute!(prob, F)     # lp::LineParameters{Tc, Tr}
		end
		domain(lp) === Dlp || throw(
			DomainError(
				domain(lp),
				"mc: inconsistent LineParameters domain across trials",
			),
		)
		size(lp.Z) == (nph, nph, nfreq) || throw(
			DimensionMismatch("mc: LineParameters dimensions changed between trials"),
		)
		size(lp.Y) == (nph, nph, nfreq) || throw(
			DimensionMismatch("mc: LineParameters dimensions changed between trials"),
		)
		if per_length
			Zscaled = lp.Z.values
			Yscaled = lp.Y.values
		else
			Zscaled = lp.Z.values .* ws.line_length
			Yscaled = lp.Y.values .* ws.line_length
		end

		@inbounds for j1 in 1:nph, j2 in 1:nph, k in 1:nfreq
			Zval = Zscaled[j1, j2, k]
			Yval = Yscaled[j1, j2, k]
			fk   = fvec[k]
			ω   = 2π * fk

			Rval = real(Zval)
			Lval = imag(Zval) / ω
			Gval = real(Yval)
			Cval = imag(Yval) / ω

			Rsamp[j1, j2, k, i] = Rval
			Lsamp[j1, j2, k, i] = Lval
			Gsamp[j1, j2, k, i] = Gval
			Csamp[j1, j2, k, i] = Cval
		end

		(i % print_step == 0) && @info "mc[Z,Y]: progress" done = i
	end

	# ─────────────────────────────────────────────────────────────────────────
	# Aggregate statistics per element (j1,j2) and per frequency k
	# ─────────────────────────────────────────────────────────────────────────

	# Measurement-valued Z,Y (nph×nph×nfreq)
	T = Complex{typeof(measurement(zero(U), zero(U)))}   # Measurement{BASE_FLOAT}
	Zmeas = Array{T}(undef, nph, nph, nfreq)
	Ymeas = Array{T}(undef, nph, nph, nfreq)

	# 3D arrays of stats for R,L,C,G
	Rstats = Array{NamedTuple, 3}(undef, nph, nph, nfreq)
	Lstats = Array{NamedTuple, 3}(undef, nph, nph, nfreq)
	Gstats = Array{NamedTuple, 3}(undef, nph, nph, nfreq)
	Cstats = Array{NamedTuple, 3}(undef, nph, nph, nfreq)

	# Optional PDFs: same 3D shape, one distribution per scalar
	Rpdf = return_pdf ? Array{LineParametersPDF{U}, 3}(undef, nph, nph, nfreq) : nothing
	Lpdf = return_pdf ? Array{LineParametersPDF{U}, 3}(undef, nph, nph, nfreq) : nothing
	Gpdf = return_pdf ? Array{LineParametersPDF{U}, 3}(undef, nph, nph, nfreq) : nothing
	Cpdf = return_pdf ? Array{LineParametersPDF{U}, 3}(undef, nph, nph, nfreq) : nothing

	@inbounds for j1 in 1:nph, j2 in 1:nph, k in 1:nfreq
		fk = fvec[k]
		ω = 2π * fk

		rvec = @view Rsamp[j1, j2, k, :]
		lvec = @view Lsamp[j1, j2, k, :]
		gvec = @view Gsamp[j1, j2, k, :]
		cvec = @view Csamp[j1, j2, k, :]

		sR = _stats(rvec)
		sL = _stats(lvec)
		sG = _stats(gvec)
		sC = _stats(cvec)

		Rstats[j1, j2, k] = sR
		Lstats[j1, j2, k] = sL
		Gstats[j1, j2, k] = sG
		Cstats[j1, j2, k] = sC

		# Z = R + j ω L, Y = G + j ω C
		Zmeas[j1, j2, k] =
			measurement(sR.mean, sR.std) +
			1im * measurement(ω * sL.mean, ω * sL.std)

		Ymeas[j1, j2, k] =
			measurement(sG.mean, sG.std) +
			1im * measurement(ω * sC.mean, ω * sC.std)

		# Optional PDFs per (j1,j2,k)
		if return_pdf
			Rpdf[j1, j2, k] = _pdf_from_hist(rvec; nbins = nbins)
			Lpdf[j1, j2, k] = _pdf_from_hist(lvec; nbins = nbins)
			Gpdf[j1, j2, k] = _pdf_from_hist(gvec; nbins = nbins)
			Cpdf[j1, j2, k] = _pdf_from_hist(cvec; nbins = nbins)
		end
	end

	# Frequency-dependent LineParameters using Measurements.jl types
	LP_meas = LineParameters(Dlp, Zmeas, Ymeas, fvec)

	@info "mc[Z,Y]: done" total = ntrials nfreq = nfreq

	stats_nt = (R = Rstats, L = Lstats, C = Cstats, G = Gstats)
	pdf_nt   = return_pdf ? (R = Rpdf, L = Lpdf, C = Cpdf, G = Gpdf) : nothing

	samples_nt = return_samples ? (R = Rsamp, L = Lsamp, C = Csamp, G = Gsamp) : nothing

	return LineParametersMC(fvec, stats_nt, pdf_nt, samples_nt, LP_meas)

end
