function mc(cbs::CableBuilderSpec;
	trials::Union{Int, Nothing} = nothing,
	distribution::Symbol        = :normal,
	seed::Union{Int, Nothing}   = nothing,
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

	@inline function _draw!(i::Int)
		des    = sample(determinize(cbs); distribution = distribution) # cbs is transformed to deterministic ranges
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

	return CableDesignMCSummary{T}(
		(R = sR, L = sL, C = sC),
		pdf_nt,
		samples_nt,
		meas,
	)
end


function mc(
	sbs::SystemBuilderSpec,
	F::EMTFormulation;
	trials::Union{Int, Nothing} = nothing,
	distribution::Symbol        = :normal,      # :uniform => Uniform(μ ± √3·σ), :normal => Normal(μ,σ)
	seed::Union{Int, Nothing}   = nothing,
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

	# number of physical phases from current mapping
	nph = count(!=(0), vcat(values.(getfield.(sbs.positions, :conn))...))

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

	# Determinize once
	sys_det = determinize(sbs)

	# ─────────────────────────────────────────────────────────────────────────
	# Monte Carlo over FULL frequency vector: one LineParameters per trial
	# ─────────────────────────────────────────────────────────────────────────
	for i in 1:ntrials
		prob = sample(sys_det; distribution = distribution)
		ws, lp = compute!(prob, F)     # lp::LineParameters{Tc, Tr}
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
	LP_meas = LineParameters(Zmeas, Ymeas, fvec)

	@info "mc[Z,Y]: done" total = ntrials nfreq = nfreq

	stats_nt = (R = Rstats, L = Lstats, C = Cstats, G = Gstats)
	pdf_nt   = return_pdf ? (R = Rpdf, L = Lpdf, C = Cpdf, G = Gpdf) : nothing

	samples_nt = return_samples ? (R = Rsamp, L = Lsamp, C = Csamp, G = Gsamp) : nothing

	return LineParametersMCSummary{U}(fvec, stats_nt, pdf_nt, samples_nt, LP_meas)

end
