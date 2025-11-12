# Piecewise-constant PDF 
struct PCPDF
	edges::Vector{Float64}   # length B+1
	dens::Vector{Float64}    # length B (area ≈ 1)
end

# Build a piecewise-constant PDF from samples
function _pdf_from_hist(x::AbstractVector{<:Real}; nbins::Int = 50)
	n = length(x)
	n == 0 && error("Empty sample set.")
	h = fit(Histogram, float.(x); nbins = nbins, closed = :left)
	edges = collect(h.edges[1])
	widths = diff(edges)
	dens = h.weights ./ (n .* widths)  # counts / (n * binwidth)
	return PCPDF(edges, dens)
end

# Density at x0 
@inline function (hp::PCPDF)(x0::Real)
	i = searchsortedlast(hp.edges, float(x0))
	(i < 1 || i >= length(hp.edges)) && return 0.0
	return hp.dens[i]
end

function mc(cbs::CableBuilderSpec;
	trials::Union{Int, Nothing} = nothing,
	distribution::Symbol = :uniform,
	seed::Union{Int, Nothing} = nothing,
	conf::Float64 = 0.95,
	dkw_eps::Float64 = 0.02,     # used only if trials === nothing (DKW sizing)
	batch::Int = 1000,
	return_samples::Bool = false,
	return_pdf::Bool = false,
	nbins::Int = 80,
)
	seed !== nothing && Random.seed!(seed)
	z = quantile(Distributions.Normal(), 0.5 + conf/2)

	n = if trials === nothing
		α = 1 - conf
		ceil(Int, log(2/α) / (2 * dkw_eps^2))
	else
		trials
	end

	μR = Vector{Float64}(undef, n)
	μL = Vector{Float64}(undef, n)
	μC = Vector{Float64}(undef, n)

	@inline function _draw!(i::Int)
		des    = sample(determinize(cbs); distribution = distribution) # cbs is transformed to deterministic ranges
		params = DataFrame(des, :baseparams).computed  # invariant
		r      = params[1];
		l      = params[2];
		c      = params[3]
		@inbounds begin
			μR[i] = float(r)
			μL[i] = float(l)
			μC[i] = float(c)
		end
		return nothing
	end

	@info "mc: starting draws" draws=n conf=conf dkw_eps=dkw_eps distribution =
		distribution === :uniform ? "Uniform(μ ± √3·σ)" : "Normal(μ, σ)"
	for i in 1:n
		_draw!(i)
		(i % batch == 0) && @info "mc: progress" done=i
	end
	@info "mc: done" total=n

	stats = function (arr::AbstractVector{<:Real})
		m  = mean(arr);
		s  = std(arr);
		N  = length(arr)
		ci = z * s / sqrt(N)
		(mean = m, std = s, min = minimum(arr), q05 = quantile(arr, 0.05),
			q50 = quantile(arr, 0.50), q95 = quantile(arr, 0.95), max = maximum(arr),
			n = N, ci_half = ci, ci_rel = ci / max(abs(m), eps()))
	end
	sR = stats(μR);
	sL = stats(μL);
	sC = stats(μC)

	summary = DataFrame(
		variable = ["R", "L", "C"],
		mean     = [sR.mean, sL.mean, sC.mean],
		std      = [sR.std, sL.std, sC.std],
		min      = [sR.min, sL.min, sC.min],
		q05      = [sR.q05, sL.q05, sC.q05],
		q50      = [sR.q50, sL.q50, sC.q50],
		q95      = [sR.q95, sL.q95, sC.q95],
		max      = [sR.max, sL.max, sC.max],
		n        = fill(n, 3),
		conf     = fill(conf, 3),
		z        = fill(z, 3),
		ci_half  = [sR.ci_half, sL.ci_half, sC.ci_half],
		ci_rel   = [sR.ci_rel, sL.ci_rel, sC.ci_rel],
	)

	if return_pdf || return_samples
		ret = (summary = summary,)
		if return_samples
			samples = DataFrame(μ_R = μR, μ_L = μL, μ_C = μC)
			ret = merge(ret, (samples = samples,))
		end
		if return_pdf
			pdfR = _pdf_from_hist(μR; nbins = nbins)
			pdfL = _pdf_from_hist(μL; nbins = nbins)
			pdfC = _pdf_from_hist(μC; nbins = nbins)
			ret = merge(ret, (pdf = (R = pdfR, L = pdfL, C = pdfC),))
		end
		return ret
	else
		return summary
	end
end

# # ─────────────────────────────────────────────────────────────────────────────
# # Monte Carlo for frequency scans (matrix Z, Y) — mirrors mc(cbs::CableBuilderSpec)
# # ─────────────────────────────────────────────────────────────────────────────
# function mc(
# 	sbs::SystemBuilderSpec,
# 	F::AbstractFormulationSet;
# 	trials::Union{Int, Nothing} = nothing,
# 	distribution::Symbol        = :uniform,      # :uniform => Uniform(μ ± √3·σ), :normal => Normal(μ,σ)
# 	seed::Union{Int, Nothing}   = nothing,
# 	conf::Float64               = 0.95,
# 	dkw_eps::Float64            = 0.02,
# 	batch::Int                  = 1000,
# 	return_samples::Bool        = false,         # heavy for big n×n×trials → off by default
# 	return_pdf::Bool            = false,         # hist-based PCPDF per element & part (Re/Im)
# 	nbins::Int                  = 80,
# )
# 	seed !== nothing && Random.seed!(seed)
# 	zcrit = quantile(Distributions.Normal(), 0.5 + conf/2)

# 	ntrials = trials === nothing ? ceil(Int, log(2/(1-conf)) / (2*dkw_eps^2)) : trials
# 	fvec = sbs.frequencies  # keep sbs frequencies verbatim

# 	# Local helper: build a single-frequency SystemBuilderSpec 
# 	_single_freq = function (s::SystemBuilderSpec, f1::Number)
# 		return typeof(s)(
# 			s.system_id,
# 			s.builder,
# 			s.positions;
# 			length      = s.length,
# 			temperature = s.temperature,
# 			earth       = s.earth,
# 			f           = [f1],
# 		)
# 	end

# 	# Reuse  mc stats kernel (mean/std/ci etc.) but on real vectors
# 	_stats = function (arr::AbstractVector{<:Real})
# 		m  = mean(arr)
# 		s  = std(arr)
# 		N  = length(arr)
# 		ci = zcrit * s / sqrt(N)
# 		return (mean = m, std = s, min = minimum(arr),
# 			q05 = quantile(arr, 0.05), q50 = quantile(arr, 0.50),
# 			q95 = quantile(arr, 0.95),
# 			max = maximum(arr), n = N, ci_half = ci, ci_rel = ci / max(abs(m), eps()))
# 	end

# 	# Labels like "Z(1,1)"
# 	_labels =
# 		(n::Int, prefix::AbstractString) ->
# 			[string(prefix, "(", i, ",", j, ")") for i ∈ 1:n, j ∈ 1:n] |> vec

# 	# Optional PDF builder (same semantics as your scalar mc)
# 	_mkpdf = (x::AbstractVector{<:Real}) -> _pdf_from_hist(x; nbins = nbins)  # you already have _pdf_from_hist

# 	out = Vector{NamedTuple}(undef, length(fvec))  # one (Z=..., Y=..., f=...) per frequency

# 	@info "mc[Z,Y]: starting" draws=ntrials conf=conf dkw_eps=dkw_eps distribution =
# 		distribution === :uniform ? "Uniform(μ ± √3·σ)" : "Normal(μ, σ)"

# 	# Frequency loop — each k is a fully independent MC on Z(·,·) and Y(·,·)
# 	for (k, fk) in pairs(fvec)
# 		# Build a single-frequency, already-determinized spec for sizing (no compute!)
# 		s1 = _single_freq(sbs, fk)
# 		ss = collapse_sbs(determinize(s1); distribution = distribution)  # your existing helper
# 		system = ss.system  # same object used later to spawn problems

# 		# Phase count from the system model (no EM solve)
# 		nph = sum(length(cable.design_data.components) for cable in system.cables)

# 		# Preallocate storage for Re/Im samples (nph×nph entries, ntrials samples)
# 		Zre = [Vector{Float64}(undef, ntrials) for _ in 1:(nph*nph)]
# 		Zim = [Vector{Float64}(undef, ntrials) for _ in 1:(nph*nph)]
# 		Yre = [Vector{Float64}(undef, ntrials) for _ in 1:(nph*nph)]
# 		Yim = [Vector{Float64}(undef, ntrials) for _ in 1:(nph*nph)]

# 		# MC draws
# 		s1d = determinize(s1)

# 		for i in 1:ntrials
# 			prob = sample(s1d; distribution = distribution)
# 			_, p = compute!(prob, F)
# 			Z, Y = p.Z[:, :, 1], p.Y[:, :, 1]
# 			@inbounds for j2 in 1:nph, j1 in 1:nph
# 				idx = (j1-1)*nph + j2
# 				Zre[idx][i] = real(Z[j1, j2]);
# 				Zim[idx][i] = imag(Z[j1, j2])
# 				Yre[idx][i] = real(Y[j1, j2]);
# 				Yim[idx][i] = imag(Y[j1, j2])
# 			end
# 			(i % batch == 0) && @info "mc[Z,Y]: progress" f=fk done=i
# 		end

# 		# Aggregate stats per element → DataFrames
# 		zlabels = _labels(nph, "Z");
# 		ylabels = _labels(nph, "Y")
# 		ZreS = map(_stats, Zre);
# 		ZimS = map(_stats, Zim)
# 		YreS = map(_stats, Yre);
# 		YimS = map(_stats, Yim)

# 		dfZ = DataFrame(
# 			variable = zlabels,
# 			mean_re  = getfield.(ZreS, :mean), mean_im  = getfield.(ZimS, :mean),
# 			std_re   = getfield.(ZreS, :std), std_im   = getfield.(ZimS, :std),
# 			min_re   = getfield.(ZreS, :min), min_im   = getfield.(ZimS, :min),
# 			q05_re   = getfield.(ZreS, :q05), q05_im   = getfield.(ZimS, :q05),
# 			q50_re   = getfield.(ZreS, :q50), q50_im   = getfield.(ZimS, :q50),
# 			q95_re   = getfield.(ZreS, :q95), q95_im   = getfield.(ZimS, :q95),
# 			max_re   = getfield.(ZreS, :max), max_im   = getfield.(ZimS, :max),
# 			n        = fill(ntrials, nph*nph),
# 			conf     = fill(conf, nph*nph),
# 			z        = fill(zcrit, nph*nph),
# 			ci_re    = getfield.(ZreS, :ci_half),
# 			ci_im    = getfield.(ZimS, :ci_half),
# 			cirel_re = getfield.(ZreS, :ci_rel),
# 			cirel_im = getfield.(ZimS, :ci_rel),
# 		)

# 		dfY = DataFrame(
# 			variable = ylabels,
# 			mean_re  = getfield.(YreS, :mean), mean_im  = getfield.(YimS, :mean),
# 			std_re   = getfield.(YreS, :std), std_im   = getfield.(YimS, :std),
# 			min_re   = getfield.(YreS, :min), min_im   = getfield.(YimS, :min),
# 			q05_re   = getfield.(YreS, :q05), q05_im   = getfield.(YimS, :q05),
# 			q50_re   = getfield.(YreS, :q50), q50_im   = getfield.(YimS, :q50),
# 			q95_re   = getfield.(YreS, :q95), q95_im   = getfield.(YimS, :q95),
# 			max_re   = getfield.(YreS, :max), max_im   = getfield.(YimS, :max),
# 			n        = fill(ntrials, nph*nph),
# 			conf     = fill(conf, nph*nph),
# 			z        = fill(zcrit, nph*nph),
# 			ci_re    = getfield.(YreS, :ci_half),
# 			ci_im    = getfield.(YimS, :ci_half),
# 			cirel_re = getfield.(YreS, :ci_rel),
# 			cirel_im = getfield.(YimS, :ci_rel),
# 		)

# 		# Optional payloads (off by default to avoid memory carnage)
# 		ret = (f = fk, Z = dfZ, Y = dfY)
# 		if return_samples || return_pdf
# 			samples =
# 				return_samples ? (
# 					Z_re = Zre, Z_im = Zim, Y_re = Yre, Y_im = Yim,
# 				) : nothing
# 			pdfs =
# 				return_pdf ?
# 				(
# 					Z_re = map(_mkpdf, Zre), Z_im = map(_mkpdf, Zim),
# 					Y_re = map(_mkpdf, Yre), Y_im = map(_mkpdf, Yim),
# 				) : nothing
# 			ret = merge(ret, (samples = samples, pdf = pdfs))
# 		end

# 		out[k] = ret
# 	end

# 	@info "mc[Z,Y]: done" total=ntrials nf=length(fvec)
# 	return out
# end
