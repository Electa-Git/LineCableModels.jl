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

	_stats = function (arr::AbstractVector{<:Real})
		m  = mean(arr);
		s  = std(arr);
		N  = length(arr)
		ci = z * s / sqrt(N)
		(mean = m, std = s, min = minimum(arr), q05 = quantile(arr, 0.05),
			q50 = quantile(arr, 0.50), q95 = quantile(arr, 0.95), max = maximum(arr),
			n = N, ci_half = ci, ci_rel = ci / max(abs(m), eps()))
	end

	sR = _stats(μR)
	sL = _stats(μL)
	sC = _stats(μC)

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

	meas = [measurement(sR.mean, sR.std), measurement(sL.mean, sL.std),
		measurement(sC.mean, sC.std)]

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
		return ret, meas
	else
		return summary, meas
	end
end

# ─────────────────────────────────────────────────────────────────────────────
# Monte Carlo for frequency scans (matrix Z, Y) — mirrors mc(cbs::CableBuilderSpec)
# ─────────────────────────────────────────────────────────────────────────────
function mc(
	sbs::SystemBuilderSpec,
	F::EMTFormulation;
	trials::Union{Int, Nothing} = nothing,
	distribution::Symbol = :uniform,      # :uniform => Uniform(μ ± √3·σ), :normal => Normal(μ,σ)
	seed::Union{Int, Nothing} = nothing,
	conf::Float64 = 0.95,
	dkw_eps::Float64 = 0.02,
	batch::Int = 1000,
	return_samples::Bool = false,         # heavy for big n×n×trials → off by default
	return_pdf::Bool = false,         # hist-based PCPDF per element & part (Real/Imag)
	nbins::Int = 80,
)

	seed !== nothing && Random.seed!(seed)
	z = quantile(Distributions.Normal(), 0.5 + conf/2)

	ntrials = trials === nothing ? ceil(Int, log(2/(1-conf)) / (2*dkw_eps^2)) : trials
	fvec = sbs.frequencies  # keep sbs frequencies verbatim

	# Local helper: build a single-frequency SystemBuilderSpec
	_single_freq = function (s::SystemBuilderSpec, f1::Number)
		return typeof(s)(
			s.system_id,
			s.builder,
			s.positions;
			length      = s.length,
			temperature = s.temperature,
			earth       = s.earth,
			f           = [f1],
		)
	end

	# Reuse mc stats kernel (mean/std/ci etc.) but on real vectors
	_stats = function (arr::AbstractVector{<:Real})
		m  = mean(arr)
		s  = std(arr)
		N  = length(arr)
		ci = z * s / sqrt(N)
		return (mean = m, std = s, min = minimum(arr),
			q05 = quantile(arr, 0.05), q50 = quantile(arr, 0.50),
			q95 = quantile(arr, 0.95),
			max = maximum(arr), n = N, ci_half = ci, ci_rel = ci / max(abs(m), eps()))
	end

	# Optional PDF builders
	_mkpdf_real = (x::AbstractVector{<:Real}) -> _pdf_from_hist(x; nbins = nbins)
	_mkpdf_pair =
		(x::AbstractVector{<:Complex}) ->
			(real = _mkpdf_real(real.(x)), imag = _mkpdf_real(imag.(x)))

	# complex Measurement eltype without hard-coding param
	nph = count(!=(0), vcat(values.(getfield.(sbs.positions, :conn))...))
	mT = typeof(measurement(0.0, 0.0))
	Zmeas = Array{Complex{mT}}(undef, nph, nph, length(fvec))
	Ymeas = Array{Complex{mT}}(undef, nph, nph, length(fvec))

	out = Vector{NamedTuple}(undef, length(fvec))  # one (Z=..., Y=..., f=...) per frequency

	@info "mc[Z,Y]: starting" draws=ntrials conf=conf dkw_eps=dkw_eps distribution =
		distribution === :uniform ? "Uniform(μ ± √3·σ)" : "Normal(μ, σ)"

	# Frequency loop — each k is a fully independent MC on Z(·,·) and Y(·,·)
	for (k, fk) in pairs(fvec)
		system = determinize(_single_freq(sbs, fk))

		# Preallocate complex sample cubes (nph×nph×ntrials)
		Zs = Array{ComplexF64}(undef, nph, nph, ntrials)
		Ys = Array{ComplexF64}(undef, nph, nph, ntrials)

		# MC draws
		for i in 1:ntrials
			prob = sample(system; distribution = distribution)
			_, p = compute!(prob, F)
			@inbounds begin
				Z = p.Z[:, :, 1];
				Y = p.Y[:, :, 1]
				Zs[:, :, i] = Z
				Ys[:, :, i] = Y
			end
			(i % batch == 0) && @info "mc[Z,Y]: progress" f=fk done=i
		end

		# Aggregate stats per element → matrices of DataFrames with rows Real/Imag
		dfZ = Matrix{DataFrame}(undef, nph, nph)
		dfY = Matrix{DataFrame}(undef, nph, nph)

		@inbounds for j1 in 1:nph, j2 in 1:nph
			# Views to avoid extra copies
			zvec = @view Zs[j1, j2, :]
			yvec = @view Ys[j1, j2, :]

			# Component-wise stats
			zr = _stats(real.(zvec));
			zi = _stats(imag.(zvec))
			yr = _stats(real.(yvec));
			yi = _stats(imag.(yvec))

			# Build two-row DataFrames (Real/Imag) with consistent schema
			dfZ[j1, j2] = DataFrame(
				part    = ["Real", "Imaginary"],
				mean    = [zr.mean, zi.mean],
				std     = [zr.std, zi.std],
				min     = [zr.min, zi.min],
				q05     = [zr.q05, zi.q05],
				q50     = [zr.q50, zi.q50],
				q95     = [zr.q95, zi.q95],
				max     = [zr.max, zi.max],
				n       = [zr.n, zi.n],
				conf    = [conf, conf],
				z       = [z, z],
				ci_half = [zr.ci_half, zi.ci_half],
				ci_rel  = [zr.ci_rel, zi.ci_rel],
			)

			dfY[j1, j2] = DataFrame(
				part    = ["Real", "Imaginary"],
				mean    = [yr.mean, yi.mean],
				std     = [yr.std, yi.std],
				min     = [yr.min, yi.min],
				q05     = [yr.q05, yi.q05],
				q50     = [yr.q50, yi.q50],
				q95     = [yr.q95, yi.q95],
				max     = [yr.max, yi.max],
				n       = [yr.n, yi.n],
				conf    = [conf, conf],
				z       = [z, z],
				ci_half = [yr.ci_half, yi.ci_half],
				ci_rel  = [yr.ci_rel, yi.ci_rel],
			)
			# Complex measurement per entry at this frequency index k
			Zmeas[j1, j2, k] =
				measurement(zr.mean, zr.std) + 1im * measurement(zi.mean, zi.std)
			Ymeas[j1, j2, k] =
				measurement(yr.mean, yr.std) + 1im * measurement(yi.mean, yi.std)
		end

		# Optional payloads
		ret = (f = fk, Z = dfZ, Y = dfY)

		if return_samples || return_pdf
			samples = return_samples ? (Z = Zs, Y = Ys) : nothing

			pdfZ = nothing;
			pdfY = nothing
			if return_pdf
				pdfZ = Matrix{NamedTuple}(undef, nph, nph)
				pdfY = Matrix{NamedTuple}(undef, nph, nph)
				@inbounds for j1 in 1:nph, j2 in 1:nph
					pdfZ[j1, j2] = _mkpdf_pair(@view Zs[j1, j2, :])
					pdfY[j1, j2] = _mkpdf_pair(@view Ys[j1, j2, :])
				end
			end

			ret = merge(ret, (samples = samples, pdf = (Z = pdfZ, Y = pdfY)))
		end

		out[k] = ret
	end

	LP = LineParameters(Zmeas, Ymeas, fvec)


	@info "mc[Z,Y]: done" total=ntrials nf=length(fvec)

	return out, LP
end