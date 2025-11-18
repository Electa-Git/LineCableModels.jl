import DataFrames: DataFrame

function DataFrame(res::LineParametersMCSummary)
	nph, _, nfreq = size(res.stats.R)
	dfs = Array{DataFrame, 3}(undef, nph, nph, nfreq)
	@inbounds for i in 1:nph, j in 1:nph, k in 1:nfreq
		r = res.stats.R[i, j, k]
		l = res.stats.L[i, j, k]
		c = res.stats.C[i, j, k]
		g = res.stats.G[i, j, k]
		dfs[i, j, k] = DataFrame(
			quantity = ["R", "L", "C", "G"],
			mean     = [r.mean, l.mean, c.mean, g.mean],
			std      = [r.std, l.std, c.std, g.std],
			min      = [r.min, l.min, c.min, g.min],
			q05      = [r.q05, l.q05, c.q05, g.q05],
			q50      = [r.q50, l.q50, c.q50, g.q50],
			q95      = [r.q95, l.q95, c.q95, g.q95],
			max      = [r.max, l.max, c.max, g.max],
			n        = [r.n, l.n, c.n, g.n],
			conf     = [r.conf, l.conf, c.conf, g.conf],
			z        = [r.z, l.z, c.z, g.z],
			ci_half  = [r.ci_half, l.ci_half, c.ci_half, g.ci_half],
			ci_rel   = [r.ci_rel, l.ci_rel, c.ci_rel, g.ci_rel],
		)
	end
	return dfs
end

function DataFrame(res::CableDesignMCSummary)
	sR, sL, sC = res.stats.R, res.stats.L, res.stats.C
	DataFrame(
		variable = ["R", "L", "C"],
		mean     = [sR.mean, sL.mean, sC.mean],
		std      = [sR.std, sL.std, sC.std],
		min      = [sR.min, sL.min, sC.min],
		q05      = [sR.q05, sL.q05, sC.q05],
		q50      = [sR.q50, sL.q50, sC.q50],
		q95      = [sR.q95, sL.q95, sC.q95],
		max      = [sR.max, sL.max, sC.max],
		ntrials  = fill(sR.n, 3),
		ci_half  = [sR.ci_half, sL.ci_half, sC.ci_half],
		ci_rel   = [sR.ci_rel, sL.ci_rel, sC.ci_rel],
	)
end
