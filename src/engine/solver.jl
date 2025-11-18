function compute!(
	problem::LineParametersProblem{T},
	formulation::EMTFormulation,
) where {T <: REALSCALAR}

	lvl = levelfrom(formulation.options.common.verbosity)
	sink =
		isnothing(formulation.options.logfile) ?
		ConsoleLogger(stderr, lvl) :
		TeeLogger(ConsoleLogger(stderr, lvl),
			FileLogger(formulation.options.logfile, lvl))
	with_logger(TimestampLogger(sink)) do

		@info "Preallocating arrays"

		ws = init_workspace(problem, formulation)
		nph, nfreq = ws.n_phases, ws.n_frequencies

		# --- full matrices are built per slice (no 3D alloc) ----------------------
		Zbuf = Matrix{Complex{T}}(undef, nph, nph)   # reordered scratch (mutated by merge_bundles!)
		Pbuf = Matrix{Complex{T}}(undef, nph, nph)
		inv_Pbuf = similar(Pbuf) # buffer to hold inv(Pbuf)

		Ztmp = Matrix{Complex{T}}(undef, nph, nph)   # raw slice coming from builders
		Ptmp = Matrix{Complex{T}}(undef, nph, nph)

		# --- index plan (constant across k) ---------------------------------------
		phase_map = ws.phase_map::Vector{Int}
		perm      = reorder_indices(phase_map)
		map_r     = phase_map[perm]                  # reordered map

		# bundle tails mask (same logic as merge_bundles!, but map-only)
		reduced_map = let m = copy(map_r), seen = Set{Int}()
			@inbounds for (i, p) in pairs(map_r)
				if p > 0 && (p in seen)
					m[i]=0
				else
					p>0 && push!(seen, p)
				end
			end
			m
		end

		# decide what Kron shall smite upon
		kron_map = if formulation.options.reduce_bundle
			if formulation.options.kron_reduction
				reduced_map                      # kill tails and keep nonzero labels
			else
				km = copy(reduced_map)           # kill only tails; keep phase-0 explicit
				@inbounds for i in eachindex(km)
					if map_r[i] == 0
						km[i] = -1
					end
				end
				km
			end
		else
			formulation.options.kron_reduction ? map_r : nothing
		end

		nkeep = kron_map === nothing ? nph : count(!=(0), kron_map)
		Zout = Array{Complex{T}, 3}(undef, nkeep, nkeep, nfreq)
		Yout = Array{Complex{T}, 3}(undef, nkeep, nkeep, nfreq)
		Mred = Matrix{Complex{T}}(undef, nkeep, nkeep) # buffer to hold Mred
		inv_Mred = similar(Mred) # buffer to hold inv(Mred)

		# tiny gather helper to avoid per-slice allocs
		@inline function _reorder_into!(dest::AbstractMatrix{Complex{T}},
			src::AbstractMatrix{Complex{T}},
			perm::AbstractVector{Int})
			n = length(perm)
			@inbounds for j in 1:n, i in 1:n
				dest[i, j] = src[perm[i], perm[j]]
			end
			return dest
		end

		# apply temperature correction if needed
		if formulation.options.temperature_correction
			ΔT = ws.temp - T₀
			@. ws.rho_cond *= 1 + ws.alpha_cond * ΔT
		end

		# pre-allocate LU factorization for admittance inversion
		I_nph = Matrix{Complex{T}}(I, nph, nph)      # identity for full size
		I_nkeep = Matrix{Complex{T}}(I, nkeep, nkeep)   # identity for reduced size

		# --- per-frequency pipeline ------------------------------------------------
		@info "Starting line parameters computation"
		for k in 1:nfreq

			compute_impedance_matrix!(Ztmp, ws, k, formulation)
			compute_admittance_matrix!(Ptmp, ws, k, formulation)

			# 1) reorder
			_reorder_into!(Zbuf, Ztmp, perm)
			_reorder_into!(Pbuf, Ptmp, perm)

			# 2) bundle reduction (in-place)
			if formulation.options.reduce_bundle
				merge_bundles!(Zbuf, map_r)
				merge_bundles!(Pbuf, map_r)
			end

			# 3) kron
			if kron_map === nothing
				symtrans!(Zbuf)
				formulation.options.ideal_transposition || line_transpose!(Zbuf)
				@views @inbounds Zout[:, :, k] .= Zbuf

				try
					F = cholesky!(Hermitian(Pbuf))               # in-place factorization
					ldiv!(inv_Pbuf, F, I_nph)                    # inv_Pbuf := P^{-1}
				catch
					F = lu!(Pbuf)                                # overwrite Pbuf with LU
					ldiv!(inv_Pbuf, F, I_nph)                    # inv_Pbuf := P^{-1}
				end
				# inv_Pbuf = pBuf
				inv_Pbuf .*= ws.jω[k]
				symtrans!(inv_Pbuf)
				formulation.options.ideal_transposition || line_transpose!(inv_Pbuf)
				@views @inbounds Yout[:, :, k] .= inv_Pbuf
			else
				kronify!(Zbuf, kron_map, Mred)
				symtrans!(Mred)
				formulation.options.ideal_transposition || line_transpose!(Mred)
				@views @inbounds Zout[:, :, k] .= Mred

				kronify!(Pbuf, kron_map, Mred)
				try
					F = cholesky!(Hermitian(Mred))
					ldiv!(inv_Mred, F, I_nkeep)
				catch
					F = lu!(Mred)
					ldiv!(inv_Mred, F, I_nkeep)
				end
				# inv_Mred = Mred
				inv_Mred .*= ws.jω[k]
				symtrans!(inv_Mred)
				formulation.options.ideal_transposition && line_transpose!(inv_Mred)

				@views @inbounds Yout[:, :, k] .= inv_Mred
			end
		end

		if !isnothing(formulation.modal_transform)
			# apply modal transformation
			_, lp = formulation.modal_transform(LineParameters(Zout, Yout, ws.freq))
		else
			lp = LineParameters(Zout, Yout, ws.freq)
		end

		@info "Line parameters computation completed successfully"
		return ws, lp
	end
end

@inline function stash!(slice_or_nothing, k::Int, src::AbstractMatrix)
	slice_or_nothing === nothing && return nothing
	@views copyto!(slice_or_nothing[:, :, k], src)
	nothing
end

# Builds an Nc×Nc earth matrix using the functors f(h, y, ρ[:,k], ε[:,k], μ[:,k], jω)
@inline function compute_earth_return_matrix!(
	E::AbstractMatrix{Complex{T}},
	cables::AbstractVector{Int},
	ws,
	k::Int,
	functor,                         # formulation.earth_impedance or .earth_admittance
) where {T}
	ρ = @view ws.rho_g[:, k]
	ε = @view ws.eps_g[:, k]
	μ = @view ws.mu_g[:, k]
	jω = ws.jω[k]

	Nc = length(cables)

	@inbounds for cj in 1:Nc
		i = cables[cj]
		for ck in 1:Nc
			j = cables[ck]
			# y: diagonal blocks use cable outer radius; off-diagonals use center distance
			yij = ws.horz_sep[i, j]
			hij = @view ws.vert[[i, j]]
			E[cj, ck] =
				cj == ck ? functor(Val(:self), hij, yij, ρ, ε, μ, jω) :
				functor(Val(:mutual), hij, yij, ρ, ε, μ, jω)
		end
	end

	return nothing
end


function compute_impedance_matrix!(
	Ztmp::AbstractMatrix{Complex{T}},
	ws,
	k::Int,
	formulation,
) where {T <: REALSCALAR}

	@inbounds fill!(Ztmp, zero(Complex{T}))
	@assert length(ws.r_ins_ext) == ws.n_phases "ws.r_ins_ext length mismatch"
	@assert length(ws.mu_ins) == ws.n_phases "ws.mu_ins length mismatch"

	Nc = ws.n_cables
	jω = ws.jω[k]

	cons_in_cable, cables = _get_cable_indices(ws)

	# Earth return impedance (Nc×Nc)
	Zext = Matrix{Complex{T}}(undef, Nc, Nc)
	compute_earth_return_matrix!(Zext, cables, ws, k, formulation.earth_impedance)
	stash!(ws.Zg, k, Zext)

	# ws.Zg[:, :, k] .= Zext # store in workspace for later use

	zinfunctor  = formulation.internal_impedance
	zinsfunctor = formulation.insulation_impedance

	@inbounds for c in 1:Nc
		cons = cons_in_cable[c];
		n = length(cons)

		for p ∈ n:-1:1
			i    = cons[p]
			rin  = ws.r_in[i]
			rex  = ws.r_ext[i]
			ρc  = ws.rho_cond[i]
			μrc = ws.mu_cond[i]

			z_outer  = zinfunctor(:outer, rin, rex, ρc, μrc, jω)
			z_inner  = (p < n) ? zinfunctor(:inner,
			ws.r_in[cons[p+1]],
			ws.r_ext[cons[p+1]],
			ws.rho_cond[cons[p+1]],
			ws.mu_cond[cons[p+1]], jω) : zero(z_outer)
			z_mutual = zinfunctor(:mutual, rin, rex, ρc, μrc, jω)

			# insulation series 
			r_ins_ext = ws.r_ins_ext[i]
			μr_ins = ws.mu_ins[i]
			z_ins = zinsfunctor(rex, r_ins_ext, μr_ins, jω)

			z_loop = z_outer + z_inner + z_ins

			if p > 1
				for a in 1:(p-1), b in 1:(p-1)
					Ztmp[cons[a], cons[b]] += (z_loop - 2*z_mutual)
				end
				for a in 1:(p-1)
					Ztmp[cons[p], cons[a]] += (z_loop - z_mutual)
					Ztmp[cons[a], cons[p]] += (z_loop - z_mutual)
				end
			end
			Ztmp[cons[p], cons[p]] += z_loop
		end

		stash!(ws.Zin, k, Ztmp)

		# self earth-return on intra-cable block
		zgself = Zext[c, c]
		for a in 1:n, b in 1:n
			Ztmp[cons[a], cons[b]] += zgself
		end
	end

	# mutual earth-return off-blocks
	@inbounds for cj in 1:(Nc-1)
		cons_j = cons_in_cable[cj];
		nj = length(cons_j)
		for ck in (cj+1):Nc
			zgmut = Zext[cj, ck]
			cons_k = cons_in_cable[ck];
			nk = length(cons_k)
			for a in 1:nj, b in 1:nk
				Ztmp[cons_j[a], cons_k[b]] += zgmut
				Ztmp[cons_k[b], cons_j[a]] += zgmut
			end
		end
	end

	stash!(ws.Z, k, Ztmp)
	return nothing
end

function compute_admittance_matrix!(
	Ptmp::AbstractMatrix{Complex{T}},
	ws,
	k::Int,
	formulation,
) where {T <: REALSCALAR}

	# Earth return (Nc×Nc)
	@inbounds fill!(Ptmp, zero(Complex{T}))
	@assert length(ws.r_ins_ext) == ws.n_phases "ws.r_ins_ext length mismatch"
	@assert length(ws.mu_ins) == ws.n_phases "ws.mu_ins length mismatch"

	Nc = ws.n_cables
	jω = ws.jω[k]

	cons_in_cable, cables = _get_cable_indices(ws)

	# Earth return admittance (Nc×Nc)
	Pext = Matrix{Complex{T}}(undef, Nc, Nc)
	compute_earth_return_matrix!(Pext, cables, ws, k, formulation.earth_admittance)
	ws.Pg[:, :, k] .= Pext # store in workspace for later use

	# --- internal Maxwell coefficients (Ametani tail-sum) -------------------------
	pinsfunctor = formulation.insulation_admittance
	@inbounds for c in 1:Nc
		cons = cons_in_cable[c]
		n = length(cons)
		if n <= 1
			continue
		end

		# gap coefficients p_g for gaps g = 1..n-1 (between cons[g] and cons[g+1])
		p = Vector{Complex{T}}(undef, n-1)
		@inbounds for g in 1:(n-1)
			i = cons[g]
			r_in_ins = ws.r_ins_in[i]     # = r_ext of conductor g
			r_ex_ins = ws.r_ins_ext[i]
			eps_ins = ws.eps_ins[i]      # eps relative
			loss_factor = ws.tan_ins[i]
			p[g] = pinsfunctor(r_in_ins, r_ex_ins, eps_ins, jω, loss_factor)
		end

		# tail sums S[k] = sum_{g=k}^{n-1} p_g, with S[n] = 0
		S = Vector{Complex{T}}(undef, n)
		S[n] = zero(Complex{T})
		@inbounds for k in (n-1):-1:1
			S[k] = p[k] + S[k+1]
		end

		# P_in[a,b] = S[max(a,b)]
		@inbounds for a in 1:n
			ia = cons[a]
			for b in 1:n
				Ptmp[ia, cons[b]] += S[max(a, b)]
			end
		end
	end
	stash!(ws.Pin, k, Ptmp)

	# stamp earth terms
	@inbounds for c in 1:Nc
		cons = cons_in_cable[c];
		n = length(cons)
		pgself = Pext[c, c]
		for a in 1:n, b in 1:n
			Ptmp[cons[a], cons[b]] += pgself
		end
	end

	@inbounds for cj in 1:(Nc-1)
		cons_j = cons_in_cable[cj];
		nj = length(cons_j)
		for ck in (cj+1):Nc
			pgmut = Pext[cj, ck]
			cons_k = cons_in_cable[ck];
			nk = length(cons_k)
			for a in 1:nj, b in 1:nk
				Ptmp[cons_j[a], cons_k[b]] += pgmut
				Ptmp[cons_k[b], cons_j[a]] += pgmut
			end
		end
	end

	stash!(ws.P, k, Ptmp)

	return nothing
end
