struct Levenberg <: AbstractTransformFormulation
	tol::BASE_FLOAT
end

const LMTOL = 1e-8  # default LM tolerance

# Convenient ctor
Levenberg(; tol::U = U(LMTOL)) where {U <: REALSCALAR} =
	Levenberg(BASE_FLOAT(tol))  # explicit downcast to Float64

get_description(
	::Levenberg,
) = "Levenberg–Marquardt (frequency-tracked eigen decomposition)"

"""
$(TYPEDSIGNATURES)

Apply Levenberg–Marquardt modal decomposition to a frequency-dependent
[`LineParameters`](@ref) object. Returns the (frequency-tracked) modal
transformation matrices and a **modal-domain** `LineParameters` holding the
**modal characteristic** impedance/admittance (diagonal per frequency).

# Arguments

- `lp`: Phase-domain line parameters (series `Z`, shunt `Y`, and `f`).
- `f::Levenberg`: Functor with solver tolerance.

# Returns

- `Tv`: Transformation matrices `T(•)` as a 3-tensor `n×n×nfreq` (columns are modes).
- `LineParameters`: Modal-domain characteristic parameters:
  - `Z.values[:,:,k] = Diagonal(Zcₖ)`, where `Zcₖ = sqrt.(diag(Zmₖ))./sqrt.(diag(Ymₖ))`.
  - `Y.values[:,:,k] = Diagonal(Ycₖ)`, with `Ycₖ = 1 ./ Zcₖ`.

# Notes

- Columns are **phase→modal** voltage transformation (same convention as your legacy code).
- Rotation `rot!` is applied per frequency to minimize the imaginary part of each column
  (Gustavsen’s scheme), stabilizing mode identity across the sweep.
"""
function (f::Levenberg)(lp::LineParameters)
	n, n2, nfreq = size(lp.Z.values)
	n == n2 || throw(DimensionMismatch("Z must be square"))
	size(lp.Y.values) == (n, n, nfreq) || throw(DimensionMismatch("Y must be n×n×nfreq"))

	# 1) Deterministic eigen/LM on nominal arrays
	Z_nom = to_nominal(lp.Z.values)
	Y_nom = to_nominal(lp.Y.values)
	f_nom = to_nominal(lp.f)
	Ti, _g_nom = _calc_transformation_matrix_LM(n, Z_nom, Y_nom, f_nom; tol = f.tol)
	_rot_min_imag!(Ti)

	# 2) Apply deterministic T to uncertain (or plain) inputs for *physical* outputs
	Zm, Ym, Zc_mod, Yc_mod, Zch, Ych =
		_calc_modal_quantities(Ti, lp.Z.values, lp.Y.values)
	Gdiag = _calc_gamma(Ti, lp.Z.values, lp.Y.values)

	# Keep your original return (Ti, modal characteristic) for compatibility,
	# but you now also have Zm, Ym, Zch, Ych, Gdiag available for downstream use.
	return Ti, LineParameters(SeriesImpedance(Zc_mod), ShuntAdmittance(Yc_mod), lp.f),
	LineParameters(SeriesImpedance(Zm), ShuntAdmittance(Ym), lp.f),
	LineParameters(SeriesImpedance(Zch), ShuntAdmittance(Ych), lp.f), Gdiag
end

#= ---------------------------------------------------------------------------
Internals
-----------------------------------------------------------------------------=#

# Propagate γ with uncertainty WITHOUT eigen():
# γ̂_k = sqrt.( diag( inv(T_k) * (Y_k*Z_k) * T_k ) )
function _calc_gamma(
	Ti::AbstractArray{Tc, 3},
	Z::AbstractArray{Tu, 3},
	Y::AbstractArray{Tu, 3},
) where {Tc <: Complex, Tu <: COMPLEXSCALAR}
	n, n2, nfreq = size(Ti)
	n == n2 || throw(DimensionMismatch("Ti must be n×n×nfreq"))
	size(Z) == size(Y) == (n, n, nfreq) || throw(DimensionMismatch("Z,Y must be n×n×nfreq"))

	# Element type follows uncertain inputs
	Tγ = promote_type(eltype(Z), eltype(Y))
	Gdiag = zeros(Tγ, n, n, nfreq)  # store as diagonal matrices for consistency

	Tk   = zeros(Tc, n, n)
	invT = zeros(Tc, n, n)

	@inbounds for k in 1:nfreq
		Tk   .= @view Ti[:, :, k]
		invT .= inv(Tk)

		S_k = @view(Y[:, :, k]) * @view(Z[:, :, k])         # Complex{Measurement} ok
		λdiag = diag(invT * S_k * Tk)
		γdiag = sqrt.(λdiag)
		@views Gdiag[:, :, k] .= Diagonal(γdiag)
	end
	return Gdiag
end

# Frequency-tracked Levenberg–Marquardt eigen solution
function _calc_transformation_matrix_LM(
	n::Int,
	Z::AbstractArray{T, 3},
	Y::AbstractArray{T, 3},
	f::AbstractVector{U};
	tol::U = LMTOL,
) where {T <: Complex, U <: Real}

	# Constants
	ε0 = U(ε₀)     # [F/m]
	μ0 = U(μ₀)

	nfreq = size(Z, 3)
	Ti = zeros(T, n, n, nfreq)
	g = zeros(T, n, n, nfreq)  # store as diagonalized in n×n×nfreq for convenience

	Zk = zeros(T, n, n)
	Yk = zeros(T, n, n)

	# k = 1 → plain eigen-decomposition seed
	Zk          .= @view Z[:, :, 1]
	Yk          .= @view Y[:, :, 1]
	S           = Yk * Zk
	E           = eigen(S)  # S*v = λ*v
	Ti[:, :, 1] .= E.vectors
	g[:, :, 1]  .= Diagonal(sqrt.(E.values)) # γ = sqrt(λ)

	# k ≥ 2 → LM tracking
	ord_sq = n^2
	for k in 2:nfreq
		Zk .= @view Z[:, :, k]
		Yk .= @view Y[:, :, k]

		S = Yk * Zk

		# Normalize as in legacy: (S / norm_val) - I
		ω = 2π * f[k]
		nrm = -(ω^2) * ε0 * μ0
		S̃ = (S ./ nrm) - I

		# Seed from previous step
		Tseed = @view Ti[:, :, k-1]
		gseed = @view g[:, :, k-1]
		λseed = (diag(gseed) .^ 2 ./ nrm) .- 1  # since S̃*T = T*Λ with Λ = λ̃ = (λ/nrm)-1

		# Build real-valued unknown vector: [Re(T); Im(T); Re(λ); Im(λ)]
		x0 = [
			vec(real(Tseed));
			vec(imag(Tseed));
			real(λseed);
			imag(λseed)
		]

		function _residual!(
			F::AbstractVector{<:R},
			x::AbstractVector{<:R},
		) where {R <: Real}
			# Unpack
			Tr = reshape(@view(x[1:ord_sq]), n, n)
			Ti_ = reshape(@view(x[(ord_sq+1):(2*ord_sq)]), n, n)

			λr = @view x[(2*ord_sq+1):(2*ord_sq+n)]
			λi = @view x[(2*ord_sq+n+1):(2*ord_sq+2n)]

			Λr = Diagonal(λr)
			Λi = Diagonal(λi)

			Sr = real(S̃);
			Si = imag(S̃)

			# Residual of S̃*T - T*Λ = 0, split into real/imag
			Rr = (Sr*Tr - Si*Ti_) - (Tr*Λr - Ti_*Λi)
			Ri = (Sr*Ti_ + Si*Tr) - (Tr*Λi + Ti_*Λr)

			F[1:ord_sq]              .= vec(Rr)
			F[(ord_sq+1):(2*ord_sq)] .= vec(Ri)

			# Column normalization constraints
			# For each column j: ||t_r||^2 - ||t_i||^2 = 1  and  t_r ⋅ t_i = 0
			c1 = sum(abs2.(Tr), dims = 1) .- sum(abs2.(Ti_), dims = 1) .- 1
			c2 = sum(Tr .* Ti_, dims = 1)
			idx = 2*ord_sq
			@inbounds for j in 1:n
				F[idx+2j-1] = c1[j]
				F[idx+2j]   = c2[j]
			end
			return nothing
		end

		sol = nlsolve(
			_residual!,
			x0;
			method = :trust_region,
			autodiff = :forward,
			xtol = tol,
			ftol = tol,
		)

		if !converged(sol)
			@warn "LM solver did not converge at k=$k, using seed eigen-decomposition fallback"
			E           = eigen(S)
			Ti[:, :, k] .= E.vectors
			g[:, :, k]  .= Diagonal(sqrt.(E.values))
			continue
		end

		x = sol.zero
		Tr = reshape(@view(x[1:ord_sq]), n, n)
		Ti_ = reshape(@view(x[(ord_sq+1):(2*ord_sq)]), n, n)
		T̂ = Tr .+ im .* Ti_

		λr = @view x[(2*ord_sq+1):(2*ord_sq+n)]
		λi = @view x[(2*ord_sq+n+1):(2*ord_sq+2n)]
		λ̃ = λr .+ im .* λi   # normalized eigenvalues

		# Undo normalization: λ = (λ̃ + 1) * nrm  ; γ = sqrt(λ)
		λ = (λ̃ .+ one(eltype(λ̃))) .* nrm
		γ = sqrt.(λ)

		Ti[:, :, k] .= T̂
		g[:, :, k]  .= Diagonal(γ)
	end

	return Ti, g
end

# In-place rotation to minimize imag part column-wise (per frequency slice)
function _rot_min_imag!(Ti::AbstractArray{T, 3}) where {T <: Complex}
	n, n2, nfreq = size(Ti)
	n == n2 || throw(DimensionMismatch("Ti must be n×n×nfreq"))
	tmp = zeros(T, n, n)
	@inbounds for k in 1:nfreq
		tmp .= @view Ti[:, :, k]
		rot!(tmp)                      # column-wise rotation in-place
		Ti[:, :, k] .= tmp
	end
	return Ti
end

# Full modal + characteristic + phase back-projection
# Returns:
#   Zm, Ym      :: n×n×nfreq   (modal-domain series/shunt matrices)
#   Zc_mod,Yc_mod :: n×n×nfreq (diagonal: per-mode characteristic)
#   Zch, Ych    :: n×n×nfreq   (phase-domain characteristic back-projected)
function _calc_modal_quantities(
	Ti::AbstractArray{Tc, 3},
	Z::AbstractArray{Tu, 3},
	Y::AbstractArray{Tu, 3},
) where {Tc <: Complex, Tu <: COMPLEXSCALAR}

	n, n2, nfreq = size(Ti)
	n == n2 || throw(DimensionMismatch("Ti must be n×n×nfreq"))
	size(Z) == size(Y) == (n, n, nfreq) || throw(DimensionMismatch("Z,Y must be n×n×nfreq"))

	Tz     = promote_type(eltype(Z), eltype(Y))  # keep uncertainties
	Zm     = zeros(Tz, n, n, nfreq)
	Ym     = zeros(Tz, n, n, nfreq)
	Zc_mod = zeros(Tz, n, n, nfreq)
	Yc_mod = zeros(Tz, n, n, nfreq)
	Zch    = zeros(Tz, n, n, nfreq)
	Ych    = zeros(Tz, n, n, nfreq)

	Tk   = zeros(Tc, n, n)
	Zk   = zeros(Tz, n, n)
	Yk   = zeros(Tz, n, n)
	invT = zeros(Tc, n, n)

	@inbounds for k in 1:nfreq
		Tk   .= @view Ti[:, :, k]
		invT .= inv(Tk)
		Zk   .= @view Z[:, :, k]
		Yk   .= @view Y[:, :, k]

		# Modal matrices (carry uncertainties)
		@views Zm[:, :, k] .= transpose(Tk) * Zk * Tk
		@views Ym[:, :, k] .= invT * Yk * transpose(invT)

		# Characteristic per-mode (diagonal) in modal domain
		zc = sqrt.(diag(@view Zm[:, :, k])) ./ sqrt.(diag(@view Ym[:, :, k]))
		@views Zc_mod[:, :, k] .= Diagonal(zc)
		@views Yc_mod[:, :, k] .= Diagonal(inv.(zc))

		# Phase-domain characteristic back-projection
		@views Zch[:, :, k] .= transpose(invT) * Zc_mod[:, :, k] * invT
		@views Ych[:, :, k] .= Tk * Yc_mod[:, :, k] * transpose(Tk)
	end
	return Zm, Ym, Zc_mod, Yc_mod, Zch, Ych
end

# column rotation to minimize imag parts
function rot!(S::AbstractMatrix{T}) where {T <: COMPLEXSCALAR}
	n, m = size(S)
	n == m || throw(DimensionMismatch("Input must be square"))
	@inbounds for j in 1:n
		col = @view S[:, j]

		# optimal angle
		num = -2 * sum(real.(col) .* imag.(col))             # real
		den = sum(real.(col) .^ 2 .- imag.(col) .^ 2)           # real
		ang = BASE_FLOAT(0.5) * atan(num, den)               # real

		s1 = cis(ang)
		s2 = cis(ang + BASE_FLOAT(pi/2))

		A = col .* s1
		B = col .* s2

		# all-real quadratic metrics
		Ar = real.(A);
		Ai = imag.(A)
		Br = real.(B);
		Bi = imag.(B)

		aaa1 = sum(Ai .^ 2)
		bbb1 = sum(Ar .* Ai)
		ccc1 = sum(Ar .^ 2)
		err1 = aaa1 * cos(ang)^2 + bbb1 * sin(2*ang) + ccc1 * sin(ang)^2   # real

		aaa2 = sum(Bi .^ 2)
		bbb2 = sum(Br .* Bi)
		ccc2 = sum(Br .^ 2)
		err2 = aaa2 * cos(ang)^2 + bbb2 * sin(2*ang) + ccc2 * sin(ang)^2   # real

		col .*= (err1 < err2 ? s1 : s2)
	end
	return S
end


# tiny helper: in-place imag (for metric term; avoids repeated allocations)
@inline function imag!(x::AbstractVector{<:Complex})
	@inbounds for i in eachindex(x)
		x[i] = imag(x[i])
	end
	return x
end
