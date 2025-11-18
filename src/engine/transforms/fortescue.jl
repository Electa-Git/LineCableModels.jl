struct Fortescue <: AbstractTransformFormulation
	tol::BASE_FLOAT
end
# Convenient constructor with default tolerance
Fortescue(; tol::BASE_FLOAT = BASE_FLOAT(1e-4)) = Fortescue(tol)
get_description(::Fortescue) = "Fortescue (symmetrical components)"

"""
$(TYPEDSIGNATURES)

Functor implementation for `Fortescue`.
"""
function (f::Fortescue)(lp::LineParameters{Tc}) where {Tc <: COMPLEXSCALAR}
	_, nph, nfreq = size(lp.Z.values)
	Tr = typeof(real(zero(Tc)))
	Tv = fortescue_F(nph, Tr)           # unitary; inverse is F'
	Z012 = similar(lp.Z.values)
	Y012 = similar(lp.Y.values)

	@inbounds for k in 1:nfreq
		Zs = symtrans(lp.Z.values[:, :, k])  # enforce reciprocity
		Ys = symtrans(lp.Y.values[:, :, k])

		Zseq = Tv * Zs * Tv'
		Yseq = Tv * Ys * Tv'

		fname = String(nameof(typeof(f)))
		offdiagZ = offdiag_ratio(Zseq)
		if offdiagZ > f.tol
			@warn "$fname: transformed Z not diagonal within tolerance, check your results" ratio =
				offdiagZ
		end
		offdiagY = offdiag_ratio(Yseq)
		if offdiagY > f.tol
			@warn "$fname: transformed Y not diagonal within tolerance, check your results" ratio =
				offdiagY
		end

		Z012[:, :, k] = Matrix(Diagonal(diag(Zseq)))
		Y012[:, :, k] = Matrix(Diagonal(diag(Yseq)))
	end
	return Tv, LineParameters(Z012, Y012, lp.f)
end

# Unitary N-point DFT (Fortescue) matrix
function fortescue_F(N::Integer, ::Type{T} = BASE_FLOAT) where {T <: REALSCALAR}
	N ≥ 1 || throw(ArgumentError("N ≥ 1"))
	θ = T(2π) / T(N)
	s = one(T) / sqrt(T(N))
	a = cis(θ)
	return s .* [a^(k * m) for k in 0:(N-1), m in 0:(N-1)]  # F; inverse is F'
end
