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

	for k in 1:nfreq
		Zs = symtrans(lp.Z.values[:, :, k])  # enforce reciprocity
		Ys = symtrans(lp.Y.values[:, :, k])

		Zseq = Tv * Zs * Tv'       # NOT inv(T)*Z*T — use unitary similarity
		Yseq = Tv * Ys * Tv'

		if offdiag_ratio(Zseq) > f.tol
			@warn "Fortescue: transformed Z not diagonal enough, check your results" ratio =
				offdiag_ratio(Zseq)
		end
		if offdiag_ratio(Yseq) > f.tol
			@warn "Fortescue: transformed Y not diagonal enough, check your results" ratio =
				offdiag_ratio(Yseq)
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
