struct Fortescue{T <: Real} <: AbstractTransformFormulation
    tol::T
end
Fortescue(; tol::Real = 1e-4) = Fortescue(float(tol))
description(::Fortescue) = "Fortescue (symmetrical components)"

"""
$(TYPEDSIGNATURES)

Functor implementation for `Fortescue`.
"""
function (f::Fortescue)(
        lp::LineParameters{
        Tc, U, PhaseDomain, Basis},
) where {Tc <: Complex, U <: Real, Basis}
    _, nph, nfreq = size(lp.Z.values)
    Tr = typeof(real(zero(Tc)))
    Tv = fortescue_F(nph, Tr)           # unitary; inverse is F'
    Z012 = similar(lp.Z.values)
    Y012 = similar(lp.Y.values)

    @inbounds for k in 1:nfreq
        Zs = reciprocity(lp.Z.values[:, :, k])  # enforce reciprocity
        Ys = reciprocity(lp.Y.values[:, :, k])

        Zseq = Tv * Zs * Tv'
        Yseq = Tv * Ys * Tv'

        fname = String(nameof(typeof(f)))
        offdiagZ = offdiagonal_ratio(Zseq)
        if offdiagZ > f.tol
            @warn "$fname: transformed Z not diagonal within tolerance, check your results" ratio = offdiagZ
        end
        offdiagY = offdiagonal_ratio(Yseq)
        if offdiagY > f.tol
            @warn "$fname: transformed Y not diagonal within tolerance, check your results" ratio = offdiagY
        end

        Z012[:, :, k] = Matrix(Diagonal(diag(Zseq)))
        Y012[:, :, k] = Matrix(Diagonal(diag(Yseq)))
    end
    return Tv, LineParameters(ModalDomain, Z012, Y012, lp.f; basis = basis(lp))
end

# Unitary N-point DFT (Fortescue) matrix
function fortescue_F(N::Integer, ::Type{T} = Float64) where {T <: Real}
    N ≥ 1 || throw(ArgumentError("N ≥ 1"))
    θ = 2 * (one(T) * π) / N
    s = one(T) / sqrt(one(T) * N)
    a = cis(θ)
    return s .* [a^(k * m) for k in 0:(N - 1), m in 0:(N - 1)]  # F; inverse is F'
end
