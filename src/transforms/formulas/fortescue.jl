assumptions(::Val{:Fortescue}) = (tolerance = 1e-4,)
description(::Formula{:Fortescue}) = "Fortescue (symmetrical components)"

# Construct the unitary Fortescue voltage and current operators.
function (::Functor{:Fortescue})(
        lp::LineParameters{
            Tc, U, PhaseDomain, Basis},
        ::NamedTuple
) where {Tc <: Complex, U <: Real, Basis}
    _, nph, nfreq = size(lp.Z.values)
    Tr = typeof(real(zero(Tc)))
    F = fortescue_F(nph, Tr)
    tensor = Array{eltype(F), 3}(undef, nph, nph, nfreq)
    @inbounds for frequency in 1:nfreq
        @views copyto!(tensor[:, :, frequency], F)
    end
    return ModalOperators(tensor, tensor)
end

# Unitary N-point DFT (Fortescue) matrix.
function fortescue_F(N::Integer, ::Type{T} = Float64) where {T <: Real}
    N ≥ 1 || throw(ArgumentError("N ≥ 1"))
    θ = 2 * (one(T) * π) / N
    s = one(T) / sqrt(one(T) * N)
    a = cis(θ)
    return s .* [a^(k * m) for k in 0:(N - 1), m in 0:(N - 1)]  # F; inverse is F'
end

:Fortescue
