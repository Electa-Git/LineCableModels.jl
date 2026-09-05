assumptions(::Val{:Fortescue}) = (tolerance = 1e-4,)
"""
$(TYPEDSIGNATURES)

**Identification.** Frequency-independent unitary symmetrical-component
transform generalized to ``N`` conductors.

**Expression.**

```math
a=e^{j2\\pi/N},\\qquad
F_{km}=\\frac1{\\sqrt N}a^{km},\\quad k,m=0,\\ldots,N-1,
\\qquad \\mathbf F^{-1}=\\mathbf F^H.
```

The same unitary operator is used for voltage and current.

**Reference.** C. L. Fortescue, “Method of Symmetrical Co-Ordinates Applied
to the Solution of Polyphase Networks,” *Transactions of the AIEE*, 37,
1027–1140, 1918.
"""
description(::Formula{:Fortescue}) =
    "Fortescue symmetrical-component transformation (1918)"

# Construct the unitary Fortescue voltage and current operators.
function modal_operators(
        ::Val{:Fortescue},
        lp::LineParameters{
            Tc, U, PhaseDomain, Basis},
        ::NamedTuple
) where {Tc <: Complex, U <: Real, Basis}
    _, nph, nfreq = size(lp.Z.values)
    Tr = typeof(real(zero(Tc)))
    F = modal_basis(Val(:Fortescue), nph, Tr)
    tensor = Array{eltype(F), 3}(undef, nph, nph, nfreq)
    @inbounds for frequency in 1:nfreq
        @views copyto!(tensor[:, :, frequency], F)
    end
    return ModalOperators(tensor, tensor)
end

# Unitary N-point DFT (Fortescue) matrix.
function modal_basis(
        ::Val{:Fortescue},
        N::Integer,
        ::Type{T} = Float64
) where {T <: Real}
    N ≥ 1 || throw(ArgumentError("N ≥ 1"))
    θ = 2 * (one(T) * π) / N
    s = one(T) / sqrt(one(T) * N)
    a = cis(θ)
    return s .* [a^(k * m) for k in 0:(N - 1), m in 0:(N - 1)]  # F; inverse is F'
end

:Fortescue
