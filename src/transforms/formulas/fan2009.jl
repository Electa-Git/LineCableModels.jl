function assumptions(::Val{:Fan2009})
    return (
        tolerance = 1e-6,
        coalescence_tolerance = 1e-10,
        history_weight = 0.3
    )
end

"""
$(TYPEDSIGNATURES)

**Identification.** Postprocessed conventional eigensolutions with optimal
mode assignment, complex-phase alignment, and unitary alignment of coalescent
eigenspaces.

**Expression.** Candidate current eigenvectors minimize the assignment cost

```math
c_{pq}=1-\\frac{|\\mathbf t_p^{H}\\mathbf t_q|}
{\\lVert\\mathbf t_p\\rVert\\lVert\\mathbf t_q\\rVert}
+w_h\\frac{|\\,|\\lambda_p-\\lambda_p^-|-|\\lambda_q-\\lambda_p|\\,|}
{\\max(|\\lambda_p|,|\\lambda_p^-|,|\\lambda_q|,\\epsilon)},
```

followed by a unitary Procrustes rotation when eigenvalues are numerically
coalescent.

**Reference.** S. Fan, Y. Li, X. Li, and L. Bi, “A Method for the Calculation
of Frequency-Dependent Transmission Line Transformation Matrices,” *IEEE
Transactions on Power Systems*, 24(2), 2009.
DOI: 10.1109/TPWRS.2009.2016381.
"""
description(::Formula{:Fan2009}) =
    "Fan et al. eigenvector-tracking transformation (2009)"

"""
$(TYPEDSIGNATURES)

Calculate current modal eigenvectors with the Fan–Li–Li–Bi postprocessing
formulation. Conventional eigensolutions at adjacent frequencies are matched
by minimizing a dimensionless cost built from their normalized inner products
and eigenvalue trends. A unitary Procrustes rotation aligns numerically
coalescent eigenspaces.

# Arguments

- `lp`: Fully coupled phase-domain line parameters, with `Z` in \\[Ω/m\\], `Y`
  in \\[S/m\\], and frequency in \\[Hz\\].
- `values`: Formula assumptions containing the modal `tolerance`,
  `coalescence_tolerance`, and nonnegative `history_weight`.

# Returns

- Frequency-dependent phase-to-modal voltage and current operators.

# Notes

The assignment is solved exactly in ``O(n^3)`` operations. Complex phase is
aligned through the Hermitian inner product; it is not restricted to a real
sign change.

# Reference

S. Fan, Y. Li, X. Li, and L. Bi, *A Method for the Calculation of
Frequency-Dependent Transmission Line Transformation Matrices*, IEEE
Transactions on Power Systems, 24(2), 2009. DOI: 10.1109/TPWRS.2009.2016381.
"""
function modal_operators(
        ::Val{:Fan2009},
        lp::LineParameters{Tc, U, PhaseDomain, Basis},
        values::NamedTuple
) where {Tc <: Complex, U <: Real, Basis}
    impedance, admittance, _ = _input(lp)
    n, _, nfrequencies = size(impedance)
    T = promote_type(eltype(impedance), eltype(admittance))
    R = typeof(real(zero(T)))
    history_weight = convert(R, values.history_weight)
    isfinite(history_weight) && history_weight >= zero(R) || throw(DomainError(
        values.history_weight,
        "history_weight must be finite and nonnegative"
    ))
    coalescence_tolerance = convert(R, values.coalescence_tolerance)
    isfinite(coalescence_tolerance) && coalescence_tolerance >= zero(R) ||
        throw(DomainError(
            values.coalescence_tolerance,
            "coalescence_tolerance must be finite and nonnegative"
        ))
    current = Array{T, 3}(undef, n, n, nfrequencies)
    product = Matrix{T}(undef, n, n)
    previous_values = Vector{T}(undef, n)
    older_values = similar(previous_values)
    previous_vectors = Matrix{T}(undef, n, n)

    @inbounds for frequency in 1:nfrequencies
        _product!(product, admittance, impedance, frequency)
        current_values, current_vectors = _seed(product)
        if frequency > 1
            _match!(
                current_values,
                current_vectors,
                previous_values,
                previous_vectors;
                older_values = frequency > 2 ? older_values : nothing,
                history_weight
            )
            matched_vectors = copy(current_vectors)
            _procrustes!(
                current_vectors,
                previous_vectors,
                current_values,
                coalescence_tolerance
            )
            for mode in 1:n
                _align!(
                    @view(current_vectors[:, mode]),
                    @view(previous_vectors[:, mode])
                )
            end
            if !_valid(product, current_values, current_vectors, values.tolerance)
                copyto!(current_vectors, matched_vectors)
                for mode in 1:n
                    _align!(
                        @view(current_vectors[:, mode]),
                        @view(previous_vectors[:, mode])
                    )
                end
            end
            copyto!(older_values, previous_values)
        end
        copyto!(previous_values, current_values)
        copyto!(previous_vectors, current_vectors)
        copyto!(@view(current[:, :, frequency]), current_vectors)
    end
    return _maps(current)
end

:Fan2009
