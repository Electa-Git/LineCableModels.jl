function assumptions(::Val{:Wedepohl1996})
    return (
        tolerance = 1e-6,
        convergence = 1e-9,
        max_iterations = 100
    )
end

"""
$(TYPEDSIGNATURES)

**Identification.** Complex Newton–Raphson eigenpair tracking, initialized
from the preceding frequency and scaled by the Frobenius norm of
``\\mathbf Y\\mathbf Z``.

**Expression.** Each mode solves

```math
\\mathbf S\\mathbf t-\\lambda\\mathbf t=0,\\qquad
\\mathbf t^T\\mathbf t-1=0,
```

with the complex eigenpair and normalization constraint advanced jointly by
Newton iteration.

**Reference.** L. M. Wedepohl, H. V. Nguyen, and G. D. Irwin,
“Frequency-Dependent Transformation Matrices for Untransposed Transmission
Lines Using Newton-Raphson Method,” *IEEE Transactions on Power Systems*,
11(3), 1996. DOI: 10.1109/59.535695.
"""
description(::Formula{:Wedepohl1996}) =
    "Wedepohl et al. Newton–Raphson modal transformation (1996)"

function newton_workspace(
        ::Val{:Wedepohl1996},
        ::Type{T},
        n::Integer
) where {T <: Complex}
    return (
        residual = Vector{T}(undef, n + 1),
        jacobian = Matrix{T}(undef, n + 1, n + 1),
        step = Vector{T}(undef, n + 1)
    )
end

function newton_step!(
        ::Val{:Wedepohl1996},
        vector::AbstractVector{T},
        value::T,
        matrix::AbstractMatrix{T},
        values::NamedTuple,
        work
) where {T <: Complex}
    n = length(vector)
    R = typeof(real(zero(T)))
    tolerance = convert(R, values.convergence)
    iterations = values.max_iterations
    iterations isa Integer && iterations > 0 || throw(DomainError(
        iterations,
        "max_iterations must be a positive integer"
    ))
    _transpose!(vector) || _unit!(vector) || return value, false

    converged = false
    @inbounds for _ in 1:iterations
        fill!(work.jacobian, zero(T))
        constraint = -one(T)
        for row in 1:n
            result = -value * vector[row]
            for column in 1:n
                result += matrix[row, column] * vector[column]
                work.jacobian[row, column] = matrix[row, column] -
                                             (row == column ? value : zero(T))
            end
            work.residual[row] = result
            work.jacobian[row, n + 1] = -vector[row]
            work.jacobian[n + 1, row] = 2 * vector[row]
            constraint += vector[row] * vector[row]
        end
        work.residual[n + 1] = constraint
        work.jacobian[n + 1, n + 1] = zero(T)
        if norm(work.residual, Inf) <= tolerance
            converged = true
            break
        end
        copyto!(work.step, work.residual)
        factorization = lu!(work.jacobian; check = false)
        issuccess(factorization) || return value, false
        ldiv!(factorization, work.step)
        all(isfinite, work.step) || return value, false
        vector .-= @view work.step[1:n]
        value -= work.step[n + 1]
        if norm(work.step, Inf) <= tolerance *
           (tolerance + max(norm(vector, Inf), abs(value)))
            converged = true
            break
        end
    end
    _unit!(vector) || return value, false
    return value, converged
end

"""
$(TYPEDSIGNATURES)

Calculate current modal eigenvectors with the Wedepohl–Nguyen–Irwin complex
Newton–Raphson formulation. Each normalized eigenpair satisfies

```math
S t - \\lambda t = 0, \\qquad t^Tt - 1 = 0,
```

and is initialized from the corresponding solution at the preceding
frequency.

# Arguments

- `lp`: Fully coupled phase-domain line parameters, with `Z` in \\[Ω/m\\], `Y`
  in \\[S/m\\], and frequency in \\[Hz\\].
- `values`: Formula assumptions containing the modal-residue `tolerance`,
  Newton `convergence`, and `max_iterations`.

# Returns

- Frequency-dependent phase-to-modal voltage and current operators.

# Notes

Every ``YZ`` slice is scaled by its Frobenius norm. The preceding physical
eigenvalue is divided by the current scale before Newton iteration, avoiding
the inconsistent seed present in the source notebook. A matched conventional
eigensolution is retained when Newton iteration fails.

# Reference

L. M. Wedepohl, H. V. Nguyen, and G. D. Irwin, *Frequency-Dependent
Transformation Matrices for Untransposed Transmission Lines Using
Newton-Raphson Method*, IEEE Transactions on Power Systems, 11(3), 1996.
DOI: 10.1109/59.535695.
"""
function modal_operators(
        ::Val{:Wedepohl1996},
        lp::LineParameters{Tc, U, PhaseDomain, Basis},
        values::NamedTuple
) where {Tc <: Complex, U <: Real, Basis}
    impedance, admittance, _ = _input(lp)
    n, _, nfrequencies = size(impedance)
    T = promote_type(eltype(impedance), eltype(admittance))
    R = typeof(real(zero(T)))
    convergence = convert(R, values.convergence)
    isfinite(convergence) && convergence > zero(R) || throw(DomainError(
        values.convergence,
        "convergence must be finite and positive"
    ))
    iterations = values.max_iterations
    iterations isa Integer && iterations > 0 || throw(DomainError(
        iterations,
        "max_iterations must be a positive integer"
    ))
    current = Array{T, 3}(undef, n, n, nfrequencies)
    product = Matrix{T}(undef, n, n)
    normalized = similar(product)
    work = newton_workspace(Val(:Wedepohl1996), T, n)
    previous_values = Vector{T}(undef, n)
    previous_vectors = Matrix{T}(undef, n, n)
    current_values = similar(previous_values)
    current_vectors = similar(previous_vectors)
    fallback_count = 0
    first_fallback = 0
    first_failed_mode = 0

    @inbounds for frequency in 1:nfrequencies
        _product!(product, admittance, impedance, frequency)
        if frequency == 1
            seed_values, seed_vectors = _seed(product)
            copyto!(previous_values, seed_values)
            copyto!(previous_vectors, seed_vectors)
        else
            scale = norm(product)
            isfinite(scale) && scale > zero(R) || throw(DomainError(
                scale,
                ":Wedepohl1996 requires a nonzero finite YZ norm"
            ))
            normalized .= product ./ scale
            copyto!(current_vectors, previous_vectors)
            failed_mode = 0
            for mode in 1:n
                reference = @view previous_vectors[:, mode]
                value, converged = newton_step!(
                    Val(:Wedepohl1996),
                    @view(current_vectors[:, mode]),
                    previous_values[mode] / scale,
                    normalized,
                    values,
                    work
                )
                current_values[mode] = value * scale
                _align!(@view(current_vectors[:, mode]), reference)
                if !converged
                    failed_mode = mode
                    break
                end
            end
            if failed_mode != 0 || !_valid(
                product,
                current_values,
                current_vectors,
                convergence
            )
                fallback_values, fallback_vectors = _fallback(
                    product,
                    previous_values,
                    previous_vectors
                )
                copyto!(current_values, fallback_values)
                copyto!(current_vectors, fallback_vectors)
                fallback_count += 1
                if iszero(first_fallback)
                    first_fallback = frequency
                    first_failed_mode = failed_mode
                end
            end
            copyto!(previous_values, current_values)
            copyto!(previous_vectors, current_vectors)
        end
        copyto!(@view(current[:, :, frequency]), previous_vectors)
    end
    if fallback_count > 0
        @warn ":Wedepohl1996 retained matched eigensolutions where Newton–Raphson did not converge" fallback_count first_fallback first_failed_mode
    end
    return _maps(current)
end

:Wedepohl1996
