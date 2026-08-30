function assumptions(::Val{:Chrysochos2014})
    return (
        tolerance = 1e-6,
        convergence = 1e-8,
        max_iterations = 100,
        damping = 1e-3
    )
end

"""
$(TYPEDSIGNATURES)

**Identification.** Levenberg–Marquardt tracking of each complex eigenpair,
initialized from the preceding frequency.

**Expression.**

```math
\\widetilde{\\mathbf S}=\\frac{\\mathbf Y\\mathbf Z}
{-\\omega^2\\mu_0\\varepsilon_0}-\\mathbf I,\\qquad
\\widetilde{\\mathbf S}\\mathbf t=\\lambda\\mathbf t,qquad
\\mathbf t^T\\mathbf t=1.
```

The complex residual is represented as a real least-squares system and solved
with a damped normal-equation step.

**Reference.** A. I. Chrysochos, T. A. Papadopoulos, and G. K. Papagiannis,
“Robust Calculation of Frequency-Dependent Transmission-Line Transformation
Matrices Using the Levenberg–Marquardt Method,” *IEEE Transactions on Power
Delivery*, 29(4), 1621–1629, 2014. DOI: 10.1109/TPWRD.2013.2284504.
"""
description(::Formula{:Chrysochos2014}) =
    "Chrysochos et al. 2014 (Levenberg–Marquardt eigenpair tracking)"

function _lmres!(
        residual::AbstractVector{R},
        x::AbstractVector{R},
        real_matrix::AbstractMatrix{R},
        imaginary_matrix::AbstractMatrix{R}
) where {R <: Real}
    n = size(real_matrix, 1)
    real_vector = @view x[1:n]
    imaginary_vector = @view x[(n + 1):(2n)]
    real_value = x[2n + 1]
    imaginary_value = x[2n + 2]
    real_constraint = -one(R)
    imaginary_constraint = zero(R)
    @inbounds for row in 1:n
        real_result = -real_value * real_vector[row] +
                      imaginary_value * imaginary_vector[row]
        imaginary_result = -imaginary_value * real_vector[row] -
                           real_value * imaginary_vector[row]
        for column in 1:n
            real_result += real_matrix[row, column] * real_vector[column] -
                           imaginary_matrix[row, column] * imaginary_vector[column]
            imaginary_result += imaginary_matrix[row, column] * real_vector[column] +
                                real_matrix[row, column] * imaginary_vector[column]
        end
        residual[row] = real_result
        residual[n + row] = imaginary_result
        real_constraint += real_vector[row]^2 - imaginary_vector[row]^2
        imaginary_constraint += 2 * real_vector[row] * imaginary_vector[row]
    end
    residual[2n + 1] = real_constraint
    residual[2n + 2] = imaginary_constraint
    return residual
end

function _lmjac!(
        jacobian::AbstractMatrix{R},
        x::AbstractVector{R},
        real_matrix::AbstractMatrix{R},
        imaginary_matrix::AbstractMatrix{R}
) where {R <: Real}
    n = size(real_matrix, 1)
    real_vector = @view x[1:n]
    imaginary_vector = @view x[(n + 1):(2n)]
    real_value = x[2n + 1]
    imaginary_value = x[2n + 2]
    fill!(jacobian, zero(R))
    @inbounds for row in 1:n
        for column in 1:n
            diagonal = row == column
            jacobian[row, column] = real_matrix[row, column] -
                                    (diagonal ? real_value : zero(R))
            jacobian[row, n + column] = -imaginary_matrix[row, column] +
                                        (diagonal ? imaginary_value : zero(R))
            jacobian[n + row, column] = imaginary_matrix[row, column] -
                                        (diagonal ? imaginary_value : zero(R))
            jacobian[n + row, n + column] = real_matrix[row, column] -
                                            (diagonal ? real_value : zero(R))
        end
        jacobian[row, 2n + 1] = -real_vector[row]
        jacobian[row, 2n + 2] = imaginary_vector[row]
        jacobian[n + row, 2n + 1] = -imaginary_vector[row]
        jacobian[n + row, 2n + 2] = -real_vector[row]
        jacobian[2n + 1, row] = 2 * real_vector[row]
        jacobian[2n + 1, n + row] = -2 * imaginary_vector[row]
        jacobian[2n + 2, row] = 2 * imaginary_vector[row]
        jacobian[2n + 2, n + row] = 2 * real_vector[row]
    end
    return jacobian
end

function _lmwork(::Type{R}, n::Integer) where {R <: Real}
    order = 2n + 2
    return (
        x = Vector{R}(undef, order),
        candidate = Vector{R}(undef, order),
        residual = Vector{R}(undef, order),
        candidate_residual = Vector{R}(undef, order),
        jacobian = Matrix{R}(undef, order, order),
        gradient = Vector{R}(undef, order),
        hessian = Matrix{R}(undef, order, order),
        system = Matrix{R}(undef, order, order),
        step = Vector{R}(undef, order),
        real_matrix = Matrix{R}(undef, n, n),
        imaginary_matrix = Matrix{R}(undef, n, n)
    )
end

function _lm!(
        vector::AbstractVector{T},
        value::T,
        values::NamedTuple,
        work
) where {T <: Complex}
    n = length(vector)
    R = typeof(real(zero(T)))
    requested_tolerance = convert(R, values.convergence)
    tolerance = max(R(100) * eps(R), requested_tolerance^2)
    iterations = values.max_iterations
    iterations isa Integer && iterations > 0 || throw(DomainError(
        iterations,
        "max_iterations must be a positive integer"
    ))
    damping = convert(R, values.damping)
    isfinite(damping) && damping > zero(R) || throw(DomainError(
        damping,
        "damping must be finite and positive"
    ))

    _transpose!(vector) || _unit!(vector) || return value, false
    @inbounds for index in 1:n
        work.x[index] = real(vector[index])
        work.x[n + index] = imag(vector[index])
    end
    work.x[2n + 1] = real(value)
    work.x[2n + 2] = imag(value)

    converged = false
    maximum_damping = inv(eps(R))
    for _ in 1:iterations
        _lmres!(
            work.residual,
            work.x,
            work.real_matrix,
            work.imaginary_matrix
        )
        residual_norm = norm(work.residual, Inf)
        if residual_norm <= tolerance
            converged = true
            break
        end
        _lmjac!(
            work.jacobian,
            work.x,
            work.real_matrix,
            work.imaginary_matrix
        )
        mul!(work.gradient, transpose(work.jacobian), work.residual)
        norm(work.gradient, Inf) <= tolerance && (converged = true; break)
        mul!(work.hessian, transpose(work.jacobian), work.jacobian)
        copyto!(work.system, work.hessian)
        diagonal_scale = one(R)
        @inbounds for index in axes(work.hessian, 1)
            diagonal_scale = max(
                diagonal_scale,
                abs(work.hessian[index, index])
            )
        end
        floor = eps(R) * diagonal_scale
        @inbounds for index in axes(work.system, 1)
            work.system[index, index] += damping *
                                         max(work.hessian[index, index], floor)
        end
        copyto!(work.step, work.gradient)
        work.step .*= -one(R)
        factorization = lu!(work.system; check = false)
        issuccess(factorization) || return value, false
        ldiv!(factorization, work.step)
        all(isfinite, work.step) || return value, false
        work.candidate .= work.x .+ work.step
        _lmres!(
            work.candidate_residual,
            work.candidate,
            work.real_matrix,
            work.imaginary_matrix
        )
        old_cost = dot(work.residual, work.residual)
        new_cost = dot(work.candidate_residual, work.candidate_residual)
        if isfinite(new_cost) && new_cost < old_cost
            copyto!(work.x, work.candidate)
            damping = max(damping / R(3), eps(R))
            if norm(work.step, Inf) <= tolerance *
               (tolerance + norm(work.x, Inf))
                converged = true
                break
            end
        else
            damping *= R(10)
            damping <= maximum_damping || return value, false
        end
    end

    @inbounds for index in 1:n
        vector[index] = complex(work.x[index], work.x[n + index])
    end
    _unit!(vector) || return value, false
    return complex(work.x[2n + 1], work.x[2n + 2]), converged
end

"""
$(TYPEDSIGNATURES)

Calculate current modal eigenvectors with the Chrysochos–Papadopoulos–
Papagiannis Levenberg–Marquardt formulation. At each frequency, the current
eigenproblem is scaled as

```math
\\widetilde{S} = \\frac{YZ}{-\\omega^2\\mu_0\\epsilon_0} - I, \\qquad
\\widetilde{S}t = \\lambda t,
```

then each complex eigenpair is solved independently through its equivalent
real residual with the constraint ``t^Tt=1``. The preceding frequency supplies
the initial eigenpair.

# Arguments

- `lp`: Fully coupled phase-domain line parameters, with `Z` in \\[Ω/m\\], `Y`
  in \\[S/m\\], and frequency in \\[Hz\\].
- `values`: Formula assumptions containing the modal-residue `tolerance`, LM
  `convergence`, `max_iterations`, and initial `damping` coefficient.

# Returns

- Frequency-dependent phase-to-modal voltage and current operators.

# Notes

The implementation uses an analytic real Jacobian and a monotone damped
normal-equation step. A conventional eigensolution matched to the preceding
frequency is retained if an iterative slice is singular or does not converge.

# Reference

A. I. Chrysochos, T. A. Papadopoulos, and G. K. Papagiannis, *Robust
Calculation of Frequency-Dependent Transmission-Line Transformation Matrices
Using the Levenberg–Marquardt Method*, IEEE Transactions on Power Delivery,
29(4), 2014. DOI: 10.1109/TPWRD.2013.2284504.
"""
function chrysochos2014(
        lp::LineParameters{Tc, U, PhaseDomain, Basis},
        values::NamedTuple
) where {Tc <: Complex, U <: Real, Basis}
    impedance, admittance, frequencies = _input(lp)
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
    damping = convert(R, values.damping)
    isfinite(damping) && damping > zero(R) || throw(DomainError(
        values.damping,
        "damping must be finite and positive"
    ))
    current = Array{T, 3}(undef, n, n, nfrequencies)
    product = Matrix{T}(undef, n, n)
    normalized = similar(product)
    work = _lmwork(R, n)

    previous_values = Vector{T}(undef, n)
    previous_vectors = Matrix{T}(undef, n, n)
    current_values = similar(previous_values)
    physical_values = similar(previous_values)
    current_vectors = similar(previous_vectors)
    fallback_count = 0
    first_fallback = 0
    validation_tolerance = max(
        R(100) * eps(R),
        convergence^2
    )

    @inbounds for frequency_index in 1:nfrequencies
        frequency = convert(R, frequencies[frequency_index])
        frequency > zero(R) || throw(DomainError(
            frequency,
            "Chrysochos2014 requires positive frequencies"
        ))
        _product!(product, admittance, impedance, frequency_index)
        unit = one(frequency)
        omega = 2 * (unit * π) * frequency
        epsilon0 = unit * 88541878128 * (unit * 10)^(-22)
        mu0 = unit * 4 * (unit * π) * (unit * 10)^(-7)
        scale = -(omega^2) * epsilon0 * mu0
        normalized .= product ./ scale
        for mode in 1:n
            normalized[mode, mode] -= one(T)
        end
        for index in eachindex(normalized)
            work.real_matrix[index] = real(normalized[index])
            work.imaginary_matrix[index] = imag(normalized[index])
        end

        if frequency_index == 1
            seed_values, seed_vectors = _seed(normalized)
            copyto!(previous_values, seed_values)
            copyto!(previous_vectors, seed_vectors)
        else
            copyto!(current_values, previous_values)
            copyto!(current_vectors, previous_vectors)
            failed_mode = 0
            for mode in 1:n
                reference = @view previous_vectors[:, mode]
                value, converged = _lm!(
                    @view(current_vectors[:, mode]),
                    previous_values[mode],
                    values,
                    work
                )
                current_values[mode] = value
                _align!(@view(current_vectors[:, mode]), reference)
                if !converged
                    failed_mode = mode
                    break
                end
            end
            for mode in 1:n
                physical_values[mode] =
                    (current_values[mode] + one(T)) * scale
            end
            if failed_mode != 0 || !_valid(
                product,
                physical_values,
                current_vectors,
                validation_tolerance
            )
                fallback_values, fallback_vectors = _fallback(
                    product,
                    previous_values,
                    previous_vectors
                )
                for mode in 1:n
                    current_values[mode] =
                        fallback_values[mode] / scale - one(T)
                end
                copyto!(current_vectors, fallback_vectors)
                fallback_count += 1
                iszero(first_fallback) && (first_fallback = frequency_index)
            end
            copyto!(previous_values, current_values)
            copyto!(previous_vectors, current_vectors)
        end
        copyto!(@view(current[:, :, frequency_index]), previous_vectors)
    end
    if fallback_count > 0
        @warn ":Chrysochos2014 retained matched eigensolutions where LM lost relative accuracy" fallback_count first_fallback
    end
    return _maps(current)
end

function (::Functor{:Chrysochos2014})(
        lp::LineParameters{Tc, U, PhaseDomain, Basis},
        values::NamedTuple
) where {Tc <: Complex, U <: Real, Basis}
    return chrysochos2014(lp, values)
end

:Chrysochos2014
