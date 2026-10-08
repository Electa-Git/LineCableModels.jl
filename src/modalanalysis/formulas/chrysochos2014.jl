
"""
$(TYPEDSIGNATURES)

**Identification.** Levenberg–Marquardt tracking of each complex eigenpair,
initialized from the preceding frequency.

**Expression.**

```math
\\widetilde{\\mathbf S}=\\frac{\\mathbf Y\\mathbf Z}
{-\\omega^2\\mu_0\\varepsilon_0}-\\mathbf I,\\qquad
\\widetilde{\\mathbf S}\\mathbf t=\\lambda\\mathbf t,\\qquad
\\mathbf t^T\\mathbf t=1.
```

The complex residual is represented as a real least-squares system and solved
with a damped normal-equation step.

**Reference.** A. I. Chrysochos, T. A. Papadopoulos, and G. K. Papagiannis,
“Robust Calculation of Frequency-Dependent Transmission-Line Transformation
Matrices Using the Levenberg–Marquardt Method,” *IEEE Transactions on Power
Delivery*, 29(4), 1621–1629, 2014. DOI: 10.1109/TPWRD.2013.2284504.

"""
function description(::Type{<:Formula{:chrysochos2014}}; compact::Bool=false)
    compact ? "Chrysochos" : "Chrysochos et al. Levenberg–Marquardt modal transformation (2014)"
end

function levenberg_marquardt_residual!(
        ::Formula{:chrysochos2014},
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

function levenberg_marquardt_jacobian!(
        ::Formula{:chrysochos2014},
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

function levenberg_marquardt_step!(
        formula::Formula{:chrysochos2014},
        vector::AbstractVector{T},
        value::T,
        iteration_options::NamedTuple,
        buffers
) where {T <: Complex}
    n = length(vector)
    R = typeof(real(zero(T)))
    requested_tolerance = convert(R, iteration_options.convergence)
    tolerance = max(R(100) * eps(R), requested_tolerance^2)
    iterations = iteration_options.max_iterations
    iterations isa Integer && iterations > 0 || throw(DomainError(
        iterations,
        "max_iterations must be a positive integer"
    ))
    damping = convert(R, iteration_options.damping)
    isfinite(damping) && damping > zero(R) || throw(DomainError(
        damping,
        "damping must be finite and positive"
    ))

    normalize_bilinear!(vector) || _unit!(vector) || return value, false, 0
    @inbounds for index in 1:n
        buffers.x[index] = real(vector[index])
        buffers.x[n + index] = imag(vector[index])
    end
    buffers.x[2n + 1] = real(value)
    buffers.x[2n + 2] = imag(value)

    converged = false
    performed = 0
    maximum_damping = inv(eps(R))
    for iteration in 1:iterations
        performed=iteration
        levenberg_marquardt_residual!(
            formula,
            buffers.residual,
            buffers.x,
            buffers.real_matrix,
            buffers.imaginary_matrix
        )
        residual_norm = norm(buffers.residual, Inf)
        if residual_norm <= tolerance
            converged = true
            break
        end
        levenberg_marquardt_jacobian!(
            formula,
            buffers.jacobian,
            buffers.x,
            buffers.real_matrix,
            buffers.imaginary_matrix
        )
        mul!(buffers.gradient, transpose(buffers.jacobian), buffers.residual)
        norm(buffers.gradient, Inf) <= tolerance && (converged = true; break)
        mul!(buffers.hessian, transpose(buffers.jacobian), buffers.jacobian)
        copyto!(buffers.system, buffers.hessian)
        diagonal_scale = one(R)
        @inbounds for index in axes(buffers.hessian, 1)
            diagonal_scale = max(
                diagonal_scale,
                abs(buffers.hessian[index, index])
            )
        end
        floor = eps(R) * diagonal_scale
        @inbounds for index in axes(buffers.system, 1)
            buffers.system[index, index] += damping *
                                         max(buffers.hessian[index, index], floor)
        end
        copyto!(buffers.step, buffers.gradient)
        buffers.step .*= -one(R)
        factorization = lu!(buffers.system; check = false)
        issuccess(factorization) || return value, false, performed
        ldiv!(factorization, buffers.step)
        all(isfinite, buffers.step) || return value, false, performed
        buffers.candidate .= buffers.x .+ buffers.step
        levenberg_marquardt_residual!(
            formula,
            buffers.candidate_residual,
            buffers.candidate,
            buffers.real_matrix,
            buffers.imaginary_matrix
        )
        old_cost = dot(buffers.residual, buffers.residual)
        new_cost = dot(buffers.candidate_residual, buffers.candidate_residual)
        if isfinite(new_cost) && new_cost < old_cost
            copyto!(buffers.x, buffers.candidate)
            damping = max(damping / R(3), eps(R))
            if norm(buffers.step, Inf) <= tolerance *
                                       (tolerance + norm(buffers.x, Inf))
                converged = true
                break
            end
        else
            damping *= R(10)
            damping <= maximum_damping || return value, false, performed
        end
    end

    @inbounds for index in 1:n
        vector[index] = complex(buffers.x[index], buffers.x[n + index])
    end
    _unit!(vector) || return value, false, performed
    return complex(buffers.x[2n + 1], buffers.x[2n + 2]), converged, performed
end

"""
$(TYPEDSIGNATURES)

Calculate current modal eigenvectors with the Chrysochos-Papadopoulos,
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

- `lp`: fully coupled phase-domain line parameters, with `Z` in \\[Ω/m\\], `Y`
  in \\[S/m\\], and frequency in \\[Hz\\].
- `values`: modal controls containing the modal-residue `tolerance`, LM
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
function initialize_buffers(::Formula{:chrysochos2014}, ::Type{T}, input,
        plan, buffers) where {T <: Complex}
    n = plan.n
    R = typeof(real(zero(T)))
    order = 2n + 2
    return merge(buffers,(
        normalized_shifted_eigenproblem=Matrix{T}(undef,n,n),
        least_squares=(
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
            imaginary_matrix = Matrix{R}(undef, n, n)),
        eigenpair_assignment=(cost=Matrix{R}(undef,n,n),u=zeros(R,n+1),v=zeros(R,n+1),
            matching=zeros(Int,n+1),way=zeros(Int,n+1),minimums=Vector{R}(undef,n+1),
            used=falses(n+1),assignment=Vector{Int}(undef,n),
            ordered_values=Vector{T}(undef,n),ordered_vectors=Matrix{T}(undef,n,n),
            residual=Vector{T}(undef,n))))
end

function decompose!(formula::Formula{:chrysochos2014}, workspace::ModalAnalysisWorkspace,
        parameters::NamedTuple, options::FormulationOptions)
    iteration_options = options.data.iteration
    impedance = workspace.input.Z
    admittance = workspace.input.Y
    frequencies = workspace.input.f
    n, _, nfrequencies = size(impedance)
    T = eltype(workspace.Ti)
    R = typeof(real(zero(T)))
    buffers = workspace.buffers
    admittance_impedance_product = buffers.admittance_impedance_product
    normalized_shifted_eigenproblem = buffers.normalized_shifted_eigenproblem
    least_squares = buffers.least_squares
    previous_eigenvalues = buffers.previous_eigenvalues
    previous_eigenvectors = buffers.previous_eigenvectors
    eigenvalues = buffers.eigenvalues
    eigenvectors = buffers.eigenvectors
    propagation_eigenvalues = buffers.propagation_eigenvalues
    convergence = convert(R, iteration_options.convergence)
    validation_tolerance = max(R(100)*eps(R), convergence^2)
    missed = workspace.diagnostics.missed_frequencies
    fallback = workspace.diagnostics.fallback_frequencies

    @inbounds for frequency_index in 1:nfrequencies
        copyto!(workspace.buffers.Zslice,@view(impedance[:,:,frequency_index]))
        copyto!(workspace.buffers.Yslice,@view(admittance[:,:,frequency_index]))
        workspace.buffers.Zslice ./= workspace.input.root_scale
        workspace.buffers.Yslice ./= workspace.input.root_scale
        Zslice=workspace.buffers.Zslice
        Yslice=workspace.buffers.Yslice
        frequency = convert(R, frequencies[frequency_index])
        frequency > zero(R) || throw(DomainError(frequency,
            "Chrysochos modal analysis requires positive frequencies"))
        mul!(admittance_impedance_product,Yslice,Zslice)
        unit = one(frequency)
        omega = 2*(unit*π)*frequency
        epsilon0 = vacuum_permittivity(typeof(frequency))
        mu0 = vacuum_permeability(typeof(frequency))
        scale = -(omega^2)*epsilon0*mu0
        normalized_shifted_eigenproblem .= admittance_impedance_product ./ scale
        for mode in 1:n
            normalized_shifted_eigenproblem[mode,mode] -= one(T)
        end
        for index in eachindex(normalized_shifted_eigenproblem)
            least_squares.real_matrix[index] = real(normalized_shifted_eigenproblem[index])
            least_squares.imaginary_matrix[index] = imag(normalized_shifted_eigenproblem[index])
        end
        if frequency_index == 1
            seed_values, seed_vectors = _seed(normalized_shifted_eigenproblem)
            copyto!(previous_eigenvalues, seed_values)
            copyto!(previous_eigenvectors, seed_vectors)
        else
            copyto!(eigenvalues, previous_eigenvalues)
            copyto!(eigenvectors, previous_eigenvectors)
            failed_mode = 0
            for mode in 1:n
                reference = @view previous_eigenvectors[:,mode]
                value, converged, iterations = levenberg_marquardt_step!(
                    formula, @view(eigenvectors[:,mode]),
                    previous_eigenvalues[mode], iteration_options, least_squares)
                workspace.diagnostics.iterations[mode,frequency_index]=iterations
                workspace.diagnostics.converged[mode,frequency_index]=converged
                eigenvalues[mode] = value
                _align!(@view(eigenvectors[:,mode]), reference)
                if !converged
                    failed_mode = mode
                    break
                end
            end
            for mode in 1:n
                propagation_eigenvalues[mode] = (eigenvalues[mode]+one(T))*scale
            end
            if failed_mode != 0 || !diagonalizes!(eigenvectors,
                    propagation_eigenvalues,admittance_impedance_product,validation_tolerance,
                    buffers.eigenpair_assignment.residual)
                push!(missed, frequency_index)
                if iteration_options.fallback === :matched
                    push!(fallback, frequency_index)
                    fallback_values, fallback_vectors = recompute_matched_eigenpairs!(
                        admittance_impedance_product,previous_eigenvalues,
                        previous_eigenvectors,buffers.eigenpair_assignment)
                    for mode in 1:n
                        eigenvalues[mode] = fallback_values[mode]/scale-one(T)
                    end
                    copyto!(eigenvectors, fallback_vectors)
                end
            end
            copyto!(previous_eigenvalues, eigenvalues)
            copyto!(previous_eigenvectors, eigenvectors)
        end
        copyto!(@view(workspace.Ti[:,:,frequency_index]),previous_eigenvectors)
        for mode in 1:n
            vector = @view workspace.Ti[:,mode,frequency_index]
            _unit!(vector) || throw(ArgumentError("current eigenvector has zero norm"))
            mul!(buffers.voltage_vector,Zslice,vector)
            divisor = norm(buffers.voltage_vector)
            isfinite(divisor) && !iszero(divisor) ||
                throw(ArgumentError("voltage eigenvector has zero or undefined norm"))
            @views workspace.Tv[:,mode,frequency_index] .= buffers.voltage_vector ./ divisor
            root = sqrt((previous_eigenvalues[mode]+one(T))*scale)
            (real(root)<0 || (iszero(real(root)) && imag(root)<0)) && (root=-root)
            workspace.roots[mode,frequency_index] = root*workspace.input.root_scale
            eigenvalue=(previous_eigenvalues[mode]+one(T))*scale
            mul!(buffers.eigenpair_assignment.residual,admittance_impedance_product,vector)
            buffers.eigenpair_assignment.residual .-= eigenvalue .* vector
            denominator=(norm(admittance_impedance_product,Inf)+abs(eigenvalue))*norm(vector,Inf)
            numerator=norm(buffers.eigenpair_assignment.residual,Inf)
            workspace.diagnostics.eigen_residual[mode,frequency_index]=
                iszero(denominator) ? (iszero(numerator) ? zero(R) : R(Inf)) :
                numerator/denominator
        end
    end
    isempty(missed) || @warn ":chrysochos2014 missed numerical targets" count=length(missed) frequencies=copy(missed) fallback_count=length(fallback)
    return workspace
end

function formulation_options(::Expression{<:Formula{:chrysochos2014}, typeof(decompose!)})
    return FormulationOptions((iteration = (
        convergence = 1e-8, max_iterations = 100, damping = 1e-3, fallback = :matched),))
end

function formulation_options(::Expression{<:Formula{:chrysochos2014}, typeof(decompose!)},
        ::Val{:iteration}, defaults::NamedTuple, supplied::NamedTuple)
    isempty(setdiff(keys(supplied), keys(defaults))) ||
        throw(ArgumentError("unknown modal iteration controls"))
    options = merge(defaults, supplied)
    for key in (:convergence, :damping)
        value = getproperty(options, key)
        value isa Real && isfinite(value) && value > 0 ||
            throw(ArgumentError("$key must be finite and positive"))
    end
    options.max_iterations isa Integer && !(options.max_iterations isa Bool) &&
    options.max_iterations > 0 ||
        throw(ArgumentError("max_iterations must be a positive integer"))
    options.fallback in (:matched, :none) ||
        throw(ArgumentError("fallback must be :matched or :none"))
    return options
end

:chrysochos2014
