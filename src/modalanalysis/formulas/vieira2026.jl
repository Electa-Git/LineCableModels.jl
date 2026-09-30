"""
$(TYPEDSIGNATURES)

Track current eigenpairs with Vieira's complex implementation of the
Chrysochos–Papadopoulos–Papagiannis Levenberg–Marquardt method:

```math
S = \\frac{YZ}{-(2\\pi f)^2\\mu_0\\epsilon_0}-I,\\qquad
r(t,\\lambda) = \\begin{bmatrix}(S-\\lambda I)t\\\\t^Tt-1\\end{bmatrix},\\qquad
J = \\begin{bmatrix}S-\\lambda I&-t\\\\2t^T&0\\end{bmatrix}.
```

Here `Z` and `Y` are per-length matrices in \\[Ω/m\\] and \\[S/m\\], and `f`
is frequency in \\[Hz\\]. Complex SVD least-squares steps minimize the Euclidean
residual norm. A two-sample predictor, greedy eigenvalue matching and clustered
eigenvector refinement track the supplied samples. Failed tracking uses a
correlation-matched eigensolution at the same frequency; no samples are inserted.
Internal vectors use the bilinear constraint; returned `Ti` and paired `Tv`
columns have unit Euclidean norm and map modal quantities to phase quantities.

`iteration` controls `convergence=1e-13` and `max_iterations=100`. `tracking`
controls `predictor_tolerance=0.25`, `eigenvalue_tolerance=1e-5`,
`cluster_tolerance=1e-3` and `order_by_velocity=true`. The latter orders all
samples by decreasing phase constant at the final frequency. Numerical misses
are diagnostics, not rejection of finite results.

Adapted from P. H. N. Vieira's `eig_levenberg_marquardt` implementation in
[Parametric-Eigenvalues-and-Vectors-of-Transmission-Lines](https://github.com/pedrohnv/Parametric-Eigenvalues-and-Vectors-of-Transmission-Lines).
The method is from A. I. Chrysochos, T. A. Papadopoulos and G. K. Papagiannis, *Robust Calculation of Frequency-Dependent Transmission-Line Transformation Matrices Using the
Levenberg–Marquardt Method*, IEEE Transactions on Power Delivery 29(4),
1621–1629 (2014), DOI: 10.1109/TPWRD.2013.2284504. The implementation uses the
package vacuum permittivity, 8.8541878128e-12 F/m, and evaluates `s = j2πf`.
"""
function description(::Type{<:Formula{:vieira2026}}; compact::Bool = false)
    return compact ? "Vieira" :
           "Vieira complex Levenberg–Marquardt modal transformation (2026)"
end

function formulation_options(::FormulaMethod{<:Formula{:vieira2026}, typeof(decompose!)})
    return FormulationOptions((
        iteration = (convergence = 1e-13, max_iterations = 100),
        tracking = (predictor_tolerance = 0.25, eigenvalue_tolerance = 1e-5,
            cluster_tolerance = 1e-3, order_by_velocity = true)))
end

function formulation_options(::FormulaMethod{<:Formula{:vieira2026}, typeof(decompose!)},
        ::Val{:iteration}, defaults::NamedTuple, supplied::NamedTuple)
    isempty(setdiff(keys(supplied), keys(defaults))) ||
        throw(ArgumentError("unknown modal iteration controls"))
    options = merge(defaults, supplied)
    options.convergence isa Real && isfinite(options.convergence) &&
    options.convergence > 0 ||
        throw(ArgumentError("convergence must be finite and positive"))
    options.max_iterations isa Integer && !(options.max_iterations isa Bool) &&
    options.max_iterations > 0 ||
        throw(ArgumentError("max_iterations must be a positive integer"))
    return options
end

function formulation_options(::FormulaMethod{<:Formula{:vieira2026}, typeof(decompose!)},
        ::Val{:tracking}, defaults::NamedTuple, supplied::NamedTuple)
    isempty(setdiff(keys(supplied), keys(defaults))) ||
        throw(ArgumentError("unknown modal tracking controls"))
    options = merge(defaults, supplied)
    for key in (:predictor_tolerance, :eigenvalue_tolerance, :cluster_tolerance)
        value = getproperty(options, key)
        value isa Real && isfinite(value) && value > 0 ||
            throw(ArgumentError("$key must be finite and positive"))
    end
    options.order_by_velocity isa Bool ||
        throw(ArgumentError("order_by_velocity must be Bool"))
    return options
end

function levenberg_marquardt_workspace(::Val{:vieira2026}, ::Type{R}, n::Integer) where {R <:
                                                                                         Real}
    T = Complex{R}
    order = n + 1
    return (
        x = Vector{T}(undef, order), candidate = Vector{T}(undef, order),
        residual = Vector{T}(undef, order), candidate_residual = Vector{T}(undef, order),
        jacobian = Matrix{T}(undef, order, order),
        system = Matrix{T}(undef, 2order, order), rhs = Vector{T}(undef, 2order),
        column_norms = Vector{R}(undef, order), step = Vector{T}(undef, order),
        projection = Vector{T}(undef, order))
end

function initialize_buffers(
        ::Val{:vieira2026}, ::Type{T}, input, invariants, buffers) where {T <: Complex}
    n = invariants.n
    R = typeof(real(zero(T)))
    return merge(buffers,
        (
            normalized_shifted_eigenproblem = Matrix{T}(undef, n, n),
            spectral_factor = Matrix{T}(undef, n, n),
            prediction = (older_vectors = Matrix{T}(undef, n, n),
                older_values = Vector{T}(undef, n), vectors = Matrix{T}(undef, n, n),
                values = Vector{T}(undef, n)),
            least_squares = levenberg_marquardt_workspace(Val(:vieira2026), R, n),
            eigenpair_assignment = (
                cost = Matrix{R}(undef, n, n), assignment = Vector{Int}(undef, n),
                labels = Vector{Int}(undef, n), stack = Vector{Int}(undef, n),
                columns = Vector{Int}(undef, n), tracks = Vector{Int}(undef, n),
                cluster_assignment = Vector{Int}(undef, n), isolated = trues(n),
                system = Matrix{T}(undef, n, n), rhs = Matrix{T}(undef, n, n),
                coefficients = Matrix{T}(undef, n, n), projection = Matrix{T}(undef, n, n)),
            eigen_residual = Vector{T}(undef, n),
            mode_order = (indices = Vector{Int}(undef, n),
                residual = Vector{Union{Nothing, R}}(undef, n),
                iterations = Vector{Union{Nothing, Int}}(undef, n),
                converged = Vector{Union{Nothing, Bool}}(undef, n))))
end

# Truncated SVD minimum-norm solve; matrix is disposable factorization storage.
function minimum_norm!(solution, matrix::AbstractMatrix{T}, rhs, projection) where {T <:
                                                                                    Complex}
    factor = svd!(matrix)
    R = typeof(real(zero(T)))
    cutoff = eps(R) * max(size(matrix)...) * first(factor.S)
    mul!(projection, adjoint(factor.U), rhs)
    @inbounds for column in axes(projection, 2), row in axes(projection, 1)

        singular = factor.S[row]
        projection[row, column] = singular > cutoff ? projection[row, column] / singular :
                                  zero(T)
    end
    mul!(solution, adjoint(factor.Vt), projection)
    return solution
end

function levenberg_marquardt_step!(::Val{:vieira2026}, vector::AbstractVector{T},
        value::T, matrix, options, work) where {T <: Complex}
    n = length(vector)
    order = n + 1
    R = typeof(real(zero(T)))
    tolerance = convert(R, options.convergence)
    copyto!(@view(work.x[1:n]), vector)
    work.x[end] = value
    eigenpair_residual!(work.residual, work.x, matrix)
    cost = real(dot(work.residual, work.residual))
    damping = zero(R)
    iterations = 0
    for iteration in 1:options.max_iterations
        sqrt(cost) < tolerance && break
        iterations = iteration
        eigenpair_jacobian!(work.jacobian, work.x, matrix)
        for column in 1:order
            work.column_norms[column] = norm(@view(work.jacobian[:, column]))
        end
        candidate_cost = cost
        while true
            rows = iszero(damping) ? order : 2order
            copyto!(@view(work.system[1:order, :]), work.jacobian)
            @views work.rhs[1:order] .= -work.residual
            if !iszero(damping)
                fill!(@view(work.system[(order + 1):end, :]), zero(T))
                fill!(@view(work.rhs[(order + 1):end]), zero(T))
                for column in 1:order
                    work.system[order + column, column] = sqrt(damping) *
                                                          work.column_norms[column]
                end
            end
            minimum_norm!(work.step, @view(work.system[1:rows, :]),
                @view(work.rhs[1:rows]), work.projection)
            work.candidate .= work.x .+ work.step
            eigenpair_residual!(work.candidate_residual, work.candidate, matrix)
            candidate_cost = real(dot(work.candidate_residual, work.candidate_residual))
            (candidate_cost < cost || damping > R(1e8)) && break
            damping = iszero(damping) ? R(1e-6) : 10damping
        end
        candidate_cost < cost || break
        copyto!(work.x, work.candidate)
        copyto!(work.residual, work.candidate_residual)
        cost = candidate_cost
        damping = damping <= R(1e-6) ? zero(R) : damping / 10
    end
    copyto!(vector, @view(work.x[1:n]))
    return work.x[end], sqrt(cost) < tolerance, iterations
end

function refine_eigenpairs!(::Val{:vieira2026}, values, vectors, previous_vectors,
        prediction, eigensystem, options, work)
    n = length(values)
    R = eltype(work.cost)
    @inbounds for column in 1:n, row in 1:n

        work.cost[row, column] = abs(values[row] - eigensystem.values[column])
    end
    greedy_assignment!(work.assignment, work.cost)
    eigenvalue_limit = R(options.eigenvalue_tolerance) *
                       max(one(R), maximum(abs, eigensystem.values))
    any(
        row -> abs(values[row] - eigensystem.values[work.assignment[row]]) >
               eigenvalue_limit, 1:n) &&
        return false

    fill!(work.labels, 0)
    clusters = 0
    for seed in 1:n
        work.labels[seed] == 0 || continue
        clusters += 1
        work.labels[seed] = clusters
        pending = 1
        work.stack[pending] = seed
        while pending > 0
            row = work.stack[pending]
            pending -= 1
            for column in 1:n
                work.labels[column] == 0 || continue
                gap = abs(eigensystem.values[row] - eigensystem.values[column])
                limit = R(options.cluster_tolerance) *
                        max(abs(1 + eigensystem.values[row]),
                    abs(1 + eigensystem.values[column]), R(1e-3))
                if gap < limit
                    work.labels[column] = clusters
                    pending += 1
                    work.stack[pending] = column
                end
            end
        end
    end
    @inbounds for mode in 1:n
        values[mode] = eigensystem.values[work.assignment[mode]]
    end
    fill!(work.isolated, true)
    for cluster in 1:clusters
        count = 0
        tracks = 0
        for mode in 1:n
            if work.labels[mode] == cluster
                count += 1
                work.columns[count] = mode
            end
            if work.labels[work.assignment[mode]] == cluster
                tracks += 1
                work.tracks[tracks] = mode
            end
        end
        count == 1 && continue
        for column in 1:count
            track = work.tracks[column]
            work.isolated[track] = false
            copyto!(@view(work.system[:, column]), @view(previous_vectors[:, track]))
            copyto!(@view(work.rhs[:, column]), @view(eigensystem.vectors[:, work.columns[column]]))
        end
        coefficients = @view work.coefficients[1:count, 1:count]
        minimum_norm!(coefficients, @view(work.system[:, 1:count]),
            @view(work.rhs[:, 1:count]), @view(work.projection[1:count, 1:count]))
        cost = @view work.cost[1:count, 1:count]
        cost .= .-abs.(coefficients)
        assignment = @view work.cluster_assignment[1:count]
        greedy_assignment!(assignment, cost)
        for column in 1:count
            track = work.tracks[column]
            source = work.columns[assignment[column]]
            copyto!(@view(vectors[:, track]), @view(eigensystem.vectors[:, source]))
            values[track] = eigensystem.values[source]
        end
    end
    largest_change = zero(R)
    for mode in 1:n
        value_change = abs(values[mode] - prediction.values[mode]) /
                       max(abs(prediction.values[mode]), R(1e-3))
        largest_change = max(largest_change, value_change)
        if work.isolated[mode]
            change = zero(R)
            magnitude = zero(R)
            for row in 1:n
                change = max(change, abs(vectors[row, mode] -
                                         prediction.vectors[row, mode]))
                magnitude = max(magnitude, abs(prediction.vectors[row, mode]))
            end
            largest_change = max(largest_change, change / magnitude)
        end
    end
    return largest_change <= R(options.predictor_tolerance)
end

function decompose!(::Val{:vieira2026}, workspace::ModalAnalysisWorkspace,
        parameters::NamedTuple, options::FormulationOptions)
    work = workspace.buffers
    input = workspace.input
    diagnostics = workspace.diagnostics
    iteration = options.data.iteration
    tracking = options.data.tracking
    prediction = work.prediction
    assignment = work.eigenpair_assignment
    n, _, nf = size(workspace.Ti)
    T = eltype(workspace.Ti)
    R = typeof(real(zero(T)))
    unit = one(R)
    epsilon0 = unit * 88541878128 * (unit * 10)^(-22)
    mu0 = unit * 4 * (unit * π) * (unit * 10)^(-7)
    for frequency in 1:nf
        copyto!(work.Zslice, @view(input.Z[:, :, frequency]))
        copyto!(work.Yslice, @view(input.Y[:, :, frequency]))
        work.Zslice ./= input.root_scale
        work.Yslice ./= input.root_scale
        mul!(work.admittance_impedance_product, work.Yslice, work.Zslice)
        omega = 2 * (unit * π) * R(input.f[frequency])
        scale = -(omega^2) * epsilon0 * mu0
        matrix = work.normalized_shifted_eigenproblem
        matrix .= work.admittance_impedance_product ./ scale
        for mode in 1:n
            matrix[mode, mode] -= one(T)
        end
        copyto!(work.spectral_factor, matrix)
        eigensystem = eigen!(work.spectral_factor)
        for mode in 1:n
            vector = @view eigensystem.vectors[:, mode]
            vector ./= sqrt(sum(value -> value * value, vector))
        end
        if frequency == 1
            copyto!(work.eigenvalues, eigensystem.values)
            copyto!(work.eigenvectors, eigensystem.vectors)
        else
            coefficient = frequency == 2 ? zero(R) :
                          clamp(
                R(input.f[frequency] - input.f[frequency - 1]) /
                R(input.f[frequency - 1] - input.f[frequency - 2]), -one(R), one(R))
            if frequency == 2
                copyto!(prediction.vectors, work.previous_eigenvectors)
                copyto!(prediction.values, work.previous_eigenvalues)
            else
                prediction.vectors .= work.previous_eigenvectors .+
                                      coefficient .* (work.previous_eigenvectors .-
                                       prediction.older_vectors)
                prediction.values .= work.previous_eigenvalues .+
                                     coefficient .*
                                     (work.previous_eigenvalues .- prediction.older_values)
            end
            copyto!(work.eigenvectors, prediction.vectors)
            for mode in 1:n
                value, converged, iterations = levenberg_marquardt_step!(Val(:vieira2026),
                    @view(work.eigenvectors[:, mode]), prediction.values[mode], matrix,
                    iteration, work.least_squares)
                work.eigenvalues[mode] = value
                diagnostics.iterations[mode, frequency] = iterations
                diagnostics.converged[mode, frequency] = converged
            end
            matched = refine_eigenpairs!(Val(:vieira2026), work.eigenvalues,
                work.eigenvectors, work.previous_eigenvectors, prediction, eigensystem, tracking, assignment)
            if !matched
                # Preserve the reference's square correlation solve and greedy assignment.
                copyto!(assignment.system, work.previous_eigenvectors)
                copyto!(assignment.coefficients, eigensystem.vectors)
                ldiv!(lu!(assignment.system), assignment.coefficients)
                assignment.cost .= .-abs.(assignment.coefficients)
                greedy_assignment!(assignment.assignment, assignment.cost)
                for mode in 1:n
                    source = assignment.assignment[mode]
                    copyto!(@view(work.eigenvectors[:, mode]), @view(eigensystem.vectors[:, source]))
                    work.eigenvalues[mode] = eigensystem.values[source]
                end
                push!(diagnostics.fallback_frequencies, frequency)
            end
            if !matched || any(!, @view(diagnostics.converged[:, frequency]))
                push!(diagnostics.missed_frequencies, frequency)
            end
            for mode in 1:n
                real(dot(@view(work.previous_eigenvectors[:, mode]),
                    @view(work.eigenvectors[:, mode]))) < 0 &&
                    (@views work.eigenvectors[:, mode] .*= -one(T))
            end
            copyto!(prediction.older_vectors, work.previous_eigenvectors)
            copyto!(prediction.older_values, work.previous_eigenvalues)
        end
        copyto!(work.previous_eigenvectors, work.eigenvectors)
        copyto!(work.previous_eigenvalues, work.eigenvalues)
        work.propagation_eigenvalues .= (work.eigenvalues .+ one(T)) .* scale
        for mode in 1:n
            vector = @view workspace.Ti[:, mode, frequency]
            copyto!(vector, @view(work.eigenvectors[:, mode]))
            _unit!(vector) ||
                throw(ArgumentError("current eigenvector has zero or undefined norm"))
            mul!(work.voltage_vector, work.Zslice, vector)
            divisor = norm(work.voltage_vector)
            isfinite(divisor) && !iszero(divisor) ||
                throw(ArgumentError("voltage eigenvector has zero or undefined norm"))
            @views workspace.Tv[:, mode, frequency] .= work.voltage_vector ./ divisor
            eigenvalue = work.propagation_eigenvalues[mode]
            root = sqrt(eigenvalue)
            (real(root) < 0 || (iszero(real(root)) && imag(root) < 0)) && (root = -root)
            workspace.roots[mode, frequency] = root * input.root_scale
            mul!(work.eigen_residual, work.admittance_impedance_product, vector)
            work.eigen_residual .-= eigenvalue .* vector
            denominator = (norm(work.admittance_impedance_product, Inf) + abs(eigenvalue)) *
                          norm(vector, Inf)
            numerator = norm(work.eigen_residual, Inf)
            diagnostics.eigen_residual[mode, frequency] = iszero(denominator) ?
                                                          (iszero(numerator) ? zero(R) :
                                                           R(Inf)) : numerator / denominator
        end
    end
    if tracking.order_by_velocity
        order = work.mode_order
        sortperm!(order.indices, @view(workspace.roots[:, end]); by = value -> -imag(value))
        for frequency in 1:nf
            copyto!(work.eigenvectors, @view(workspace.Ti[:, :, frequency]))
            copyto!(work.previous_eigenvectors, @view(workspace.Tv[:, :, frequency]))
            copyto!(work.eigenvalues, @view(workspace.roots[:, frequency]))
            copyto!(order.residual, @view(diagnostics.eigen_residual[:, frequency]))
            copyto!(order.iterations, @view(diagnostics.iterations[:, frequency]))
            copyto!(order.converged, @view(diagnostics.converged[:, frequency]))
            for mode in 1:n
                source = order.indices[mode]
                copyto!(@view(workspace.Ti[:, mode, frequency]), @view(work.eigenvectors[:, source]))
                copyto!(@view(workspace.Tv[:, mode, frequency]), @view(work.previous_eigenvectors[:, source]))
                workspace.roots[mode, frequency] = work.eigenvalues[source]
                diagnostics.eigen_residual[mode, frequency] = order.residual[source]
                diagnostics.iterations[mode, frequency] = order.iterations[source]
                diagnostics.converged[mode, frequency] = order.converged[source]
            end
        end
    end
    isempty(diagnostics.missed_frequencies) ||
        @warn ":vieira2026 missed numerical targets" frequencies=copy(diagnostics.missed_frequencies) fallback_count=length(diagnostics.fallback_frequencies)
    return workspace
end

:vieira2026
