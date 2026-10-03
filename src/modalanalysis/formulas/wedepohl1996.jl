"""
$(TYPEDSIGNATURES)

Track current eigenpairs by Newton–Raphson refinement of the normalized
admittance-impedance product:

```math
A = YZ/\\|YZ\\|_2,\\qquad
r(t,\\lambda) = \\begin{bmatrix}(A-\\lambda I)t\\\\t^Tt-1\\end{bmatrix},\\qquad
J = \\begin{bmatrix}A-\\lambda I&-t\\\\2t^T&0\\end{bmatrix}.
```

Here `Z` and `Y` are per-length matrices in \\[Ω/m\\] and \\[S/m\\]. Each
eigenpair starts from the previous stored frequency. The seed is ordered by
decreasing attenuation. Failed iterations or duplicate eigenpairs use a direct
eigendecomposition at the same frequency, greedily matched by column overlap.
The calculation uses the existing physical evaluation and frequency samples.

`iteration.convergence=1e-9` bounds the largest absolute Newton correction in
the normalized problem. `iteration.max_iterations=60` limits each eigenpair.
Misses are recorded and warned about, without discarding the fallback result.

Adapted from `eig_newton` in the supplied UniversalLineModel `modal.jl`, based on
L. M. Wedepohl, H. V. Nguyen and G. D. Irwin, *Frequency-dependent transformation
matrices for untransposed transmission lines using Newton-Raphson method*,
IEEE Transactions on Power Systems 11(3), 1538–1546 (1996),
DOI: 10.1109/59.535695. The selection identifier is `:wedepohl1996`.
"""
function description(::Type{<:Formula{:wedepohl1996}}; compact::Bool = false)
    compact ? "Wedepohl" :
    "Wedepohl–Nguyen–Irwin Newton–Raphson modal transformation (1996)"
end

function formulation_options(::FormulaMethod{<:Formula{:wedepohl1996}, typeof(decompose!)})
    FormulationOptions((iteration = (convergence = 1e-9, max_iterations = 60),))
end

function formulation_options(::FormulaMethod{<:Formula{:wedepohl1996}, typeof(decompose!)},
        ::Val{:iteration}, defaults::NamedTuple, supplied::NamedTuple)
    isempty(setdiff(keys(supplied), keys(defaults))) ||
        throw(ArgumentError("unknown modal iteration controls"))
    options=merge(defaults, supplied)
    options.convergence isa Real && isfinite(options.convergence) &&
    options.convergence>0 ||
        throw(ArgumentError("convergence must be finite and positive"))
    options.max_iterations isa Integer && !(options.max_iterations isa Bool) &&
    options.max_iterations>0 ||
        throw(ArgumentError("max_iterations must be a positive integer"))
    return options
end

function initialize_buffers(
        ::Val{:wedepohl1996}, ::Type{T}, input, invariants, buffers) where {T <: Complex}
    n=invariants.n
    R=typeof(real(zero(T)))
    return merge(buffers,
        (
            normalized_eigenproblem = Matrix{T}(undef, n, n),
            spectral_factor = Matrix{T}(undef, n, n),
            newton = (x = Vector{T}(undef, n+1), residual = Vector{T}(undef, n+1),
                jacobian = Matrix{T}(undef, n+1, n+1)),
            eigenpair_assignment = (
                cost = Matrix{R}(undef, n, n), indices = Vector{Int}(undef, n)),
            eigen_residual = Vector{T}(undef, n)))
end

function newton_eigenpair!(
        vector::AbstractVector{T}, value::T, matrix, options, work) where {T <: Complex}
    n=length(vector)
    R=typeof(real(zero(T)))
    _unit!(vector) || return value, false, 0
    copyto!(@view(work.x[1:n]), vector)
    work.x[end]=value
    correction=R(Inf)
    iterations=0
    for iteration in 1:options.max_iterations
        iterations=iteration
        eigenpair_jacobian!(work.jacobian, work.x, matrix)
        eigenpair_residual!(work.residual, work.x, matrix)
        factor=lu!(work.jacobian; check = false)
        issuccess(factor) || break
        ldiv!(factor, work.residual)
        correction=maximum(abs, work.residual)
        work.x .-= work.residual
        all(isfinite, work.x) || break
        correction<=R(options.convergence) && break
    end
    copyto!(vector, @view(work.x[1:n]))
    valid=all(isfinite, work.x) && _unit!(vector)
    return work.x[end], valid && correction<=R(options.convergence), iterations
end

function decompose!(::Val{:wedepohl1996}, workspace::ModalAnalysisWorkspace,
        parameters::NamedTuple, options::FormulationOptions)
    work=workspace.buffers
    input=workspace.input
    diagnostics=workspace.diagnostics
    n, _, nf=size(workspace.Ti)
    T=eltype(workspace.Ti)
    R=typeof(real(zero(T)))
    for frequency in 1:nf
        copyto!(work.Zslice, @view(input.Z[:, :, frequency]))
        copyto!(work.Yslice, @view(input.Y[:, :, frequency]))
        work.Zslice ./= input.root_scale
        work.Yslice ./= input.root_scale
        matrix=work.admittance_impedance_product
        mul!(matrix, work.Yslice, work.Zslice)
        use_newton_result=frequency>1
        if frequency>1
            copyto!(work.spectral_factor, matrix)
            scale=first(svdvals!(work.spectral_factor))
            iszero(scale) && (scale=one(R))
            work.normalized_eigenproblem .= matrix ./ scale
            copyto!(work.eigenvectors, work.previous_eigenvectors)
            for mode in 1:n
                value, converged, iterations=newton_eigenpair!(
                    @view(work.eigenvectors[:, mode]), work.previous_eigenvalues[mode]/scale,
                    work.normalized_eigenproblem, options.data.iteration, work.newton)
                work.eigenvalues[mode]=value*scale
                diagnostics.iterations[mode, frequency]=iterations
                diagnostics.converged[mode, frequency]=converged
                use_newton_result &= converged
            end
            if use_newton_result
                for first_mode in 1:n, second_mode in (first_mode + 1):n

                    overlap=abs(dot(@view(work.eigenvectors[:, first_mode]),
                        @view(work.eigenvectors[:, second_mode])))
                    gap=abs(work.eigenvalues[first_mode]-work.eigenvalues[second_mode])
                    if one(R)-overlap<R(1e-8) &&
                       gap<=R(1e-8)*abs(work.eigenvalues[first_mode])
                        use_newton_result=false
                    end
                end
            end
        end
        if !use_newton_result
            copyto!(work.spectral_factor, matrix)
            eigensystem=eigen!(work.spectral_factor)
            assignment=work.eigenpair_assignment
            if frequency==1
                sortperm!(assignment.indices, eigensystem.values; by = value->-real(sqrt(value)))
            else
                for column in 1:n, row in 1:n

                    assignment.cost[row, column]=-abs(dot(
                        @view(work.previous_eigenvectors[:, row]),
                        @view(eigensystem.vectors[:, column])))
                end
                greedy_assignment!(assignment.indices, assignment.cost)
                push!(diagnostics.missed_frequencies, frequency)
                push!(diagnostics.fallback_frequencies, frequency)
            end
            for mode in 1:n
                source=assignment.indices[mode]
                work.eigenvalues[mode]=eigensystem.values[source]
                copyto!(@view(work.eigenvectors[:, mode]), @view(eigensystem.vectors[:, source]))
            end
        end
        for mode in 1:n
            vector=@view(work.eigenvectors[:, mode])
            _unit!(vector) ||
                throw(ArgumentError("current eigenvector has zero or undefined norm"))
            if frequency>1 &&
               real(dot(@view(work.previous_eigenvectors[:, mode]), vector))<0
                vector .*= -one(T)
            end
            copyto!(@view(workspace.Ti[:, mode, frequency]), vector)
            mul!(work.voltage_vector, work.Zslice, vector)
            divisor=norm(work.voltage_vector)
            isfinite(divisor) && !iszero(divisor) ||
                throw(ArgumentError("voltage eigenvector has zero or undefined norm"))
            @views workspace.Tv[:, mode, frequency] .= work.voltage_vector ./ divisor
            value=work.eigenvalues[mode]
            root=sqrt(value)
            (real(root)<0 || (iszero(real(root)) && imag(root)<0)) && (root=-root)
            workspace.roots[mode, frequency]=root*input.root_scale
            mul!(work.eigen_residual, matrix, vector)
            work.eigen_residual .-= value .* vector
            denominator=(norm(matrix, Inf)+abs(value))*norm(vector, Inf)
            numerator=norm(work.eigen_residual, Inf)
            diagnostics.eigen_residual[mode, frequency]=iszero(denominator) ?
                                                        (iszero(numerator) ? zero(R) :
                                                         R(Inf)) : numerator/denominator
        end
        copyto!(work.previous_eigenvectors, work.eigenvectors)
        copyto!(work.previous_eigenvalues, work.eigenvalues)
    end
    isempty(diagnostics.missed_frequencies) ||
        @warn ":wedepohl1996 missed numerical targets" frequencies=copy(diagnostics.missed_frequencies) fallback_count=length(diagnostics.fallback_frequencies)
    return workspace
end

:wedepohl1996
