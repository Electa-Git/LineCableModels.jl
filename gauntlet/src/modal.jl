function _assignment(scores::AbstractMatrix{<:Real})
    n = LinearAlgebra.checksquare(scores)
    # Hungarian maximum-weight assignment in O(n^3).
    u = zeros(Float64, n + 1)
    v = zeros(Float64, n + 1)
    matched = zeros(Int, n + 1)
    path = zeros(Int, n + 1)
    for row in 1:n
        matched[1] = row
        column0 = 1
        minimums = fill(Inf, n + 1)
        used = falses(n + 1)
        while true
            used[column0] = true
            row0 = matched[column0]
            delta = Inf
            column1 = 1
            for column in 2:(n + 1)
                used[column] && continue
                cost = -Float64(scores[row0, column - 1]) - u[row0 + 1] - v[column]
                if cost < minimums[column]
                    minimums[column] = cost
                    path[column] = column0
                end
                if minimums[column] < delta
                    delta = minimums[column]
                    column1 = column
                end
            end
            for column in 1:(n + 1)
                if used[column]
                    u[matched[column] + 1] += delta
                    v[column] -= delta
                else
                    minimums[column] -= delta
                end
            end
            column0 = column1
            iszero(matched[column0]) && break
        end
        while true
            column1 = path[column0]
            matched[column0] = matched[column1]
            column0 = column1
            column0 == 1 && break
        end
    end
    assignment = zeros(Int, n)
    for column in 2:(n + 1)
        assignment[matched[column]] = column - 1
    end
    return assignment
end

function _correlations(actual::AbstractMatrix, reference::AbstractMatrix)
    n = LinearAlgebra.checksquare(actual)
    size(reference) == (n, n) || throw(DimensionMismatch(
        "modal transforms must have equal square dimensions",
    ))
    scores = Matrix{Float64}(undef, n, n)
    for target in 1:n, observed in 1:n

        left = @view reference[:, target]
        right = @view actual[:, observed]
        scores[target, observed] = abs(dot(left, right)) /
                                   max(norm(left) * norm(right), eps(Float64))
    end
    return scores
end

function _clusters(values::AbstractVector, tolerance::Tolerance)
    remaining = Set(eachindex(values))
    result = Vector{Vector{Int}}()
    while !isempty(remaining)
        seed = minimum(remaining)
        group = sort!([index
                       for index in remaining
                       if isapprox(
            values[index], values[seed]; rtol = tolerance.rtol, atol = tolerance.atol
        )])
        foreach(index -> delete!(remaining, index), group)
        push!(result, group)
    end
    return result
end

function _align_slice(actual, reference, assignment, groups)
    aligned = similar(actual)
    for group in groups
        observed = actual[:, assignment[group]]
        target = reference[:, group]
        if length(group) == 1
            phase = dot(@view(target[:, 1]), @view(observed[:, 1]))
            aligned[:, group[1]] .= @view(observed[:, 1]) .* exp(-im * angle(phase))
        else
            factor = svd(observed' * target)
            aligned[:, group] .= observed * (factor.U * factor.Vt)
        end
    end
    return aligned
end

function _modal_alignment(actual::Modes, reference::Modes, tolerance::Tolerance)
    actual.transform === nothing && throw(ArgumentError(
        "native modal transformation is unavailable",
    ))
    reference.transform === nothing && throw(ArgumentError(
        "reference modal transformation is unavailable",
    ))
    indices = _frequency_indices(actual.frequency, reference.frequency)
    observed = actual.transform[:, :, indices]
    target = reference.transform
    size(observed) == size(target) || throw(DimensionMismatch(
        "native and reference modal transforms have different dimensions",
    ))
    n, _, count = size(target)
    aligned = similar(observed)
    assignments = Vector{Vector{Int}}(undef, count)
    groups = Vector{Vector{Vector{Int}}}(undef, count)
    residual = zeros(Float64, count)
    for frequency in 1:count
        assignment = _assignment(_correlations(
            @view(observed[:, :, frequency]), @view(target[:, :, frequency])
        ))
        assignments[frequency] = assignment
        modal_values = reference.propagation === nothing ?
                       collect(1.0:n) : reference.propagation[:, frequency]
        groups[frequency] = _clusters(modal_values, tolerance)
        aligned[:, :, frequency] .= _align_slice(
            @view(observed[:, :, frequency]),
            @view(target[:, :, frequency]),
            assignment,
            groups[frequency]
        )
        residual[frequency] = norm(
            @view(aligned[:, :, frequency]) - @view(target[:, :, frequency])
        ) / max(norm(@view(target[:, :, frequency])), eps(Float64))
    end
    return (; aligned, assignments, groups, residual, indices)
end

function _native_modes(actual::EN.LineParameters, line_length::Real)
    impedance = Array(Z(actual))
    admittance = Array(Y(actual))
    n, _, count = size(impedance)
    transform = similar(impedance)
    propagation = Matrix{eltype(impedance)}(undef, n, count)
    for frequency in 1:count
        decomposition = eigen(
            @view(admittance[:, :, frequency]) *
            @view(impedance[:, :, frequency])
        )
        transform[:, :, frequency] .= decomposition.vectors
        propagation[:, frequency] .= decomposition.values
    end
    _, _, _, _, _, characteristic = EN.Transforms.modal_quantities(
        transform, impedance, admittance
    )
    modal_h = similar(propagation)
    phase_h = similar(transform)
    for frequency in 1:count
        modal_h[:, frequency] .= exp.(
            .-sqrt.(propagation[:, frequency]) .* line_length
        )
        matrix = @view transform[:, :, frequency]
        phase_h[:, :, frequency] .= matrix * Diagonal(modal_h[:, frequency]) * inv(matrix)
    end
    return Modes(
        collect(frequencies(actual));
        transform,
        propagation,
        characteristic,
        phase_propagation = phase_h,
        modal_propagation = modal_h
    )
end

_line_length(case::Case) = case.problem.system.line_length

function _aligned_parameters(actual::EN.LineParameters, ports)
    isempty(ports) && return actual
    order = getfield.(ports, :phase)
    sort(order) == collect(1:size(Z(actual), 1)) || throw(ArgumentError(
        "modal reference ports do not map one-to-one onto native phases",
    ))
    return EN.LineParameters(
        LineCableModels.PhaseDomain,
        Array(Z(actual))[order, order, :],
        Array(Y(actual))[order, order, :],
        collect(frequencies(actual));
        basis = basis(actual)
    )
end

function _modal_state(actual::EN.LineParameters, line_length::Real, ports = Port[])
    aligned = _aligned_parameters(actual, ports)
    return ModalState(aligned, _native_modes(aligned, line_length))
end

function compare_modes(actual, reference::Nothing, check, tolerance)
    return Comparison(
        check,
        Unavailable(),
        Tolerance(Float64(tolerance.rtol), Float64(tolerance.atol)),
        nothing,
        "modal reference is unavailable"
    )
end

function compare_modes(
        actual::EN.LineParameters,
        reference::Modes,
        check::ModalCheck,
        tolerance
)
    return compare_modes(_modal_state(actual, one(eltype(frequencies(actual)))),
        reference, check, tolerance)
end

function _modal_metrics(values, target, frequency, tolerance)
    observed = ndims(values) == 2 ? reshape(values, size(values, 1), 1, :) : values
    expected = ndims(target) == 2 ? reshape(target, size(target, 1), 1, :) : target
    size(observed) == size(expected) || throw(DimensionMismatch(
        "native and reference modal channels must have equal dimensions",
    ))
    difference = observed .- expected
    absolute_error = abs.(difference)
    denominator = max.(abs.(expected), tolerance.atol, eps(Float64))
    relative_error = absolute_error ./ denominator
    worst = argmax(relative_error)
    row, column, frequency_index = Tuple(worst)
    frobenius = Float64[norm(@view difference[:, :, index]) /
                        max(norm(@view expected[:, :, index]), tolerance.atol, eps(Float64))
                        for index in axes(expected, 3)]
    return Metrics(
        Float64(maximum(absolute_error; init = 0.0)),
        Float64(maximum(relative_error; init = 0.0)),
        frobenius,
        Float64(frequency[frequency_index]),
        (row, column),
        0.0,
        0.0,
        0.0,
        0.0
    )
end

function _matched_values(values, alignment)
    source = values[:, alignment.indices]
    result = similar(source)
    for frequency in axes(result, 2)
        result[:, frequency] .= source[alignment.assignments[frequency], frequency]
    end
    return result
end

function _modal_comparison(check, metrics, tolerance, detail)
    return Comparison(
        check,
        _verdict(metrics, tolerance),
        Tolerance(Float64(tolerance.rtol), Float64(tolerance.atol)),
        metrics,
        detail
    )
end

function _fit_comparison(check, metrics, tolerance, detail)
    verdict = _verdict(metrics, tolerance)
    return Comparison(
        check,
        verdict isa Pass ? verdict : Diagnostic(),
        Tolerance(Float64(tolerance.rtol), Float64(tolerance.atol)),
        metrics,
        verdict isa Pass ? detail :
        detail *
        "; imported coefficients disagree with the detailed fitted channel"
    )
end

function compare_modes(
        state::ModalState,
        reference::Modes,
        check::ModalCheck{:transform},
        tolerance
)
    native = state.modes
    alignment = _modal_alignment(native, reference, tolerance)
    difference = alignment.aligned .- reference.transform
    max_rel = maximum(alignment.residual)
    worst = argmax(alignment.residual)
    metrics = Metrics(
        maximum(abs, difference), max_rel, alignment.residual,
        reference.frequency[worst], (1, 1), 0.0, 0.0, 0.0, max_rel
    )
    return _modal_comparison(
        check,
        metrics,
        tolerance,
        "modal columns matched by O(n³) assignment, phase alignment, and " *
        "degenerate-subspace Procrustes alignment"
    )
end

function compare_modes(
        state::ModalState,
        reference::Modes,
        check::ModalCheck{:propagation},
        tolerance
)
    reference.propagation === nothing && return compare_modes(
        state, nothing, check, tolerance
    )
    native = state.modes
    alignment = _modal_alignment(native, reference, tolerance)
    metrics = _modal_metrics(
        _matched_values(native.propagation, alignment),
        reference.propagation,
        reference.frequency,
        tolerance
    )
    return _modal_comparison(
        check, metrics, tolerance, "compared propagation quantities after modal assignment"
    )
end

function _phase_modal_comparison(state, reference, check, tolerance, field::Symbol)
    target = getfield(reference, field)
    target === nothing && return compare_modes(state, nothing, check, tolerance)
    native = state.modes
    indices = _frequency_indices(native.frequency, reference.frequency)
    observed = getfield(native, field)
    observed = ndims(observed) == 3 ? observed[:, :, indices] : observed[:, indices]
    metrics = _modal_metrics(observed, target, reference.frequency, tolerance)
    return _modal_comparison(
        check,
        metrics,
        tolerance,
        "compared calculated reference $(replace(String(field), '_' => ' '))"
    )
end

function compare_modes(state::ModalState, reference::Modes,
        check::ModalCheck{:characteristic}, tolerance)
    _phase_modal_comparison(state, reference, check, tolerance, :characteristic)
end
function compare_modes(state::ModalState, reference::Modes,
        check::ModalCheck{:phase_propagation}, tolerance)
    _phase_modal_comparison(state, reference, check, tolerance, :phase_propagation)
end

function compare_modes(
        state::ModalState,
        reference::Modes,
        check::ModalCheck{:modal_propagation},
        tolerance
)
    reference.modal_propagation === nothing && return compare_modes(
        state, nothing, check, tolerance
    )
    native = state.modes
    alignment = _modal_alignment(native, reference, tolerance)
    observed = _matched_values(native.modal_propagation, alignment)
    observed = observed[axes(reference.modal_propagation, 1), :]
    metrics = _modal_metrics(
        observed,
        reference.modal_propagation,
        reference.frequency,
        tolerance
    )
    return _modal_comparison(
        check,
        metrics,
        tolerance,
        "compared modal propagation functions after modal assignment"
    )
end

function compare_modes(
        state::ModalState,
        reference::Modes,
        check::ModalCheck{:reconstruction},
        tolerance
)
    reference.transform === nothing && return compare_modes(
        state, nothing, check, tolerance
    )
    native = state.modes
    alignment = _modal_alignment(native, reference, tolerance)
    residual = Float64[]
    for (offset, index) in enumerate(alignment.indices)
        transform = @view alignment.aligned[:, :, offset]
        product = @view(Y(state.parameters)[:, :, index]) *
                  @view(Z(state.parameters)[:, :, index])
        diagonal = Diagonal(diag(inv(transform) * product * transform))
        push!(residual, norm(product * transform - transform * diagonal) /
                        max(norm(product), eps(Float64)))
    end
    worst = argmax(residual)
    maximum_residual = maximum(residual)
    metrics = Metrics(
        maximum_residual, maximum_residual, residual,
        reference.frequency[worst], (1, 1), 0.0, 0.0, 0.0, maximum_residual
    )
    return _modal_comparison(
        check,
        metrics,
        tolerance,
        "checked phase-domain reconstruction after modal alignment"
    )
end

function compare_fit(fit::Nothing, modes, check, tolerance)
    return Comparison(
        check,
        Unavailable(),
        Tolerance(Float64(tolerance.rtol), Float64(tolerance.atol)),
        nothing,
        "vector-fit reference is unavailable"
    )
end

function compare_fit(fit::Fit, modes::Nothing, check, tolerance)
    return Comparison(
        check,
        Unavailable(),
        Tolerance(Float64(tolerance.rtol), Float64(tolerance.atol)),
        nothing,
        "detailed fitted channels are unavailable"
    )
end

function compare_fit(fit::Fit, modes::Modes, check::FitCheck{:Yc}, tolerance)
    modes.fitted_characteristic === nothing && return compare_fit(
        fit, nothing, check, tolerance
    )
    evaluated = evaluate(fit, modes.frequency)
    metrics = _metrics(
        evaluated.characteristic,
        modes.fitted_characteristic,
        modes.frequency,
        tolerance
    )
    return _fit_comparison(
        check, metrics, tolerance, "evaluated imported characteristic-admittance fit"
    )
end

function compare_fit(fit::Fit, modes::Modes, check::FitCheck{:H}, tolerance)
    modes.fitted_phase_propagation === nothing && return compare_fit(
        fit, nothing, check, tolerance
    )
    evaluated = evaluate(fit, modes.frequency)
    metrics = _metrics(
        evaluated.propagation,
        modes.fitted_phase_propagation,
        modes.frequency,
        tolerance
    )
    return _fit_comparison(
        check, metrics, tolerance, "evaluated imported propagation fit"
    )
end
