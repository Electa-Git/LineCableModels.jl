_quantity(::MatrixCheck{:Z}, result) = Array(Z(result))
_quantity(::MatrixCheck{:Y}, result) = Array(Y(result))
_quantity(::ModalCheck{:Z}, result) = Array(Z(result))
_quantity(::ModalCheck{:Y}, result) = Array(Y(result))

function _frequency_indices(actual, reference)
    indices = Int[]
    for frequency in reference
        index = findfirst(value -> isapprox(value, frequency; rtol = 1e-7, atol = 0), actual)
        index === nothing && throw(DimensionMismatch(
            "actual result has no sample at reference frequency $frequency Hz",
        ))
        push!(indices, index)
    end
    return indices
end

function _property_residual(values, operation)
    result = 0.0
    for index in axes(values, 3)
        result = max(result, operation(@view values[:, :, index]))
    end
    return result
end

_symmetry(matrix) = norm(matrix - transpose(matrix)) / max(norm(matrix), eps(Float64))
_reciprocity(matrix) = _symmetry(matrix)
function _positive_real(matrix)
    hermitian = Hermitian((matrix + matrix') / 2)
    negative = max(0.0, -minimum(real.(eigvals(hermitian))))
    return negative / max(norm(hermitian), eps(Float64))
end

function _metrics(actual, expected, frequency, tolerance::Tolerance)
    size(actual) == size(expected) || throw(DimensionMismatch(
        "actual and reference matrix tensors must have equal dimensions",
    ))
    difference = actual .- expected
    absolute_error = abs.(difference)
    denominator = max.(abs.(expected), tolerance.atol)
    relative_error = absolute_error ./ max.(denominator, eps(Float64))
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
        _property_residual(actual, _symmetry),
        _property_residual(actual, _reciprocity),
        _property_residual(actual, _positive_real),
        0.0
    )
end

function _verdict(metrics::Metrics, tolerance::Tolerance)
    acceptable = metrics.max_abs <= tolerance.atol || metrics.max_rel <= tolerance.rtol
    return acceptable ? Pass() : Fail()
end

"""Compare one native result matrix against a reference result."""
function compare(
        actual::EN.LineParameters,
        expected::EN.LineParameters,
        check::MatrixCheck,
        tolerance::Tolerance;
        ports = Port[]
)
    basis(actual) == basis(expected) || throw(ArgumentError(
        "actual and reference bases differ: $(basis(actual)) and $(basis(expected))",
    ))
    indices = _frequency_indices(frequencies(actual), frequencies(expected))
    values = _quantity(check, actual)
    order = if isempty(ports)
        collect(axes(values, 1))
    else
        phase = getfield.(ports, :phase)
        length(phase) == size(values, 1) || throw(DimensionMismatch(
            "$(length(phase)) reference ports cannot align a $(size(values, 1))-port result",
        ))
        sort(phase) == collect(1:size(values, 1)) || throw(ArgumentError(
            "reference ports do not map one-to-one onto native phase identifiers",
        ))
        phase
    end
    observed = values[order, order, indices]
    target = _quantity(check, expected)
    metrics = _metrics(observed, target, frequencies(expected), tolerance)
    return Comparison(
        check,
        _verdict(metrics, tolerance),
        Tolerance(Float64(tolerance.rtol), Float64(tolerance.atol)),
        metrics,
        "compared $(size(target, 1))×$(size(target, 2)) matrices at " *
        "$(size(target, 3)) frequencies"
    )
end

function _sequence_matrix(transform, ports, n)
    transform === nothing && throw(ArgumentError(
        "the reference sequence transform is unavailable",
    ))
    q = LinearAlgebra.checksquare(transform)
    length(ports) == n || throw(DimensionMismatch(
        "$(length(ports)) explicit ports cannot describe a $n-port sequence matrix",
    ))
    groups = Vector{Vector{Int}}()
    names = String[]
    for (index, port) in enumerate(ports)
        position = findfirst(==(port.component), names)
        if position === nothing
            push!(names, port.component)
            push!(groups, [index])
        else
            push!(groups[position], index)
        end
    end
    chunks = Vector{Pair{Vector{Int}, Bool}}()
    for (name, group) in zip(names, groups)
        transforms = !startswith(name, "ground:") && length(group) % q == 0
        if transforms
            for first_index in 1:q:length(group)
                columns = group[first_index:(first_index + q - 1)]
                push!(chunks, columns => true)
            end
        else
            for column in group
                push!(chunks, [column] => false)
            end
        end
    end
    matrix = zeros(ComplexF64, n, n)
    row = 1
    for (columns, transforms) in chunks
        if transforms
            rows = row:(row + q - 1)
            matrix[rows, columns] .= transform
            row += q
        else
            matrix[row, only(columns)] = 1
            row += 1
        end
    end
    return matrix
end

function compare(
        actual::EN.LineParameters,
        expected::EN.LineParameters,
        check::ModalCheck{Q},
        tolerance::Tolerance;
        transform = nothing,
        ports = Port[]
) where {Q}
    indices = _frequency_indices(frequencies(actual), frequencies(expected))
    observed = _quantity(MatrixCheck{Q}(), actual)[:, :, indices]
    matrix = _sequence_matrix(transform, ports, size(observed, 1))
    transformed = similar(observed)
    for index in axes(observed, 3)
        transformed[:, :, index] .= matrix * observed[:, :, index] * matrix'
    end
    target = _quantity(MatrixCheck{Q}(), expected)
    metrics = _metrics(transformed, target, frequencies(expected), tolerance)
    return Comparison(
        check,
        _verdict(metrics, tolerance),
        Tolerance(Float64(tolerance.rtol), Float64(tolerance.atol)),
        metrics,
        "applied the recorded reference sequence transform to aligned port groups"
    )
end

function compare(
        actual::EN.LineParameters,
        ::Nothing,
        check::Check,
        tolerance::Tolerance
)
    return Comparison(
        check,
        Unavailable(),
        Tolerance(Float64(tolerance.rtol), Float64(tolerance.atol)),
        nothing,
        "reference channel is unavailable"
    )
end

function compare(
        actual::EN.LineParameters,
        check::PhysicalCheck{Q},
        tolerance::Tolerance
) where {Q}
    values = Q === :positive_real ? (Array(Z(actual)), Array(Y(actual))) :
             (Array(Z(actual)), Array(Y(actual)))
    operation = Q === :symmetry ? _symmetry :
                Q === :reciprocity ? _reciprocity : _positive_real
    residual = maximum(_property_residual(value, operation) for value in values)
    metrics = Metrics(
        residual, residual, [residual], first(frequencies(actual)), (1, 1),
        Q === :symmetry ? residual : 0.0,
        Q === :reciprocity ? residual : 0.0,
        Q === :positive_real ? residual : 0.0,
        0.0
    )
    return Comparison(
        check,
        residual <= max(tolerance.atol, tolerance.rtol) ? Pass() : Fail(),
        Tolerance(Float64(tolerance.rtol), Float64(tolerance.atol)),
        metrics,
        "$(Q) residual"
    )
end

function _diagnostic(comparison::Comparison)
    return Comparison(
        comparison.check,
        comparison.verdict isa Unavailable ? comparison.verdict : Diagnostic(),
        comparison.tolerance,
        comparison.metrics,
        comparison.detail
    )
end

const _TOLERANCE_CONFIG = TOML.parsefile(joinpath(
    _GAUNTLET_ROOT, "config", "tolerances.toml"
))

_tolerance_key(::MatrixCheck{Q}) where {Q} = ("matrix", String(Q))
_tolerance_key(::ModalCheck{:Z}) = ("sequence", "Z")
_tolerance_key(::ModalCheck{:Y}) = ("sequence", "Y")
_tolerance_key(::ModalCheck{Q}) where {Q} = ("modal", String(Q))
_tolerance_key(::FitCheck{Q}) where {Q} = ("fit", String(Q))
_tolerance_key(::PhysicalCheck) = ("physical", nothing)

function _check_key(check::Check)
    join(
        filter(value -> value !== nothing && !isempty(value), _tolerance_key(check)), "."
    )
end

function _tolerance(check::Check)
    section, quantity = _tolerance_key(check)
    values = _TOLERANCE_CONFIG[section]
    quantity === nothing || (values = values[quantity])
    return Tolerance(Float64(values["rtol"]), Float64(values["atol"]))
end

function _tolerance(check::FitCheck, case::Case)
    section, quantity = _tolerance_key(check)
    values = _TOLERANCE_CONFIG[section][quantity]
    family = lowercase(String(nameof(typeof(case.family))))
    selected = get(values, family, values)
    return Tolerance(Float64(selected["rtol"]), Float64(selected["atol"]))
end

_tolerance(check::Check, ::Case) = _tolerance(check)

function _tolerance(check::Check, overrides::AbstractDict)
    key = _check_key(check)
    return get(overrides, key, _tolerance(check))
end

function _case_comparison(case::Case, actual, check::MatrixCheck{:Z})
    return compare(
        actual,
        case.reference.phase,
        check,
        _tolerance(check);
        ports = case.reference.ports
    )
end

function _case_comparison(case::Case, actual, check::MatrixCheck{:Y})
    return compare(
        actual,
        case.reference.phase,
        check,
        _tolerance(check);
        ports = case.reference.ports
    )
end

function _case_comparison(case::Case, actual, check::ModalCheck{:Z})
    return compare(
        actual, case.reference.sequence, check, _tolerance(check);
        transform = case.reference.sequence_transform,
        ports = case.reference.ports
    )
end

function _case_comparison(case::Case, actual, check::ModalCheck{:Y})
    return compare(
        actual, case.reference.sequence, check, _tolerance(check);
        transform = case.reference.sequence_transform,
        ports = case.reference.ports
    )
end

function _case_comparison(case::Case, actual, check::PhysicalCheck)
    return compare(actual, check, _tolerance(check))
end

function _case_comparison(case::Case, _, check::PhysicalCheck{:terminal})
    terminal = case.reference.terminal
    terminal === nothing && return Comparison(
        check,
        Unavailable(),
        _tolerance(check),
        nothing,
        "terminal-response reference is unavailable"
    )
    finite = all(values -> values === nothing || all(isfinite, values),
        (terminal.open, terminal.short))
    return Comparison(
        check,
        finite ? Pass() : Fail(),
        _tolerance(check),
        nothing,
        "validated imported open/short channels; empty source outputs: " *
        (isempty(terminal.empty) ? "none" : join(terminal.empty, ", "))
    )
end

function _case_comparison(case::Case, actual, check::ModalCheck)
    return compare_modes(actual, case.reference.modes, check, _tolerance(check))
end

function _case_comparison(case::Case, actual, check::ModalCheck, state::ModalState)
    return compare_modes(state, case.reference.modes, check, _tolerance(check))
end

function _case_comparison(case::Case, actual, check::Check, _)
    _case_comparison(case, actual, check)
end

function _case_comparison(case::Case, actual, check::ModalCheck{:Z}, ::ModalState)
    _case_comparison(case, actual, check)
end
function _case_comparison(case::Case, actual, check::ModalCheck{:Y}, ::ModalState)
    _case_comparison(case, actual, check)
end

function _case_comparison(case::Case, _, check::FitCheck)
    return compare_fit(
        case.reference.fit, case.reference.modes, check, _tolerance(check, case)
    )
end

_native_check(::MatrixCheck) = true
_native_check(::ModalCheck) = true
_native_check(::PhysicalCheck{:terminal}) = false
_native_check(::PhysicalCheck) = true
_native_check(::FitCheck) = false

function _diagnostic_policy(::ExactOnly, case::Case, check::Check)
    case.fidelity isa Approximate && _native_check(check)
end
_diagnostic_policy(::AllCases, ::Case, ::Check) = false

function _execution_failure(case::Case, error)
    check = PhysicalCheck{:execution}()
    comparison = Comparison(
        check,
        Fail(),
        _tolerance(check),
        nothing,
        sprint(showerror, error)
    )
    return Trial(case, nothing, Comparison[comparison], nothing)
end

"""Run one case through the native solver and retain complete diagnostics."""
function gauntlet(
        case::Case;
        performance::Bool = false,
        checks = case.checks,
        policy::Policy = ExactOnly(),
        tolerances = Dict{String, Tolerance{Float64}}()
)
    if case.fidelity isa Rejected
        comparison = Comparison(
            PhysicalCheck{:reference}(),
            ReferenceRejected(),
            Tolerance(0.0, 0.0),
            nothing,
            "$(case.provenance.datasource) rejected this source case; no Julia failure is inferred"
        )
        return Trial(case, nothing, Comparison[comparison], nothing)
    end
    actual = try
        compute!(case.problem, case.formulation)
    catch error
        error isa InterruptException && rethrow()
        return _execution_failure(case, error)
    end
    comparisons = Comparison[]
    modal_state = any(
        check -> check isa ModalCheck &&
                 !(check isa Union{ModalCheck{:Z}, ModalCheck{:Y}}),
        checks) ?
                  _modal_state(actual, _line_length(case), case.reference.ports) : nothing
    for check in checks
        comparison = try
            default = _tolerance(check)
            override = get(tolerances, _check_key(check), default)
            # Comparison methods obtain the frozen default internally. A suite
            # override replaces it by rerunning the public dispatch where needed.
            result = _case_comparison(case, actual, check, modal_state)
            override == default ? result :
            _with_tolerance(case, actual, check, override, modal_state)
        catch error
            error isa InterruptException && rethrow()
            push!(comparisons,
                Comparison(
                    check,
                    Fail(),
                    _tolerance(check, tolerances),
                    nothing,
                    sprint(showerror, error)
                ))
            continue
        end
        push!(comparisons, _diagnostic_policy(policy, case, check) ?
                           _diagnostic(comparison) : comparison)
    end
    performance_record = performance ? benchmark(case) : nothing
    performance_record === nothing || push!(
        comparisons, _performance_comparison(performance_record)
    )
    return Trial(case, actual, comparisons, performance_record)
end

function _with_tolerance(case::Case, actual, check::MatrixCheck, tolerance)
    return compare(
        actual,
        case.reference.phase,
        check,
        tolerance;
        ports = case.reference.ports
    )
end

function _with_tolerance(case::Case, actual, check::MatrixCheck, tolerance, _)
    _with_tolerance(case, actual, check, tolerance)
end

function _with_tolerance(case::Case, actual, check::ModalCheck, tolerance)
    return check isa Union{ModalCheck{:Z}, ModalCheck{:Y}} ?
           compare(
        actual,
        case.reference.sequence,
        check,
        tolerance;
        transform = case.reference.sequence_transform,
        ports = case.reference.ports
    ) : compare_modes(actual, case.reference.modes, check, tolerance)
end

function _with_tolerance(
        case::Case,
        actual,
        check::ModalCheck,
        tolerance,
        state::Union{Nothing, ModalState}
)
    check isa Union{ModalCheck{:Z}, ModalCheck{:Y}} && return _with_tolerance(
        case, actual, check, tolerance
    )
    state === nothing && return compare_modes(
        actual, case.reference.modes, check, tolerance
    )
    return compare_modes(state, case.reference.modes, check, tolerance)
end

function _with_tolerance(case::Case, actual, check::FitCheck, tolerance)
    compare_fit(case.reference.fit, case.reference.modes, check, tolerance)
end
function _with_tolerance(::Case, actual, check::PhysicalCheck, tolerance)
    compare(actual, check, tolerance)
end
function _with_tolerance(case::Case, actual, check::PhysicalCheck{:terminal}, tolerance)
    _case_comparison(case, actual, check)
end

function _with_tolerance(case::Case, actual, check::FitCheck, tolerance, _)
    _with_tolerance(case, actual, check, tolerance)
end
function _with_tolerance(case::Case, actual, check::PhysicalCheck, tolerance, _)
    _with_tolerance(case, actual, check, tolerance)
end

function _port_tensor(trial::Trial, quantity::Symbol)
    values = quantity === :Z ? Array(Z(trial.actual)) : Array(Y(trial.actual))
    ports = trial.case.reference.ports
    phase = getfield.(ports, :phase)
    sort(phase) == collect(1:size(values, 1)) || throw(ArgumentError(
        "case $(trial.case.id) has no one-to-one native port map",
    ))
    return values[phase, phase, :]
end

function _kron(values, keep, quantity::Symbol)
    eliminate = setdiff(axes(values, 1), keep)
    isempty(eliminate) && return values[keep, keep, :]
    result = similar(values, length(keep), length(keep), size(values, 3))
    for frequency in axes(values, 3)
        matrix = @view values[:, :, frequency]
        if quantity === :Y
            potential = inv(matrix)
            reduced = potential[keep, keep] -
                      potential[keep, eliminate] *
                      (potential[eliminate, eliminate] \
                       potential[eliminate, keep])
            result[:, :, frequency] .= inv(reduced)
        else
            result[:, :, frequency] .= matrix[keep, keep] -
                                       matrix[keep, eliminate] *
                                       (matrix[eliminate, eliminate] \
                                        matrix[eliminate, keep])
        end
    end
    return result
end

function _kron_comparison(retained::Trial, reduced::Trial, quantity::Symbol)
    retained_ids = getfield.(retained.case.reference.ports, :id)
    reduced_ids = getfield.(reduced.case.reference.ports, :id)
    keep = [something(findfirst(==(id), retained_ids), 0) for id in reduced_ids]
    any(iszero, keep) && throw(ArgumentError(
        "reduced ports do not form a subset of retained ports",
    ))
    reference = _port_tensor(retained, quantity)
    expected = _kron(reference, keep, quantity)
    observed = _port_tensor(reduced, quantity)
    indices = _frequency_indices(
        frequencies(retained.actual), frequencies(reduced.actual)
    )
    expected = expected[:, :, indices]
    tolerance = _tolerance(MatrixCheck{quantity}())
    metrics = _metrics(
        observed, expected, collect(frequencies(reduced.actual)), tolerance
    )
    check = PhysicalCheck{Symbol("kron_$(quantity)")}()
    return Comparison(
        check,
        retained.case.fidelity isa Exact && reduced.case.fidelity isa Exact ?
        _verdict(metrics, tolerance) : Diagnostic(),
        tolerance,
        metrics,
        "compared the reduced case with a Schur reduction of cohort " *
        retained.case.provenance.cohort
    )
end

function _append_kron!(trials)
    cohorts = Dict{String, Vector{Int}}()
    for (index, trial) in pairs(trials)
        trial.actual === nothing && continue
        push!(get!(cohorts, trial.case.provenance.cohort, Int[]), index)
    end
    for indices in values(cohorts)
        retained = findfirst(index -> trials[index].case.reduction isa Retained, indices)
        reduced = findfirst(index -> trials[index].case.reduction isa Reduced, indices)
        retained === nothing && continue
        reduced === nothing && continue
        retained_trial = trials[indices[retained]]
        reduced_trial = trials[indices[reduced]]
        for quantity in (:Z, :Y)
            comparison = try
                _kron_comparison(retained_trial, reduced_trial, quantity)
            catch error
                error isa InterruptException && rethrow()
                Comparison(
                    PhysicalCheck{Symbol("kron_$(quantity)")}(),
                    Diagnostic(),
                    _tolerance(MatrixCheck{quantity}()),
                    nothing,
                    sprint(showerror, error)
                )
            end
            push!(reduced_trial.comparisons, comparison)
        end
    end
    return trials
end

"""Run every selected case and return one report."""
function gauntlet(suite::Suite)
    started = now()
    trials = Trial[]
    for (index, id) in enumerate(suite.ids)
        (index == 1 || index % 25 == 0) && @info("Running Gauntlet suite",
            suite=suite.name,
            case=index,
            total=length(suite.ids),
            id)
        case = suite.dataset[id]
        push!(trials,
            gauntlet(
                case;
                performance = id in suite.performance,
                checks = suite.checks === nothing ? case.checks : suite.checks,
                policy = suite.policy,
                tolerances = suite.tolerances
            ))
    end
    _append_kron!(trials)
    return Report(suite, trials, started, now())
end
