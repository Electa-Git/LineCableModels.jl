const REFERENCE_MIN_FREQUENCY = 0.1

function _loggrid(first_frequency, last_frequency, count)
    10.0 .^ range(
        log10(Float64(first_frequency));
        stop = log10(Float64(last_frequency)),
        length = count
    )
end

function reference_grid(frequencies::AbstractVector)
    values = Float64.(frequencies)
    length(values) >= 101 || throw(ArgumentError(
        "a cross-backend reference grid requires at least 101 frequencies",
    ))
    issorted(values) || throw(ArgumentError(
        "cross-backend frequencies must be sorted",
    ))
    all(>(0), values) || throw(DomainError(
        values,
        "cross-backend frequencies must be positive"
    ))
    length(unique(values)) == length(values) || throw(ArgumentError(
        "cross-backend frequencies must be unique",
    ))

    first_frequency = max(first(values), REFERENCE_MIN_FREQUENCY)
    last(values) > first_frequency || throw(ArgumentError(
        "the cross-backend frequency range must extend above $first_frequency Hz",
    ))
    return _loggrid(first_frequency, last(values), 101)
end

function reference_case(case_id::Symbol)
    nominal = load_case(case_id)
    selected = reference_grid(nominal.nominal_problem.frequencies)
    selected == nominal.nominal_problem.frequencies && return nominal
    return load_case(
        case_id;
        variation = ExactOverrides(frequencies = selected)
    )
end
