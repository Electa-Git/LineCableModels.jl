_calibration_bucket(id::AbstractString) = parse(Int, first(id, 2); base = 16) % 5

function _round125(value::Real)
    value <= 0 && return 0.0
    exponent = floor(Int, log10(value))
    scale = 10.0^exponent
    mantissa = value / scale
    rounded = mantissa <= 1 ? 1.0 : mantissa <= 2 ? 2.0 : mantissa <= 5 ? 5.0 : 10.0
    return rounded * scale
end

function _fit_calibration(
        dataset::Dataset,
        check::FitCheck,
        family::Type{<:Family};
        atol::Real = 1e-12
)
    records = NamedTuple[]
    tolerance = Tolerance(1e-3, Float64(atol))
    for id in keys(dataset.cases)
        case = dataset[id]
        case.family isa family || continue
        comparison = compare_fit(
            case.reference.fit, case.reference.modes, check, tolerance
        )
        comparison.metrics === nothing && continue
        push!(records, (;
            id,
            bucket = _calibration_bucket(id),
            error = comparison.metrics.max_rel
        ))
    end
    calibration = filter(record -> record.bucket != 4, records)
    validation = filter(record -> record.bucket == 4, records)
    isempty(calibration) && return nothing
    isempty(validation) && return nothing
    calibration_max = maximum(getfield.(calibration, :error))
    validation_max = maximum(getfield.(validation, :error))
    threshold = _round125(2calibration_max)
    return (;
        check = _check_name(check),
        family = lowercase(String(nameof(family))),
        calibration_count = length(calibration),
        validation_count = length(validation),
        calibration_max,
        validation_max,
        threshold,
        capped = threshold <= 1e-3,
        held_out_pass = validation_max <= threshold
    )
end

"""Derive source-internal fit thresholds without changing frozen configuration."""
function calibrate_fit(dataset::Dataset)
    result = NamedTuple[]
    for check in (FitCheck{:Yc}(), FitCheck{:H}())
        for family in (Coax, Overhead, Mixed, Pipe)
            record = _fit_calibration(dataset, check, family)
            record === nothing || push!(result, record)
        end
    end
    return result
end

function _polar_log_interpolate(values, frequency, target)
    upper = searchsortedfirst(frequency, target)
    1 < upper <= length(frequency) || throw(DomainError(
        target, "interpolation frequency lies outside the detailed source range"
    ))
    lower = upper - 1
    weight = (log(target) - log(frequency[lower])) /
             (log(frequency[upper]) - log(frequency[lower]))
    first_value = @view values[:, :, lower]
    last_value = @view values[:, :, upper]
    first_phase = angle.(first_value)
    phase_step = mod.(angle.(last_value) .- first_phase .+ π, 2π) .- π
    magnitude = exp.(
        (1 - weight) .* log.(abs.(first_value)) .+
        weight .* log.(abs.(last_value))
    )
    return magnitude .* cis.(first_phase .+ weight .* phase_step)
end
