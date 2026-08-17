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
        corpus::Corpus,
        check::FitCheck,
        family::Type{<:Family};
        atol::Real = 1e-12
)
    records = NamedTuple[]
    tolerance = Tolerance(1e-3, Float64(atol))
    for id in keys(corpus.cases)
        case = corpus[id]
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
function calibrate_fit(corpus::Corpus)
    result = NamedTuple[]
    for check in (FitCheck{:Yc}(), FitCheck{:H}())
        for family in (Coax, Overhead, Mixed, Pipe)
            record = _fit_calibration(corpus, check, family)
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

"""Cross-check ordinary 60 Hz matrices against interpolated detailed evidence."""
function ordinary_consistency(raw::AbstractString)
    root = joinpath(abspath(raw), "cases")
    isdir(root) || throw(ArgumentError("raw corpus has no cases directory: $raw"))
    impedance = Float64[]
    admittance = Float64[]
    unavailable = 0
    for id in readdir(root)
        directory = joinpath(root, id)
        names = readdir(directory)
        any(endswith(lowercase(name), "_zm.out") for name in names) || continue
        files = Dict(name => joinpath(directory, name) for name in names)
        phase, _ = _phase_reference(files)
        ordinary = _ordinary(files)
        f = frequencies(phase)
        if ordinary === nothing || ordinary.Z === nothing ||
           !(first(f) < ordinary.frequency < last(f))
            unavailable += 1
            continue
        end
        z = _polar_log_interpolate(Array(Z(phase)), f, ordinary.frequency)
        y = _polar_log_interpolate(Array(Y(phase)), f, ordinary.frequency)
        push!(impedance, norm(z - ordinary.Z) / max(norm(ordinary.Z), eps()))
        push!(admittance, norm(y - ordinary.Y) / max(norm(ordinary.Y), eps()))
    end
    return (;
        method = "polar interpolation on log frequency",
        checked = length(impedance),
        unavailable,
        impedance_max = maximum(impedance; init = 0.0),
        impedance_median = isempty(impedance) ? 0.0 : median(impedance),
        admittance_max = maximum(admittance; init = 0.0),
        admittance_median = isempty(admittance) ? 0.0 : median(admittance)
    )
end
