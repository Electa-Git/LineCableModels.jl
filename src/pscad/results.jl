struct DetailedTable
    frequency::Vector{Float64}
    values::Matrix{Float64}
end

const PSCAD_FREQUENCY_RTOL = 5.0e-8

function _pscad_number(text::AbstractString; required::Bool = true)
    value = replace(strip(text), 'D' => 'E', 'd' => 'e')
    value = replace(value, r"^([+-])\." => s"\g<1>0.")
    value = replace(value, r"^\." => "0.")
    value = replace(
        value,
        r"^([+-]?(?:\d+(?:\.\d*)?|\.\d+))([+-]\d+)$" => s"\g<1>E\g<2>"
    )
    parsed = tryparse(Float64, value)
    parsed === nothing && required && throw(ArgumentError("invalid PSCAD numeric field $(repr(text))"))
    return parsed
end

function _read_detailed(path::AbstractString)
    lines = collect(eachline(path))
    header = findfirst(line -> occursin("LOG10(FN)", line), lines)
    header === nothing && throw(ArgumentError(
        "missing PSCAD detailed-output header in $path",
    ))
    rows = [_pscad_number.(split(strip(line)))
            for line in @view(lines[(header + 1):end]) if !isempty(strip(line))]
    isempty(rows) && throw(ArgumentError("PSCAD detailed output contains no data: $path"))
    width = length(first(rows))
    width >= 3 || throw(DimensionMismatch("PSCAD detailed output needs data columns"))
    all(row -> length(row) == width, rows) || throw(DimensionMismatch(
        "PSCAD detailed output has inconsistent row widths: $path",
    ))
    table = reduce(vcat, permutedims.(rows))
    all(isfinite, table) || throw(DomainError(
        path,
        "PSCAD detailed output contains nonfinite data"
    ))
    return DetailedTable(table[:, 2], table[:, 3:end])
end

function _combine_polar(magnitude::DetailedTable, phase::DetailedTable)
    magnitude.frequency == phase.frequency || throw(ArgumentError(
        "PSCAD magnitude and phase frequencies differ",
    ))
    size(magnitude.values) == size(phase.values) || throw(DimensionMismatch(
        "PSCAD magnitude and phase table shapes differ",
    ))
    flat = magnitude.values .* cispi.(phase.values ./ 180)
    dimension = isqrt(size(flat, 2))
    dimension^2 == size(flat, 2) || throw(DimensionMismatch(
        "PSCAD detailed table does not contain a square matrix",
    ))
    result = Array{ComplexF64, 3}(undef, dimension, dimension, size(flat, 1))
    for frequency in axes(flat, 1), row in 1:dimension, column in 1:dimension
        result[row, column, frequency] = flat[frequency, (row - 1) * dimension + column]
    end
    return result
end

function _result_path(directory::AbstractString, name::AbstractString)
    path = joinpath(directory, name)
    isfile(path) || throw(ArgumentError("required PSCAD result is missing: $path"))
    return path
end

function _case_frequencies(observed::AbstractVector, expected::AbstractVector)
    length(observed) == length(expected) || throw(DimensionMismatch(
        "PSCAD emitted $(length(observed)) frequencies; expected $(length(expected))",
    ))
    for index in eachindex(observed, expected)
        isapprox(
            observed[index], expected[index];
            rtol = max(PSCAD_FREQUENCY_RTOL, 8Float64(eps(float(one(eltype(expected)))))),
            atol = 0.0
        ) || throw(ArgumentError(
            "PSCAD frequency at row $index is $(observed[index]); " *
            "expected $(expected[index]) within the PSCAD text-output resolution",
        ))
    end
    return Float64.(expected)
end

function read_pscad_result(
        output_dir::AbstractString,
        expected_frequencies::AbstractVector,
        expected_size::NTuple{3, Int}
)
    isdir(output_dir) ||
        throw(ArgumentError("PSCAD output directory is missing: $output_dir"))
    zm = _read_detailed(_result_path(output_dir, "result_zm.out"))
    zp = _read_detailed(_result_path(output_dir, "result_zp.out"))
    ym = _read_detailed(_result_path(output_dir, "result_ym.out"))
    yp = _read_detailed(_result_path(output_dir, "result_yp.out"))
    observed_frequencies = zm.frequency
    zp.frequency == observed_frequencies && ym.frequency == observed_frequencies &&
    yp.frequency == observed_frequencies || throw(ArgumentError(
        "PSCAD detailed outputs use different frequency samples",
    ))
    frequencies_value = _case_frequencies(observed_frequencies, expected_frequencies)
    impedance = _combine_polar(zm, zp)
    admittance = _combine_polar(ym, yp)
    size(impedance) == expected_size || throw(DimensionMismatch(
        "PSCAD impedance has size $(size(impedance)); expected $expected_size",
    ))
    size(admittance) == expected_size || throw(DimensionMismatch(
        "PSCAD admittance has size $(size(admittance)); expected $expected_size",
    ))
    return LineParameters(
        PhaseDomain,
        impedance,
        admittance,
        frequencies_value;
        basis = :pul,
        details = ComputationDetails(; native_frequencies = copy(observed_frequencies))
    )
end

# Native headings delimit both complete reports and logs interrupted by Tline.
function _pscad_fit_section(text, heading, ending)
    return match(Regex(raw"(?ms)^[ \t]*" * heading * raw"[^\r\n]*\n(.*?)(?=^[ \t]*(?:" * ending * raw")|\z)"), text)
end

"""
$(TYPEDSIGNATURES)

Read the curve-fit and passivity reports from a native Tline log.

# Arguments

- `path`: Local line-constants log file.

# Returns

A named tuple with these fields. Absent sections and values are `nothing`.

- `admittance`: `before_dc` and `after_dc` error records, and the reported `poles`
  count. Each error record has `maximum_error_percent` and `rms_error_percent`.
- `propagation`: The same error records, `poles` per delay group, and
  `poles_source` (`:final_fit` or `:last_hmode`).
- `maximum_residue_pole_ratio`: Largest ratio of residue to pole in the final fit.
- `scan`: `band` with `lower_Hz` and `upper_Hz`, violation `ranges` with those
  limits and `smallest_eigenvalue`, the overall `smallest_eigenvalue`, native
  `outcome` and `report` text. A clean scan has empty ranges.
- `enforcement`: Native `outcome`, `initial_smallest_eigenvalue`,
  `final_smallest_eigenvalue`, overall `smallest_eigenvalue`,
  `maximum_abs_error_percent`, `ranges`, `report` and `iterations`.
  Each inner iteration has `main` and `inner` indices, `frequencies` \\[Hz\\] and
  `eigenvalues`. An iteration without a table has empty sample vectors.
- `reports`: `fitting`, `before_dc` and `after_dc` text, with LF line endings.

The `before_dc` errors describe the last fit before any DC correction, including
runs with DC correction disabled. Enforcement minima describe the first and last
reported eigenvalue tables. The outcome states whether enforcement succeeded.
Errors use percent. Frequencies use Hz. Eigenvalues and ratios of residue to pole
keep the native normalization.
"""
function read_pscad_fit_report(path::AbstractString)
    text = replace(read(path, String), "\r\n" => "\n")
    number = raw"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
    fitting = _pscad_fit_section(text, "CURVE FITTING OUTPUT", raw"LCP solve constants complete\.")
    fit_text = fitting === nothing ? "" : fitting.captures[1]
    correction = match(r"(?m)^[ \t]*Applying DC Correction[ \t]*$", fit_text)
    before, after = correction === nothing ? (fit_text, "") :
        Tuple(split(fit_text, correction.match; limit = 2))
    reports = (fitting = fitting === nothing ? nothing : String(strip(fit_text)),
        before_dc = correction === nothing ? nothing : String(strip(before)),
        after_dc = correction === nothing ? nothing : String(strip(after)))

    yc_rows = map((before, after)) do section
        yc = _pscad_fit_section(section, "Fitting Characteristic Admittance Yc:", "Fitting Propagation Function H:")
        yc === nothing && return nothing
        return match(Regex(raw"(?m)^[ \t]*Yc:[^\n]*\n[ \t]*Fitting Error[ \t]+of Poles[ \t]*\n[ \t-]+\n[ \t]*(" *
            number * raw")[ \t]*%[ \t]+(\d+)[ \t]*$"), yc.captures[1])
    end
    yc_errors = map(yc_rows) do row
        row === nothing ? nothing :
            (maximum_error_percent = _pscad_number(row[1]), rms_error_percent = nothing)
    end
    yc_final = yc_rows[2] === nothing ? yc_rows[1] : yc_rows[2]
    admittance = yc_final === nothing ? nothing :
        (before_dc = yc_errors[1], after_dc = yc_errors[2], poles = parse(Int, yc_final[2]))

    h_before, h_after, h_poles, h_source, ratio = nothing, nothing, nothing, nothing, nothing
    propagation_section = _pscad_fit_section(before, "Fitting Propagation Function H:", "Applying DC Correction")
    if propagation_section !== nothing
        attempts = collect(eachmatch(Regex(raw"(?m)^[ \t]*Hphase:[^\n]*\n[ \t]*Fitting Error[ \t]+RMS Error[ \t]+Residue/Pole Ratio[ \t]*\n[ \t-]+\n[ \t]*(" *
            number * raw")[ \t]*%[ \t]+(" * number * raw")[ \t]*%[ \t]+(" * number * raw")"), propagation_section[1]))
        if !isempty(attempts)
            row = last(attempts)
            h_before = (maximum_error_percent = _pscad_number(row[1]), rms_error_percent = _pscad_number(row[2]))
            ratio = _pscad_number(row[3])
        end
        modal = collect(eachmatch(r"(?ms)^[ \t]*Hmode:[^\n]*\n(.*?)(?=^[ \t]*Hphase:|\z)", propagation_section[1]))
        if !isempty(modal)
            counts = [parse(Int, row[1]) for row in eachmatch(Regex(raw"(?m)^[ \t]*\d+[ \t]+" * number *
                raw"[ \t]*%[ \t]+(\d+)[ \t]+" * number * raw"[ \t]+ms[ \t]*$"), last(modal)[1])]
            if !isempty(counts)
                h_poles, h_source = counts, :last_hmode
            end
        end
    end
    final = _pscad_fit_section(fit_text, "Final Fit:", "Maximum Residue/Pole Ratio:")
    if final !== nothing
        row = match(Regex(raw"(?m)^[ \t]*(" * number * raw")[ \t]*%[ \t]+(" * number * raw")[ \t]*%[ \t]+(\d+)[ \t]+" *
            number * raw"[ \t]+ms[ \t]*$"), final[1])
        if row !== nothing
            h_poles = [parse(Int, row[3]); [parse(Int, item[1]) for item in eachmatch(
                Regex(raw"(?m)^[ \t]*(\d+)[ \t]+" * number * raw"[ \t]+ms[ \t]*$"), final[1])]]
            h_source = :final_fit
            if correction !== nothing && final.offset > correction.offset
                h_after = (maximum_error_percent = _pscad_number(row[1]), rms_error_percent = _pscad_number(row[2]))
            end
        end
    end
    final_ratio = match(Regex(raw"(?m)^[ \t]*Maximum Residue/Pole Ratio:[ \t]*(" * number * raw")[ \t]*$"), fit_text)
    final_ratio === nothing || (ratio = _pscad_number(final_ratio[1]))
    propagation = h_before === nothing && h_after === nothing && h_poles === nothing ? nothing :
        (before_dc = h_before, after_dc = h_after, poles = h_poles, poles_source = h_source)

    enforcement_heading = raw"Passivity Enforcement[ \t]*$"
    scan_section = _pscad_fit_section(text, "Checking for Passivity Violations Between", enforcement_heading * raw"|LCP solve constants complete\.")
    enforcement_section = _pscad_fit_section(text, enforcement_heading, raw"LCP solve constants complete\.")
    range_pattern = Regex(raw"(?m)^[ \t]*From[ \t]+(" * number * raw")[ \t]+Hz[ \t]+to[ \t]+(" *
        number * raw")[ \t]+Hz[ \t]+(" * number * raw")[ \t]*$")
    passivity_ranges = map((scan_section, enforcement_section)) do section
        section === nothing && return nothing
        return [(lower_Hz = _pscad_number(row[1]), upper_Hz = _pscad_number(row[2]),
            smallest_eigenvalue = _pscad_number(row[3])) for row in eachmatch(range_pattern, section[1])]
    end
    scan = nothing
    if scan_section !== nothing
        band = match(Regex(raw"Checking for Passivity Violations Between[ \t]+(" * number * raw")[ \t]+Hz[ \t]+and[ \t]+(" *
            number * raw")[ \t]+Hz"), scan_section.match)
        outcomes = collect(eachmatch(r"(?m)^[ \t]*Passivity violations detected\. See help for further details[ \t]*$", scan_section[1]))
        ranges = passivity_ranges[1]
        scan = (band = band === nothing ? nothing : (lower_Hz = _pscad_number(band[1]), upper_Hz = _pscad_number(band[2])),
            ranges, smallest_eigenvalue = isempty(ranges) ? nothing : minimum(row.smallest_eigenvalue for row in ranges),
            outcome = isempty(outcomes) ? nothing : String(strip(last(outcomes).match)), report = String(strip(scan_section.match)))
    end
    enforcement = nothing
    if enforcement_section !== nothing
        body = enforcement_section[1]
        # cspell:ignore sucessfully
        outcomes = collect(eachmatch(r"(?m)^[ \t]*(?:Constrained solution failed|There are no passivity violations\. Passivity is sucessfully enforced)[ \t]*$", body))
        table_pattern = Regex(raw"(?m)^[ \t]*Frequency \[Hz\][ \t]+Negative eigenvalue minimas[ \t]*\n[ \t.]+\n((?:[ \t]*" *
            number * raw"[ \t]+" * number * raw"[ \t]*\n)+)")
        point_pattern = Regex(raw"(?m)^[ \t]*(" * number * raw")[ \t]+(" * number * raw")[ \t]*$")
        tables = collect(eachmatch(table_pattern, body))
        minima = [minimum(_pscad_number(row[2]) for row in eachmatch(point_pattern, table[1])) for table in tables]
        markers = collect(eachmatch(r"(?m)^[ \t]*(Main|Inner) Iteration[ \t]+#[ \t]*(\d+)[ \t]*$", body))
        iterations = NamedTuple{(:main, :inner, :frequencies, :eigenvalues), Tuple{Union{Nothing, Int}, Int, Vector{Float64}, Vector{Float64}}}[]
        main = nothing
        for (index, marker) in enumerate(markers)
            if marker[1] == "Main"
                main = parse(Int, marker[2])
                continue
            end
            stop = index == length(markers) ? lastindex(body) : prevind(body, markers[index + 1].offset)
            table = match(table_pattern, SubString(body, marker.offset, stop))
            points = table === nothing ? Tuple{Float64, Float64}[] :
                [(_pscad_number(row[1]), _pscad_number(row[2])) for row in eachmatch(point_pattern, table[1])]
            push!(iterations, (main, inner = parse(Int, marker[2]), frequencies = Float64[first(row) for row in points],
                eigenvalues = Float64[last(row) for row in points]))
        end
        error_match = match(Regex(raw"(?m)^[ \t]*Maximum abs\. error[ \t]+(" * number * raw")[ \t]*%[ \t]*$"), body)
        enforcement = (outcome = isempty(outcomes) ? nothing : String(strip(last(outcomes).match)),
            initial_smallest_eigenvalue = isempty(minima) ? nothing : first(minima),
            final_smallest_eigenvalue = isempty(minima) ? nothing : last(minima),
            smallest_eigenvalue = isempty(minima) ? (isempty(passivity_ranges[2]) ? nothing :
                minimum(row.smallest_eigenvalue for row in passivity_ranges[2])) : minimum(minima),
            maximum_abs_error_percent = error_match === nothing ? nothing : _pscad_number(error_match[1]),
            ranges = passivity_ranges[2], iterations = isempty(markers) ? nothing : iterations,
            report = String(strip(enforcement_section.match)))
    end
    return (; admittance, propagation, maximum_residue_pole_ratio = ratio, scan, enforcement, reports)
end

"""
$(TYPEDEF)

Report the earliest nonfinite channel sample in PSCAD transient output.
The first channel in INF order identifies simultaneous nonfinite samples.

$(TYPEDFIELDS)
"""
struct NonfiniteOutputError <: Exception
    "Local path of the numbered OUT file containing the sample."
    output_path::String
    "PGB channel name from the INF file."
    channel::String
    "Time of the first nonfinite sample \\[s\\]."
    time::Float64
end

function Base.showerror(io::IO, error::NonfiniteOutputError)
    print(io, "PSCAD output has nonfinite samples in channel ", repr(error.channel),
        " at t=", error.time, " s: ", error.output_path)
end

"""
    read_pscad_transient(directory::AbstractString; stem::AbstractString)

Read PSCAD INF channel names and numbered OUT files with time followed by at
most ten channels. Return `(time, channels)`, with time in seconds and channels
keyed by their PGB names in the recorded units. Accept Fortran exponents.
Throw `NonfiniteOutputError` for the earliest nonfinite channel sample across
all numbered files. Invalid time coordinates and missing output remain errors.
"""
function read_pscad_transient(directory::AbstractString; stem::AbstractString)
    occursin(r"^[A-Za-z0-9][A-Za-z0-9_]{0,19}$", stem) ||
        throw(ArgumentError("invalid PSCAD output stem $stem"))
    names = [entry.desc for entry in _inf_entries(_result_path(directory, stem * ".inf"))]
    allunique(names) || throw(ArgumentError("PSCAD channel names must be unique"))
    time = Float64[]
    channels = Dict{String, Vector{Float64}}()
    divergence = nothing
    for file_index in 1:cld(length(names), 10)
        offset = 10(file_index - 1)
        count = min(10, length(names) - offset)
        path = _result_path(directory, stem * "_" * lpad(string(file_index), 2, '0') * ".out")
        rows = Vector{Float64}[]
        for line in eachline(path)
            fields = split(strip(line))
            isempty(fields) && continue
            _pscad_number(first(fields); required = false) === nothing && continue
            push!(rows, _pscad_number.(fields))
        end
        isempty(rows) && throw(ArgumentError("PSCAD output contains no samples: $path"))
        all(row -> length(row) == count + 1, rows) ||
            throw(ArgumentError("PSCAD output has an incorrect column count: $path"))
        table = reduce(vcat, permutedims.(rows))
        times = table[:, 1]
        all(isfinite, times) || throw(ArgumentError("PSCAD output has nonfinite times: $path"))
        all(>(0), diff(times)) ||
            throw(ArgumentError("PSCAD output time must be strictly increasing: $path"))
        if file_index == 1
            time = times
        else
            times == time || throw(ArgumentError("PSCAD output time columns differ: $path"))
        end
        invalid = findfirst(!isfinite, transpose(@view table[:, 2:end]))
        if invalid !== nothing
            channel, sample = Tuple(invalid)
            if divergence === nothing || times[sample] < divergence.time
                divergence = NonfiniteOutputError(String(path), names[offset + channel], times[sample])
            end
        end
        for index in 1:count
            channels[names[offset + index]] = table[:, index + 1]
        end
    end
    divergence === nothing || throw(divergence)
    return (; time, channels)
end

function _pscad_response_channels!(voltage, current, record, coordinates,
        assignments::Tuple, end_index)
    channels = record.channels
    end_name = end_index == 1 ? "sending" : "receiving"
    assigned = Set(assignment.terminal for assignment in assignments)
    for (terminal, label) in enumerate(coordinates)
        vname, iname = "$end_name:$label:v", "$end_name:$label:i"
        haskey(channels, vname) || throw(ArgumentError("missing PSCAD channel $vname"))
        voltage[end_index, terminal, :] = channels[vname]
        if label in assigned
            haskey(channels, iname) || throw(ArgumentError("missing PSCAD channel $iname"))
            current[end_index, terminal, :] = channels[iname]
        end
    end
    return nothing
end

function _pscad_response_channels!(voltage, current, record, coordinates,
        ::Matched, end_index)
    channels = record.channels
    end_name = end_index == 1 ? "sending" : "receiving"
    quantities = ((voltage, "v"), (current, "i"))
    for (terminal, label) in enumerate(coordinates), (array, kind) in quantities
        name = "$end_name:$label:$kind"
        haskey(channels, name) || throw(ArgumentError("missing PSCAD channel $name"))
        array[end_index, terminal, :] = channels[name]
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Build a `TransientResponse` from the reader's `(time, channels)` record in the
problem's terminal order. Arrays use end × terminal × time, sending first.
Require a voltage channel for every terminal at both ends and a current channel
for every assigned terminal and matched-end branch. Unassigned open terminals
have zero current because their circuits are open.
Missing channels throw `ArgumentError` naming the channel.
The exported PGB channels already use volts \\[V\\] and amperes \\[A\\].
"""
function read_pscad_transient(record::NamedTuple{(:time, :channels)},
        problem::TransientAnalysisProblem{<:LineParametersProblem})
    coordinates = terminal_labels(problem.line.system)
    voltage = zeros(Float64, 2, length(coordinates), length(record.time))
    current = zeros(Float64, size(voltage))
    _pscad_response_channels!(voltage, current, record, coordinates, problem.sending, 1)
    _pscad_response_channels!(voltage, current, record, coordinates, problem.receiving, 2)
    return Engine.TransientAnalysis.TransientResponse(copy(record.time), voltage, current,
        coordinates, ComputationDetails(; backend = :pscad))
end
