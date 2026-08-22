struct DetailedTable
    frequency::Vector{Float64}
    values::Matrix{Float64}
end

const PSCAD_FREQUENCY_RTOL = 5.0e-8

function _pscad_number(text::AbstractString)
    value = replace(strip(text), 'D' => 'E', 'd' => 'e')
    value = replace(value, r"^([+-])\." => s"\g<1>0.")
    value = replace(value, r"^\." => "0.")
    value = replace(
        value,
        r"^([+-]?(?:\d+(?:\.\d*)?|\.\d+))([+-]\d+)$" => s"\g<1>E\g<2>"
    )
    parsed = tryparse(Float64, value)
    parsed === nothing && throw(ArgumentError("invalid PSCAD numeric field $(repr(text))"))
    return parsed
end

function _float_row(line::AbstractString)
    return _pscad_number.(split(strip(line)))
end

function _read_detailed(path::AbstractString)
    lines = collect(eachline(path))
    header = findfirst(line -> occursin("LOG10(FN)", line), lines)
    header === nothing && throw(ArgumentError(
        "missing PSCAD detailed-output header in $path",
    ))
    rows = [_float_row(line)
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
            rtol = PSCAD_FREQUENCY_RTOL,
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
        basis = :per_length
    )
end
