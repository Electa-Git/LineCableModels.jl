module PSCADIO

using ..Gauntlet: DelayGroup, Fit, Modes, RationalColumn
using LinearAlgebra
using Printf
using SHA

export Block,
       Detailed,
       Ordinary,
       ParsedInput,
       block,
       combine_polar,
       combine_polar_pairs,
       read_detailed,
       read_fit,
       read_input,
       read_ordinary,
       terminal_values

"""Represent one nested PSCAD line-constants input section."""
mutable struct Block
    name::String
    value::Any
    fields::Dict{String, Any}
    children::Vector{Block}
end

function Block(name::AbstractString; value = nothing)
    Block(String(name), value, Dict{String, Any}(), Block[])
end

"""Hold a parsed PSCAD line/cable input and its source kind."""
struct ParsedInput
    kind::Symbol
    root::Block
end

"""Hold one PSCAD detailed frequency table."""
struct Detailed{T <: Real}
    frequency::Vector{T}
    values::Matrix{T}
    label::String
end

"""Convert one detailed table of interleaved magnitude/degree columns."""
function terminal_values(table::Detailed{T}) where {T}
    iseven(size(table.values, 2)) || throw(DimensionMismatch(
        "terminal output requires interleaved magnitude/phase columns",
    ))
    magnitude = @view table.values[:, 1:2:end]
    phase = @view table.values[:, 2:2:end]
    return magnitude .* cispi.(phase ./ T(180))
end

"""Combine interleaved calculated/fitted magnitude and phase columns."""
function combine_polar_pairs(
        magnitude::Detailed{T},
        phase::Detailed{T};
        matrix::Bool = true
) where {T}
    magnitude.frequency == phase.frequency || throw(DimensionMismatch(
        "magnitude and phase frequencies differ",
    ))
    size(magnitude.values) == size(phase.values) || throw(DimensionMismatch(
        "magnitude and phase table shapes differ",
    ))
    iseven(size(magnitude.values, 2)) || throw(DimensionMismatch(
        "calculated/fitted PSCAD tables require paired columns",
    ))
    calculated = Detailed(
        magnitude.frequency,
        magnitude.values[:, 1:2:end],
        magnitude.label
    )
    calculated_phase = Detailed(
        phase.frequency,
        phase.values[:, 1:2:end],
        phase.label
    )
    fitted = Detailed(
        magnitude.frequency,
        magnitude.values[:, 2:2:end],
        magnitude.label
    )
    fitted_phase = Detailed(
        phase.frequency,
        phase.values[:, 2:2:end],
        phase.label
    )
    return (
        combine_polar(calculated, calculated_phase; matrix),
        combine_polar(fitted, fitted_phase; matrix)
    )
end

"""Hold ordinary phase and sequence matrices at one frequency."""
struct Ordinary{T <: Real}
    frequency::T
    Z::Union{Nothing, Matrix{Complex{T}}}
    Y::Union{Nothing, Matrix{Complex{T}}}
    sequence_transform::Union{Nothing, Matrix{Complex{T}}}
    sequence_Z::Union{Nothing, Matrix{Complex{T}}}
    sequence_Y::Union{Nothing, Matrix{Complex{T}}}
end

function _number(text::AbstractString)
    value = replace(strip(text), 'D' => 'E', 'd' => 'e')
    normalized = replace(value, r"^([+-])\." => s"\g<1>0.")
    normalized = replace(normalized, r"^\." => "0.")
    # Some PSCAD writers omit the exponent marker when the formatted field
    # exhausts its width, for example `0.15942974-142`.  This is ordinary
    # Fortran scientific notation, not subtraction.
    normalized = replace(
        normalized,
        r"^([+-]?(?:\d+(?:\.\d*)?|\.\d+))([+-]\d+)$" => s"\g<1>E\g<2>"
    )
    parsed_int = tryparse(Int, normalized)
    parsed_int === nothing || return parsed_int
    parsed_float = tryparse(Float64, normalized)
    return parsed_float === nothing ? nothing : parsed_float
end

function _value(text::AbstractString)
    stripped = strip(text)
    isempty(stripped) && return ""
    pieces = split(stripped)
    parsed = [_number(piece) for piece in pieces]
    if all(!isnothing, parsed)
        values = something.(parsed)
        return length(values) == 1 ? only(values) : values
    end
    return stripped
end

function _significant_lines(path::AbstractString)
    return [strip(line)
            for line in eachline(path)
            if !isempty(strip(line)) && !startswith(strip(line), "!")]
end

"""Read a nested PSCAD `.cli` or `.tli` input file."""
function read_input(path::AbstractString)
    extension = lowercase(splitext(path)[2])
    extension in (".cli", ".tli") || throw(ArgumentError(
        "PSCAD line input must use .cli or .tli",
    ))
    lines = _significant_lines(path)
    root = Block("root")
    stack = Block[root]
    pending_name = nothing
    pending_value = nothing
    index = 1
    while index <= length(lines)
        line = lines[index]
        if line == "{"
            pending_name === nothing && throw(ArgumentError(
                "unnamed block at line $index in $path",
            ))
            child = Block(pending_name; value = pending_value)
            push!(last(stack).children, child)
            push!(stack, child)
            pending_name = nothing
            pending_value = nothing
        elseif line == "}"
            length(stack) > 1 || throw(ArgumentError(
                "unmatched closing block at line $index in $path",
            ))
            pop!(stack)
        elseif endswith(line, ':')
            pending_name = strip(line[1:(end - 1)])
            pending_value = nothing
        elseif occursin('=', line)
            key, raw = split(line, '='; limit = 2)
            key = strip(key)
            value = _value(raw)
            if index < length(lines) && lines[index + 1] == "{"
                pending_name = key
                pending_value = value
            else
                haskey(last(stack).fields, key) && throw(ArgumentError(
                    "duplicate field '$key' in PSCAD block $(last(stack).name)",
                ))
                last(stack).fields[key] = value
            end
        else
            throw(ArgumentError("unrecognized PSCAD input line '$line' in $path"))
        end
        index += 1
    end
    length(stack) == 1 || throw(ArgumentError("unterminated PSCAD block in $path"))
    has_tower = any(child -> child.name == "Line Constants Tower", root.children)
    has_coax = any(child -> child.name == "Coax Cable", root.children)
    kind = any(child -> child.name == "Pipe-Type Cable", root.children) ? :pipe :
           has_tower && has_coax ? :mixed :
           extension == ".tli" || has_tower ? :overhead : :coax
    return ParsedInput(kind, root)
end

"""Return all immediate children of `parent` with the requested name."""
function block(parent::Block, name::AbstractString)
    [child for child in parent.children if child.name == name]
end

function block(parent::Block, name::AbstractString, index::Integer)
    matches = block(parent, name)
    checkbounds(matches, index)
    return matches[index]
end

function _float_row(line::AbstractString)
    values = Float64[]
    for field in split(strip(line))
        value = _number(field)
        value isa Real || throw(ArgumentError("invalid PSCAD numeric field '$field'"))
        push!(values, Float64(value))
    end
    return values
end

"""Read one detailed PSCAD frequency table."""
function read_detailed(path::AbstractString; allow_empty::Bool = false)
    lines = collect(eachline(path))
    header = findfirst(line -> occursin("LOG10(FN)", line), lines)
    header === nothing && throw(ArgumentError("missing detailed-output header in $path"))
    rows = Vector{Vector{Float64}}()
    for line in @view(lines[(header + 1):end])
        isempty(strip(line)) && continue
        push!(rows, _float_row(line))
    end
    if isempty(rows)
        allow_empty && return nothing
        throw(ArgumentError("detailed output contains no data: $path"))
    end
    width = length(first(rows))
    all(row -> length(row) == width, rows) || throw(DimensionMismatch(
        "detailed PSCAD output has inconsistent row widths: $path",
    ))
    width >= 3 || throw(DimensionMismatch("detailed PSCAD output needs data columns"))
    table = reduce(vcat, permutedims.(rows))
    all(isfinite, table) ||
        throw(DomainError(path, "detailed output contains nonfinite data"))
    return Detailed(table[:, 2], table[:, 3:end], String(strip(lines[header])))
end

"""Combine PSCAD magnitude and degree-phase detailed tables."""
function combine_polar(magnitude::Detailed{T}, phase::Detailed{T}; matrix = true) where {T}
    magnitude.frequency == phase.frequency || throw(DimensionMismatch(
        "magnitude and phase frequencies differ",
    ))
    size(magnitude.values) == size(phase.values) || throw(DimensionMismatch(
        "magnitude and phase table shapes differ",
    ))
    values = magnitude.values .* cispi.(phase.values ./ T(180))
    matrix || return permutedims(values)
    n = isqrt(size(values, 2))
    n * n == size(values, 2) || throw(DimensionMismatch(
        "detailed table does not contain a square matrix",
    ))
    result = Array{Complex{T}, 3}(undef, n, n, size(values, 1))
    for frequency in axes(values, 1), row in 1:n, column in 1:n
        result[row, column, frequency] = values[frequency, (row - 1) * n + column]
    end
    return result
end

function _complex(text::AbstractString)
    pieces = split(strip(text), ',')
    length(pieces) == 2 || throw(ArgumentError("invalid PSCAD complex value '$text'"))
    real_part = _number(first(pieces))
    imag_part = _number(last(pieces))
    real_part isa Real && imag_part isa Real || throw(ArgumentError(
        "invalid PSCAD complex value '$text'",
    ))
    return complex(Float64(real_part), Float64(imag_part))
end

function _matrix_after(lines, marker)
    start = findfirst(line -> occursin(marker, line), lines)
    start === nothing && return nothing
    rows = Vector{Vector{ComplexF64}}()
    for line in @view(lines[(start + 1):end])
        stripped = strip(line)
        isempty(stripped) && (!isempty(rows) && break; continue)
        occursin(',', stripped) || (!isempty(rows) && break; continue)
        push!(rows, _complex.(split(stripped)))
    end
    isempty(rows) && return nothing
    width = length(first(rows))
    all(row -> length(row) == width, rows) || throw(DimensionMismatch(
        "ordinary PSCAD matrix '$marker' has inconsistent rows",
    ))
    length(rows) == width || throw(DimensionMismatch(
        "ordinary PSCAD matrix '$marker' is not square",
    ))
    return reduce(vcat, permutedims.(rows))
end

"""Read ordinary PSCAD phase and sequence matrices."""
function read_ordinary(path::AbstractString)
    lines = collect(eachline(path))
    frequency_line = findfirst(line -> occursin("PHASE DOMAIN DATA @", line), lines)
    if frequency_line === nothing
        frequency_line = findfirst(line -> occursin("SEQUENCE COMPONENT DATA @", line), lines)
    end
    frequency_line === nothing && throw(ArgumentError(
        "ordinary PSCAD output contains no phase or sequence frequency: $path",
    ))
    capture = match(r"@\s*([+\-0-9.Ee]+)\s*Hz", lines[frequency_line])
    capture === nothing && throw(ArgumentError("cannot parse ordinary output frequency"))
    frequency = parse(Float64, capture.captures[1])
    return Ordinary(
        frequency,
        _matrix_after(lines, "SERIES IMPEDANCE MATRIX (Z)"),
        _matrix_after(lines, "SHUNT ADMITTANCE MATRIX (Y)"),
        _matrix_after(lines, "SEQUENCE TRANSFORM MATRIX"),
        _matrix_after(lines, "SEQUENCE IMPEDANCE MATRIX (Zsq)"),
        _matrix_after(lines, "SEQUENCE ADMITTANCE MATRIX (Ysq)")
    )
end

mutable struct NumberStream
    lines::Vector{String}
    index::Int
end

function _next_number!(stream::NumberStream)
    while stream.index <= length(stream.lines)
        line = strip(stream.lines[stream.index])
        stream.index += 1
        isempty(line) && continue
        startswith(line, '!') && continue
        value = _number(first(split(line)))
        value === nothing && continue
        return Float64(value)
    end
    throw(EOFError())
end

_next_complex!(stream::NumberStream) = complex(_next_number!(stream), _next_number!(stream))

function _seek_after!(stream::NumberStream, phrase::AbstractString)
    found = findnext(line -> occursin(phrase, line), stream.lines, stream.index)
    found === nothing && throw(ArgumentError("missing '$phrase' in PSCAD fit"))
    stream.index = found + 1
    return stream
end

"""Read PSCAD `.clo` or `.tlo` rational fit coefficients."""
function read_fit(path::AbstractString; frequency_range = (1e-3, 1e7))
    extension = lowercase(splitext(path)[2])
    extension in (".clo", ".tlo") || throw(ArgumentError(
        "PSCAD fitted model must use .clo or .tlo",
    ))
    lines = collect(eachline(path))
    stream = NumberStream(lines, 1)
    _seek_after!(stream, "N.o. conductors")
    n = Int(_next_number!(stream))
    n > 0 || throw(DomainError(n, "fit conductor count must be positive"))

    columns = RationalColumn{Float64}[]
    for _ in 1:n
        order = Int(_next_number!(stream))
        poles = Vector{ComplexF64}(undef, order)
        residues = Matrix{ComplexF64}(undef, n, order)
        for pole in 1:order
            poles[pole] = _next_complex!(stream)
            for row in 1:n
                residues[row, pole] = _next_complex!(stream)
            end
        end
        push!(columns, RationalColumn(poles, residues))
    end
    constant = Matrix{Float64}(undef, n, n)
    for column in 1:n, row in 1:n

        constant[row, column] = _next_number!(stream)
    end

    _seek_after!(stream, "Fitting parameters for the propagation function")
    groups_count = Int(_next_number!(stream))
    delays = [_next_number!(stream) for _ in 1:groups_count]
    orders = Vector{Int}(undef, groups_count)
    poles = Vector{Vector{ComplexF64}}(undef, groups_count)
    for group in 1:groups_count
        orders[group] = Int(_next_number!(stream))
        poles[group] = [_next_complex!(stream) for _ in 1:orders[group]]
    end
    groups = DelayGroup{Float64}[]
    for group in 1:groups_count
        residues = Array{ComplexF64, 3}(undef, n, n, orders[group])
        for column in 1:n, row in 1:n, pole in 1:orders[group]
            residues[row, column, pole] = _next_complex!(stream)
        end
        push!(groups, DelayGroup(delays[group], poles[group], residues))
    end
    return Fit(columns, constant, groups, frequency_range)
end

end
