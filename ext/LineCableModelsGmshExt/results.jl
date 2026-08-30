const FEM_RAW_HEADER = [
    "frequency_index",
    "frequency_hz",
    "response_terminal",
    "basis_terminal",
    "real",
    "imaginary"
]

const FEM_COMPLETE_HEADER = [
    "frequency_count",
    "terminal_count",
    "expected_rows",
    "z_rows",
    "p_rows",
    "success"
]

function _nonempty_lines(path::String, run::FEMRun)
    isfile(path) || _fem_error(
        :results,
        basename(path),
        :file,
        "required FEM output is missing";
        run_directory = run.path
    )
    return filter(!isempty, strip.(readlines(path)))
end

function _parse_integer(token::AbstractString, path::String, line::Int, run::FEMRun)
    value = tryparse(Int, token)
    value === nothing && _fem_error(
        :results,
        basename(path),
        :row,
        "invalid integer at line $line: $(repr(token))";
        run_directory = run.path
    )
    return value
end

function _parse_real(
        ::Type{T}, token::AbstractString, path::String, line::Int, run::FEMRun
) where {T <: Real}
    value = tryparse(T, token)
    value === nothing && _fem_error(
        :results,
        basename(path),
        :row,
        "invalid floating-point value at line $line: $(repr(token))";
        run_directory = run.path
    )
    isfinite(value) || _fem_error(
        :results,
        basename(path),
        :row,
        "non-finite floating-point value at line $line";
        run_directory = run.path
    )
    return value
end

function _parse_raw_matrix(
        ::Type{T},
        path::String,
        frequencies::Vector{T},
        terminal_count::Int,
        run::FEMRun
) where {T <: Real}
    lines = _nonempty_lines(path, run)
    isempty(lines) && _fem_error(
        :results, basename(path), :header, "raw output is empty";
        run_directory = run.path
    )
    split(first(lines), '\t'; keepempty = true) == FEM_RAW_HEADER || _fem_error(
        :results,
        basename(path),
        :header,
        "malformed raw header: $(first(lines))";
        run_directory = run.path
    )
    expected = length(frequencies) * terminal_count * terminal_count
    length(lines) - 1 == expected || _fem_error(
        :results,
        basename(path),
        :cardinality,
        "expected $expected data rows, found $(length(lines) - 1)";
        run_directory = run.path
    )
    matrix = Array{Complex{T}, 3}(
        undef, terminal_count, terminal_count, length(frequencies)
    )
    seen = Set{NTuple{3, Int}}()
    for (line_index, line) in zip(2:length(lines), @view(lines[2:end]))
        columns = split(line, '\t'; keepempty = true)
        length(columns) == 6 || _fem_error(
            :results,
            basename(path),
            :row,
            "line $line_index has $(length(columns)) columns instead of 6";
            run_directory = run.path
        )
        frequency_index = _parse_integer(columns[1], path, line_index, run)
        response = _parse_integer(columns[3], path, line_index, run)
        basis = _parse_integer(columns[4], path, line_index, run)
        frequency_index in eachindex(frequencies) || _fem_error(
            :results, basename(path), :frequency_index,
            "frequency index $frequency_index is out of range";
            run_directory = run.path
        )
        response in 1:terminal_count || _fem_error(
            :results, basename(path), :response_terminal,
            "response terminal $response is out of range";
            run_directory = run.path
        )
        basis in 1:terminal_count || _fem_error(
            :results, basename(path), :basis_terminal,
            "basis terminal $basis is out of range";
            run_directory = run.path
        )
        frequency = _parse_real(T, columns[2], path, line_index, run)
        isapprox(
            frequency,
            frequencies[frequency_index];
            rtol = 16eps(T),
            atol = zero(T)
        ) || _fem_error(
            :results,
            basename(path),
            :frequency_hz,
            "frequency $frequency does not match input " *
            "$(frequencies[frequency_index]) at index $frequency_index";
            run_directory = run.path
        )
        key = (response, basis, frequency_index)
        key in seen && _fem_error(
            :results,
            basename(path),
            :indices,
            "duplicate row for response/basis/frequency $key";
            run_directory = run.path
        )
        push!(seen, key)
        real_part = _parse_real(T, columns[5], path, line_index, run)
        imaginary_part = _parse_real(T, columns[6], path, line_index, run)
        matrix[response, basis, frequency_index] = complex(real_part, imaginary_part)
    end
    length(seen) == expected || _fem_error(
        :results,
        basename(path),
        :indices,
        "one or more response/basis/frequency tuples are missing";
        run_directory = run.path
    )
    return matrix
end

function _validate_completion(
        path::String,
        frequency_count::Int,
        terminal_count::Int,
        run::FEMRun
)
    lines = _nonempty_lines(path, run)
    length(lines) == 2 || _fem_error(
        :results,
        basename(path),
        :completion,
        "completion marker must contain exactly one data row";
        run_directory = run.path
    )
    split(lines[1], '\t'; keepempty = true) == FEM_COMPLETE_HEADER || _fem_error(
        :results,
        basename(path),
        :header,
        "malformed completion header";
        run_directory = run.path
    )
    columns = split(lines[2], '\t'; keepempty = true)
    length(columns) == 6 || _fem_error(
        :results, basename(path), :completion, "malformed completion row";
        run_directory = run.path
    )
    values = [_parse_integer(column, path, 2, run) for column in columns]
    expected_rows = frequency_count * terminal_count * terminal_count
    values == [
        frequency_count,
        terminal_count,
        expected_rows,
        expected_rows,
        expected_rows,
        1
    ] || _fem_error(
        :results,
        basename(path),
        :completion,
        "completion marker reports an incomplete or failed scan: $values";
        run_directory = run.path
    )
    return nothing
end

function _expected_map_paths(run::FEMRun, frequency_count::Int, terminal_count::Int)
    return [joinpath(
                run.path,
                "maps",
                @sprintf("%s_f%04d_b%04d.pos", quantity, frequency, basis)
            )
            for frequency in 1:frequency_count
            for basis in 1:terminal_count
            for quantity in FEM_FIELD_QUANTITIES]
end

function _validate_maps(
        run::FEMRun,
        frequency_count::Int,
        terminal_count::Int,
        enabled::Bool
)
    expected = enabled ? _expected_map_paths(run, frequency_count, terminal_count) :
               String[]
    missing = filter(!isfile, expected)
    isempty(missing) || _fem_error(
        :results,
        "field_maps",
        :files,
        "missing $(length(missing)) expected field-map files";
        run_directory = run.path
    )
    if !enabled
        maps_directory = joinpath(run.path, "maps")
        emitted = isdir(maps_directory) ?
                  filter(
            path -> endswith(path, ".pos"), readdir(maps_directory; join = true)
        ) : String[]
        isempty(emitted) || _fem_error(
            :results,
            "field_maps",
            :files,
            "field maps were emitted while plot_field_maps=false";
            run_directory = run.path
        )
    end
    return expected
end

function _parse_scan(
        run::FEMRun,
        model::FEMResolvedModel{T},
        formulation::LineCableModelsFEM
) where {T <: Real}
    raw = joinpath(run.path, "raw")
    terminal_count = length(model.terminal_ids)
    frequency_count = length(model.problem.frequencies)
    _validate_completion(
        joinpath(raw, "scan_complete.tsv"),
        frequency_count,
        terminal_count,
        run
    )
    Z = _parse_raw_matrix(
        T,
        joinpath(raw, "Z.tsv"),
        model.problem.frequencies,
        terminal_count,
        run
    )
    P = _parse_raw_matrix(
        T,
        joinpath(raw, "P.tsv"),
        model.problem.frequencies,
        terminal_count,
        run
    )
    maps = _validate_maps(
        run,
        frequency_count,
        terminal_count,
        formulation.execution.plot_field_maps
    )
    return FEMScan(Z, P, maps)
end

function _line_parameters(
        run::FEMRun,
        model::FEMResolvedModel{T},
        formulation::LineCableModelsFEM,
        execution::NamedTuple,
        scan::FEMScan{T}
) where {T <: Real}
    reduced = Engine.reduce_primitive_matrices(
        scan.Z,
        scan.P,
        model.problem.system.connection_order,
        formulation.options
    )
    inversion = Engine.potential_to_admittance(reduced.P; diagnostics = true)
    Z = reduced.Z
    Y = inversion.Y
    basis = :pul
    if execution.output_basis === Val(:total)
        Z = Z .* model.problem.system.line_length
        Y = Y .* model.problem.system.line_length
        basis = :total
    end
    keep_run = formulation.execution.keep_run_directory
    record = FEMRunRecord(
        completed,
        keep_run ? run.path : nothing,
        run.mesh_source,
        run.mesh_fingerprint,
        run.getdp_invocations,
        keep_run ? scan.map_paths : String[]
    )
    trace = execution.trace === Val(true) ?
            (
        Z_primitive = scan.Z,
        P_primitive = scan.P,
        phase_map = copy(model.problem.system.connection_order)
    ) : nothing
    details = (
        fem = (
        run = record,
        terminal_ids = copy(model.terminal_ids),
        reduced_phase_map = reduced.phase_map,
        inversion_residuals = inversion.residuals,
        condition_numbers = inversion.condition_numbers,
        primitive = trace
    ),
    )
    return LineParameters(
        PhaseDomain,
        SeriesImpedance(Z; basis),
        ShuntAdmittance(Y; basis),
        model.problem.frequencies,
        details
    )
end

function _merge_maps!(paths::Vector{String})
    for path in paths
        gmsh.merge(path)
    end
    return nothing
end
