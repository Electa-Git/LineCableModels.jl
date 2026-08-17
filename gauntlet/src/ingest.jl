function _json_dict(path::AbstractString)
    return JSON3.read(read(path, String), Dict{String, Any})
end

_native_path(path::AbstractString) = replace(path, '\\' => '/')

function _sha256(path::AbstractString)
    return open(path, "r") do io
        bytes2hex(SHA.sha256(io))
    end
end

function _verify(path::AbstractString, expected::AbstractString)
    isfile(path) || throw(ArgumentError("manifest artifact is missing: $path"))
    observed = _sha256(path)
    observed == lowercase(expected) || throw(ArgumentError(
        "artifact checksum mismatch for $path: expected $expected, got $observed",
    ))
    return path
end

function _block_dict(block::PSCADIO.Block)
    return Dict{String, Any}(
        "name" => block.name,
        "value" => block.value,
        "fields" => Dict{String, Any}(block.fields),
        "children" => _block_dict.(block.children)
    )
end

function _line_dict(parameters)
    parameters === nothing && return nothing
    return Dict{String, Any}(
        "frequency" => collect(frequencies(parameters)),
        "Z" => Array(Z(parameters)),
        "Y" => Array(Y(parameters)),
        "basis" => String(basis(parameters))
    )
end

function _fit_dict(fit::Nothing)
    return nothing
end

function _fit_dict(fit::Fit)
    return Dict{String, Any}(
        "columns" => [Dict{String, Any}(
             "poles" => column.poles,
             "residues" => column.residues
         ) for column in fit.columns],
        "constant" => fit.constant,
        "groups" => [Dict{String, Any}(
             "delay" => group.delay,
             "poles" => group.poles,
             "residues" => group.residues
         ) for group in fit.groups],
        "frequency_range" => collect(fit.frequency_range)
    )
end

function _modes_dict(modes::Nothing)
    return nothing
end

function _modes_dict(modes::Modes)
    return Dict{String, Any}(
        "frequency" => modes.frequency,
        "transform" => modes.transform,
        "propagation" => modes.propagation,
        "characteristic" => modes.characteristic,
        "fitted_characteristic" => modes.fitted_characteristic,
        "phase_propagation" => modes.phase_propagation,
        "modal_propagation" => modes.modal_propagation,
        "fitted_phase_propagation" => modes.fitted_phase_propagation
    )
end

function _terminal_dict(terminal::Nothing)
    return nothing
end

function _terminal_dict(terminal::Terminal)
    return Dict{String, Any}(
        "frequency" => terminal.frequency,
        "open" => terminal.open,
        "short" => terminal.short,
        "empty" => String.(terminal.empty)
    )
end

function _port_dict(port::Port)
    return Dict{String, Any}(
        "id" => port.id,
        "cable" => port.cable,
        "component" => port.component,
        "phase" => port.phase
    )
end

function _provenance_dict(provenance::Provenance)
    return Dict{String, Any}(
        "source" => String(provenance.source),
        "version" => provenance.version,
        "campaign" => provenance.campaign,
        "definition" => provenance.definition,
        "source_sha256" => provenance.source_sha256,
        "specification_sha256" => provenance.specification_sha256,
        "case_sha256" => provenance.case_sha256,
        "elapsed_seconds" => provenance.elapsed_seconds,
        "artifact_sha256" => provenance.artifact_sha256,
        "cohort" => provenance.cohort
    )
end

function _assumption_dict(assumption::Assumption)
    return Dict{String, Any}(
        "subject" => String(assumption.subject),
        "detail" => assumption.detail
    )
end

function _artifact_paths(case_root, manifest)
    result = Dict{String, String}()
    hashes = Dict{String, String}()
    for artifact in manifest["artifacts"]
        relative = _native_path(String(artifact["path"]))
        path = _verify(joinpath(case_root, relative), String(artifact["sha256"]))
        name = basename(path)
        result[name] = path
        hashes[name] = String(artifact["sha256"])
    end
    return result, hashes
end

function _one_file(files, extensions)
    matches = sort!([path
                     for (name, path) in files
                     if lowercase(splitext(name)[2]) in extensions])
    isempty(matches) && return nothing
    return first(matches)
end

function _suffix_file(files, suffix)
    matches = [path for (name, path) in files if endswith(lowercase(name), suffix)]
    isempty(matches) && return nothing
    length(matches) == 1 || throw(ArgumentError("multiple PSCAD files end with $suffix"))
    return only(matches)
end

function _polar(files, magnitude_suffix, phase_suffix; matrix = true)
    magnitude_path = _suffix_file(files, magnitude_suffix)
    phase_path = _suffix_file(files, phase_suffix)
    magnitude_path === nothing && phase_path === nothing && return nothing, nothing
    (magnitude_path === nothing) == (phase_path === nothing) || throw(ArgumentError(
        "PSCAD magnitude/phase pair is incomplete: $magnitude_suffix, $phase_suffix",
    ))
    magnitude = PSCADIO.read_detailed(magnitude_path; allow_empty = true)
    phase = PSCADIO.read_detailed(phase_path; allow_empty = true)
    magnitude === nothing && phase === nothing && return nothing, nothing
    (magnitude === nothing) == (phase === nothing) || throw(ArgumentError(
        "one PSCAD detailed polar channel is empty: $magnitude_suffix, $phase_suffix",
    ))
    return magnitude, PSCADIO.combine_polar(magnitude, phase; matrix)
end

function _polar_pairs(files, magnitude_suffix, phase_suffix; matrix = true)
    magnitude_path = _suffix_file(files, magnitude_suffix)
    phase_path = _suffix_file(files, phase_suffix)
    magnitude_path === nothing && phase_path === nothing && return nothing, nothing, nothing
    (magnitude_path === nothing) == (phase_path === nothing) || throw(ArgumentError(
        "PSCAD calculated/fitted pair is incomplete: $magnitude_suffix, $phase_suffix",
    ))
    magnitude = PSCADIO.read_detailed(magnitude_path; allow_empty = true)
    phase = PSCADIO.read_detailed(phase_path; allow_empty = true)
    magnitude === nothing && phase === nothing && return nothing, nothing, nothing
    (magnitude === nothing) == (phase === nothing) || throw(ArgumentError(
        "one PSCAD calculated/fitted channel is empty: $magnitude_suffix, $phase_suffix",
    ))
    calculated, fitted = PSCADIO.combine_polar_pairs(magnitude, phase; matrix)
    return magnitude, calculated, fitted
end

function _ordinary(files)
    candidates = [path
                  for (name, path) in files
                  if lowercase(splitext(name)[2]) == ".out" &&
        !occursin(
            r"_(zm|zp|ym|yp|timp|tipp|lamdam|lamdap|ycmp|ycpp|hmp|hpp|hmm|hpm|gdio|gdis|open|short|zs|ys)\.out$"i,
            name)]
    isempty(candidates) && return nothing
    return PSCADIO.read_ordinary(first(sort!(candidates)))
end

function _phase_reference(files)
    z_table, impedance = _polar(files, "_zm.out", "_zp.out")
    y_table, admittance = _polar(files, "_ym.out", "_yp.out")
    if impedance !== nothing || admittance !== nothing
        (impedance === nothing) == (admittance === nothing) || throw(ArgumentError(
            "detailed PSCAD phase output requires both Z and Y",
        ))
        z_table.frequency == y_table.frequency || throw(DimensionMismatch(
            "detailed PSCAD Z and Y frequencies differ",
        ))
        return EN.LineParameters(
            impedance, admittance, z_table.frequency; basis = :per_length
        ),
        z_table.frequency
    end
    ordinary = _ordinary(files)
    ordinary === nothing && return nothing, Float64[]
    (ordinary.Z === nothing) == (ordinary.Y === nothing) || throw(ArgumentError(
        "ordinary PSCAD phase output requires both Z and Y",
    ))
    ordinary.Z === nothing && return nothing, Float64[]
    z = reshape(ordinary.Z, size(ordinary.Z)..., 1)
    y = reshape(ordinary.Y, size(ordinary.Y)..., 1)
    return EN.LineParameters(z, y, [ordinary.frequency]; basis = :per_length),
    [ordinary.frequency]
end

function _sequence_reference(files)
    ordinary = _ordinary(files)
    ordinary === nothing && return nothing, nothing
    (ordinary.sequence_Z === nothing) == (ordinary.sequence_Y === nothing) ||
        throw(ArgumentError("ordinary PSCAD sequence output requires both Z and Y"))
    ordinary.sequence_Z === nothing && return nothing, ordinary.sequence_transform
    z = reshape(ordinary.sequence_Z, size(ordinary.sequence_Z)..., 1)
    y = reshape(ordinary.sequence_Y, size(ordinary.sequence_Y)..., 1)
    parameters = EN.LineParameters(
        LineCableModels.ModalDomain,
        z,
        y,
        [ordinary.frequency];
        basis = :per_length
    )
    return parameters, ordinary.sequence_transform
end

function _modal_reference(files)
    transform_table, transform = _polar(files, "_timp.out", "_tipp.out")
    lambda_table, propagation = _polar(files, "_lamdam.out", "_lamdap.out"; matrix = false)
    yc_table, characteristic, fitted_characteristic = _polar_pairs(
        files, "_ycmp.out", "_ycpp.out"
    )
    h_table, phase_propagation, fitted_phase_propagation = _polar_pairs(
        files, "_hmp.out", "_hpp.out"
    )
    hm_table, modal_propagation, _ = _polar_pairs(
        files, "_hmm.out", "_hpm.out"; matrix = false
    )
    tables = filter(!isnothing, (
        transform_table, lambda_table, yc_table, h_table, hm_table
    ))
    isempty(tables) && return nothing
    frequency = first(tables).frequency
    all(table -> table.frequency == frequency, tables) || throw(DimensionMismatch(
        "PSCAD modal output frequencies differ",
    ))
    return Modes(
        frequency;
        transform,
        propagation,
        characteristic,
        fitted_characteristic,
        phase_propagation,
        modal_propagation,
        fitted_phase_propagation
    )
end

function _empty_outputs(files)
    result = Symbol[]
    for name in (:open, :short, :zs, :ys)
        path = _suffix_file(files, "_$(name).out")
        path === nothing && continue
        table = PSCADIO.read_detailed(path; allow_empty = true)
        table === nothing && push!(result, name)
    end
    return result
end

function _terminal_reference(files)
    open_path = _suffix_file(files, "_gdio.out")
    short_path = _suffix_file(files, "_gdis.out")
    empty = _empty_outputs(files)
    open_path === nothing && short_path === nothing && isempty(empty) && return nothing
    open_table = open_path === nothing ? nothing :
                 PSCADIO.read_detailed(open_path; allow_empty = true)
    short_table = short_path === nothing ? nothing :
                  PSCADIO.read_detailed(short_path; allow_empty = true)
    tables = filter(!isnothing, (open_table, short_table))
    frequency = isempty(tables) ? Float64[] : first(tables).frequency
    all(table -> table.frequency == frequency, tables) || throw(DimensionMismatch(
        "PSCAD open- and short-circuit frequencies differ",
    ))
    return Terminal(
        frequency;
        open = open_table === nothing ? nothing : PSCADIO.terminal_values(open_table),
        short = short_table === nothing ? nothing : PSCADIO.terminal_values(short_table),
        empty
    )
end

function _fit_reference(files, input)
    path = _one_file(files, (".clo", ".tlo"))
    path === nothing && return nothing
    occursin("N.o. conductors", read(path, String)) || return nothing
    options = only(PSCADIO.block(input.root, "Frequency Dep. (Phase) Model Options"))
    range = (
        Float64(_field(options, "Curve Fitting Start Frequency")),
        Float64(_field(options, "Curve Fitting End Frequency"))
    )
    return PSCADIO.read_fit(path; frequency_range = range)
end

function _is_prefix(short_path::AbstractString, long_path::AbstractString)
    filesize(long_path) > filesize(short_path) || return false
    open(short_path, "r") do short
        open(long_path, "r") do long
            while !eof(short)
                read(short, UInt8) == read(long, UInt8) || return false
            end
        end
    end
    return true
end

function _fit_candidates(case_root::AbstractString, artifact::AbstractString)
    name = basename(artifact)
    candidates = String[]
    for (root, _, names) in walkdir(case_root), candidate_name in names

        candidate_name == name || continue
        candidate = joinpath(root, candidate_name)
        candidate == artifact && continue
        _is_prefix(artifact, candidate) || continue
        try
            PSCADIO.read_fit(candidate)
            push!(candidates, candidate)
        catch
            # A longer file is useful only when the complete fit grammar parses.
        end
    end
    return candidates
end

function _recover_truncated_fit!(files, hashes, case_root, assumptions)
    path = _one_file(files, (".clo", ".tlo"))
    path === nothing && return false
    occursin("N.o. conductors", read(path, String)) || return false
    try
        PSCADIO.read_fit(path)
        return false
    catch error
        error isa EOFError || rethrow()
    end
    candidates = _fit_candidates(case_root, path)
    unique_hashes = unique(_sha256.(candidates))
    length(unique_hashes) == 1 || throw(ArgumentError(
        "truncated PSCAD fit $(basename(path)) has $(length(unique_hashes)) " *
        "distinct complete recovery candidates",
    ))
    recovered = only(filter(candidate -> _sha256(candidate) == only(unique_hashes), candidates))
    name = basename(path)
    hashes["$name.manifest-truncated"] = hashes[name]
    hashes[name] = _sha256(recovered)
    files[name] = recovered
    push!(assumptions,
        Assumption(
            :artifact_recovery,
            "The manifest fit artifact was truncated. Gauntlet recovered the complete " *
            "same-run file whose bytes begin with the complete manifest artifact; both " *
            "source hashes are retained."
        ))
    return true
end

function _reference(files, input, ports)
    phase, _ = _phase_reference(files)
    sequence, sequence_transform = _sequence_reference(files)
    modes = _modal_reference(files)
    fit = _fit_reference(files, input)
    terminal = _terminal_reference(files)
    matrix_size = phase === nothing ?
                  modes === nothing || modes.transform === nothing ? length(ports) :
                  size(modes.transform, 1) : size(Z(phase), 1)
    selected_ports = if length(ports) == matrix_size
        ports
    else
        active = filter(port -> port.phase > 0, ports)
        length(active) == matrix_size || throw(DimensionMismatch(
            "PSCAD matrix size $matrix_size cannot be aligned with $(length(ports)) " *
            "native ports ($(length(active)) active)",
        ))
        active
    end
    return Reference(;
        phase,
        sequence,
        sequence_transform,
        modes,
        fit,
        terminal,
        ports = selected_ports
    )
end

function _copy_raw(files, destination)
    mkpath(destination)
    copied = Dict{String, String}()
    for (name, source) in sort!(collect(files); by = first)
        lowercase(splitext(name)[2]) in (".cli", ".tli", ".out", ".clo", ".tlo") ||
            continue
        target = joinpath(destination, name)
        cp(source, target; force = true)
        copied[name] = _sha256(target)
    end
    return copied
end

function _case_toml(path, payload)
    open(path, "w") do io
        TOML.print(io, payload; sorted = true)
    end
    return path
end

function _record_provenance(record, manifest, hashes)
    formulations = record["variant"]["formulations"]
    cohort_parts = String[
        String(record["source_sha256"]),
        String(record["line"]["definition"]),
        string(record["variant"]["resistivity_ohm_m"])
    ]
    append!(cohort_parts, (
        "$(key)=$(formulations[key])" for key in sort!(collect(keys(formulations)))
    ))
    cohort = bytes2hex(SHA.sha256(join(cohort_parts, '\0')))
    return Provenance(
        :pscad,
        String(record["pscad_version"]),
        "PSCAD-reference-v1",
        String(record["line"]["definition"]),
        String(record["source_sha256"]),
        String(record["specification_sha256"]),
        String(record["case_id"]),
        get(manifest, "elapsed_seconds", nothing),
        hashes,
        cohort
    )
end

function _canonical_records(dataset)
    canonical_paths = Dict{String, String}()
    aliases = Dict{String, String}()
    for group in dataset["duplicate_effective_case_groups"]
        id = String(group["case_id"])
        canonical_paths[id] = _native_path(String(group["canonical"]))
        for alias in group["aliases"]
            aliases[_native_path(String(alias))] = id
        end
    end
    records = Dict{String, Dict{String, Any}}()
    for raw_record in dataset["records"]
        record = Dict{String, Any}(String(key) => value for (key, value) in raw_record)
        id = String(record["case_id"])
        path = _native_path(String(record["path"]))
        if !haskey(records, id) || get(canonical_paths, id, path) == path
            records[id] = record
        end
    end
    return records, aliases
end

function _write_success(
        record,
        manifest,
        files,
        hashes,
        case_root,
        raw_root,
        normalized_root
)
    id = String(record["case_id"])
    input_path = _one_file(files, (".cli", ".tli"))
    input_path === nothing && throw(ArgumentError("successful case $id has no line input"))
    input = PSCADIO.read_input(input_path)
    family = _family(input)
    reduction = _case_reduction(record["variant"]["reduction"])
    problem, ports = materialize(family, input, reduction, id)
    assumptions_value = _assumptions(input, family)
    recovered = _recover_truncated_fit!(files, hashes, case_root, assumptions_value)
    reference_value = _reference(files, input, ports)
    earth_name = String(get(record["variant"]["formulations"], "EarthForm", "UNKNOWN"))
    earth_name == "DIRECT_NUMERICAL_INTEGRATION" || push!(
        assumptions_value,
        Assumption(
            :formulation,
            "PSCAD $earth_name has no like-named native implementation and is diagnostic."
        )
    )
    fidelity = earth_name == "DIRECT_NUMERICAL_INTEGRATION" ?
               _fidelity(family, input) : Approximate()
    provenance_value = _record_provenance(record, manifest, hashes)

    raw_case = joinpath(raw_root, "cases", id)
    copied_hashes = _copy_raw(files, raw_case)
    raw_record = joinpath(raw_case, "case.toml")
    _case_toml(raw_record,
        Dict(
            "schema_version" => 1,
            "id" => id,
            "family" => lowercase(String(nameof(typeof(family)))),
            "fidelity" => lowercase(String(nameof(typeof(fidelity)))),
            "reduction" => lowercase(String(nameof(typeof(reduction)))),
            "source_path" => String(record["path"]),
            "source" => String(manifest["source"]),
            "line" => manifest["line"],
            "variant" => manifest["variant"],
            "mutations" => manifest["mutations"],
            "provenance" => _provenance_dict(provenance_value),
            "artifact_sha256" => copied_hashes
        ))

    record_payload = Dict{String, Any}(
        "schema_version" => 1,
        "id" => id,
        "family" => lowercase(String(nameof(typeof(family)))),
        "fidelity" => lowercase(String(nameof(typeof(fidelity)))),
        "reduction" => lowercase(String(nameof(typeof(reduction)))),
        "input_kind" => String(input.kind),
        "input" => _block_dict(input.root),
        "reference" => Dict{String, Any}(
            "phase" => _line_dict(reference_value.phase),
            "sequence" => _line_dict(reference_value.sequence),
            "sequence_transform" => reference_value.sequence_transform,
            "modes" => _modes_dict(reference_value.modes),
            "fit" => _fit_dict(reference_value.fit),
            "terminal" => _terminal_dict(reference_value.terminal),
            "ports" => _port_dict.(reference_value.ports)
        ),
        "assumptions" => _assumption_dict.(assumptions_value),
        "provenance" => _provenance_dict(provenance_value)
    )
    target = joinpath(normalized_root, "cases", "$id.jld2")
    mkpath(dirname(target))
    JLD2.jldopen(target, "w"; compress = true) do file
        file["record"] = record_payload
    end
    # Loading the complete object here proves that each normalized successful
    # case still materializes through the current LineCableModels API.
    _load_success(target)
    return (; target, raw_record, recovered)
end

function _trimmed_messages(manifest)
    messages = String[]
    for key in ("pscad_messages", "pscad_output")
        value = get(manifest, key, nothing)
        value === nothing && continue
        text = value isa AbstractString ? String(value) : JSON3.write(value)
        lines = split(text, '\n')
        append!(messages, last(lines, min(40, length(lines))))
    end
    return join(filter(!isempty, strip.(messages)), '\n')
end

function _rejection_family(files)
    input_path = _one_file(files, (".cli", ".tli"))
    input_path === nothing && return Coax(), nothing
    input = PSCADIO.read_input(input_path)
    return _family(input), input_path
end

function _write_rejection(record, manifest, files, hashes, raw_root, normalized_root)
    id = String(record["case_id"])
    family, input_path = _rejection_family(files)
    reduction = _case_reduction(record["variant"]["reduction"])
    provenance_value = _record_provenance(record, manifest, hashes)
    payload = Dict{String, Any}(
        "schema_version" => 1,
        "id" => id,
        "family" => lowercase(String(nameof(typeof(family)))),
        "reduction" => lowercase(String(nameof(typeof(reduction)))),
        "reason" => _trimmed_messages(manifest),
        "source" => String(manifest["source"]),
        "line" => manifest["line"],
        "variant" => manifest["variant"],
        "mutations" => manifest["mutations"],
        "assumptions" => Any[],
        "provenance" => _provenance_dict(provenance_value)
    )
    raw_case = joinpath(raw_root, "rejections", id)
    mkpath(raw_case)
    input_path === nothing ||
        cp(input_path, joinpath(raw_case, basename(input_path)); force = true)
    raw_record = _case_toml(joinpath(raw_case, "rejection.toml"), payload)
    target = joinpath(normalized_root, "rejections", "$id.toml")
    mkpath(dirname(target))
    _case_toml(target, payload)
    return (; target, raw_record)
end

function _write_schema(root)
    _case_toml(joinpath(root, "schema.toml"),
        Dict(
            "schema_version" => 1,
            "frequency_unit" => "Hz",
            "impedance_unit" => "ohm/m",
            "admittance_unit" => "S/m",
            "angle_unit" => "rad",
            "complex_storage" => "ComplexF64 arrays",
            "missing_channels" => "absent fields"
        ))
end

"""Ingest an authoritative PSCAD campaign into raw and normalized corpora."""
function ingest(
        ::PSCAD,
        source::AbstractString,
        destination::AbstractString;
        ids = nothing
)
    source_root = abspath(source)
    destination_root = abspath(destination)
    dataset_path = joinpath(source_root, "dataset_manifest.json")
    isfile(dataset_path) ||
        throw(ArgumentError("missing dataset_manifest.json in $source_root"))
    dataset = _json_dict(dataset_path)
    get(dataset, "schema_version", 0) == 1 || throw(ArgumentError(
        "unsupported PSCAD dataset schema $(get(dataset, "schema_version", nothing))",
    ))
    records, aliases = _canonical_records(dataset)
    length(records) == Int(dataset["unique_effective_case_count"]) ||
        throw(DimensionMismatch("canonical record count differs from the source manifest"))
    full_campaign = ids === nothing
    if !full_campaign
        selected = Set(String.(ids))
        missing_ids = setdiff(selected, keys(records))
        isempty(missing_ids) || throw(KeyError(
            "selected case identifiers are absent from the source: $(join(sort!(collect(missing_ids)), ", "))",
        ))
        records = Dict(id => records[id] for id in selected)
        aliases = Dict(alias => target for (alias, target) in aliases if target in selected)
    end

    raw_root = joinpath(destination_root, "raw")
    normalized_root = joinpath(destination_root, "normalized")
    ispath(raw_root) && throw(ArgumentError("raw destination already exists: $raw_root"))
    ispath(normalized_root) && throw(ArgumentError(
        "normalized destination already exists: $normalized_root",
    ))
    mkpath(raw_root)
    mkpath(normalized_root)

    success_index = Dict{String, String}()
    rejection_index = Dict{String, String}()
    raw_success_index = Dict{String, String}()
    raw_rejection_index = Dict{String, String}()
    family_counts = Dict{String, Int}()
    detailed = 0
    sparse = 0
    recovered_fits = 0
    failures = Dict{String, String}()
    ordered_ids = sort!(collect(keys(records)))
    for (position, id) in enumerate(ordered_ids)
        (position == 1 || position % 25 == 0) &&
            @info("Ingesting PSCAD corpus", case=position, total=length(ordered_ids), id)
        record = records[id]
        try
            case_root = joinpath(source_root, _native_path(String(record["path"])))
            manifest_path = joinpath(case_root, "manifest.json")
            manifest = _json_dict(manifest_path)
            String(manifest["case_id"]) == id || throw(ArgumentError(
                "case identifier differs between dataset and case manifest: $id",
            ))
            files, hashes = _artifact_paths(case_root, manifest)
            status = String(record["status"])
            if status == "success"
                written = _write_success(
                    record,
                    manifest,
                    files,
                    hashes,
                    case_root,
                    raw_root,
                    normalized_root
                )
                success_index[id] = relpath(written.target, normalized_root)
                raw_success_index[id] = relpath(written.raw_record, raw_root)
                recovered_fits += written.recovered
                input_path = _one_file(files, (".cli", ".tli"))
                kind = String(PSCADIO.read_input(input_path).kind)
                family_counts[kind] = get(family_counts, kind, 0) + 1
                if _suffix_file(files, "_zm.out") === nothing
                    sparse += 1
                else
                    detailed += 1
                end
            elseif status == "solver_rejected"
                written = _write_rejection(
                    record, manifest, files, hashes, raw_root, normalized_root
                )
                rejection_index[id] = relpath(written.target, normalized_root)
                raw_rejection_index[id] = relpath(written.raw_record, raw_root)
            else
                throw(ArgumentError("unknown PSCAD case status '$status'"))
            end
        catch error
            failures[id] = sprint(showerror, error)
        end
    end
    if !isempty(failures)
        _case_toml(joinpath(destination_root, "ingestion_failures.toml"),
            Dict(
                "schema_version" => 1,
                "failures" => failures
            ))
        throw(ErrorException(
            "PSCAD ingestion rejected $(length(failures)) canonical cases; " *
            "see $(joinpath(destination_root, "ingestion_failures.toml"))",
        ))
    end

    expected_success = full_campaign ? 869 :
                       count(
        record -> String(record["status"]) == "success", values(records)
    )
    expected_rejections = full_campaign ? 199 :
                          count(
        record -> String(record["status"]) == "solver_rejected", values(records)
    )
    length(success_index) == expected_success || throw(DimensionMismatch(
        "expected $expected_success canonical successes, found $(length(success_index))",
    ))
    length(rejection_index) == expected_rejections || throw(DimensionMismatch(
        "expected $expected_rejections rejections, found $(length(rejection_index))",
    ))
    if full_campaign
        detailed == 689 || throw(DimensionMismatch(
            "expected 689 detailed cases, found $detailed",
        ))
        sparse == 180 || throw(DimensionMismatch(
            "expected 180 sparse cases, found $sparse",
        ))
        family_counts == Dict(
            "coax" => 288,
            "overhead" => 240,
            "mixed" => 269,
            "pipe" => 72
        ) || throw(DimensionMismatch(
            "successful family counts differ from the authoritative campaign",
        ))
    end

    index = Dict{String, Any}(
        "schema_version" => 1,
        "source_manifest_sha256" => _sha256(dataset_path),
        "cases" => success_index,
        "rejections" => rejection_index,
        "excluded" => dataset["excluded_definitions"],
        "family_counts" => family_counts,
        "detailed_cases" => detailed,
        "sparse_cases" => sparse,
        "recovered_fits" => recovered_fits
    )
    _case_toml(joinpath(normalized_root, "index.toml"), index)
    raw_index = copy(index)
    raw_index["artifact_kind"] = "raw"
    raw_index["cases"] = raw_success_index
    raw_index["rejections"] = raw_rejection_index
    _case_toml(joinpath(raw_root, "index.toml"), raw_index)
    alias_payload = Dict("aliases" => aliases)
    _case_toml(joinpath(normalized_root, "aliases.toml"), alias_payload)
    _case_toml(joinpath(raw_root, "aliases.toml"), alias_payload)
    _write_schema(normalized_root)
    for name in ("dataset_manifest.json", "benchmark.yml", "official_sources.yml", "README.md")
        source_path = joinpath(source_root, name)
        isfile(source_path) || continue
        cp(source_path, joinpath(raw_root, name); force = true)
    end

    # Final load proves the index and every lazy path are self-consistent.
    corpus = Corpus(normalized_root)
    length(corpus) == length(records) || throw(DimensionMismatch(
        "normalized corpus length differs from the canonical source",
    ))
    return (;
        raw = raw_root,
        normalized = normalized_root,
        successes = length(success_index),
        rejections = length(rejection_index),
        aliases = length(aliases),
        detailed,
        sparse,
        recovered_fits,
        families = family_counts
    )
end
