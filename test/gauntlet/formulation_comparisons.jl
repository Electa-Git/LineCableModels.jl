include(joinpath(@__DIR__, "fem_catalogue.jl"))

const COMPARISON_SCHEMA_VERSION = 2
const PSCAD_REFERENCE_ROOT = joinpath(
    pkgdir(LineCableModels),
    ".linecablemodels",
    "gauntlet",
    "references",
    "pscad"
)
const COMPARISON_ROOT = joinpath(
    pkgdir(LineCableModels),
    ".linecablemodels",
    "gauntlet",
    "comparisons",
    "formulation_corpus"
)

function tree_digest(root)
    paths = sort!([joinpath(directory, file)
                   for (directory, _, files) in walkdir(root)
                   for file in files
                   if endswith(file, ".jl")])
    contents = UInt8[]
    for path in paths
        append!(contents, read(path))
        push!(contents, 0x00)
    end
    return bytes2hex(sha256(contents))
end

const ENGINE_SHA256 = tree_digest(joinpath(pkgdir(LineCableModels), "src"))

function file_digest(paths)
    contents = UInt8[]
    for path in sort!(abspath.(collect(paths)))
        append!(contents, read(path))
        push!(contents, 0x00)
    end
    return bytes2hex(sha256(contents))
end

const COMPARISON_SOURCE_SHA256 = file_digest((
    @__FILE__,
    joinpath(@__DIR__, "fem_catalogue.jl"),
    joinpath(@__DIR__, "reference_grid.jl")
))

function reference_parameters(document)
    LineParameters(
        PhaseDomain,
        document["Z"],
        document["Y"],
        document["frequencies"];
        basis = :pul
    )
end

function reference_record(backend, formula, field, path, document)
    return (;
        backend,
        formula,
        field,
        path,
        sha256 = bytes2hex(sha256(read(path))),
        parameters = reference_parameters(document)
    )
end

function references(case_id, model)
    records = NamedTuple[]
    fem_path, fem_document, _ = load_reference(case_id, model)
    push!(records, reference_record(
        :fem,
        :LineCableModelsFEM,
        :all,
        fem_path,
        fem_document
    ))

    root = joinpath(PSCAD_REFERENCE_ROOT, string(case_id))
    isdir(root) || error("PSCAD references are missing for $case_id: $root")
    for directory in sort!(filter(isdir, readdir(root; join = true)))
        path = joinpath(directory, "reference.jld2")
        isfile(path) || continue
        document = JLD2.load(path)
        document["schema_version"] == 1 || error(
            "unsupported PSCAD reference schema at $path",
        )
        document["status"] === :complete || error(
            "incomplete PSCAD reference at $path",
        )
        document["case_source_sha256"] == model.source_sha256 || error(
            "$case_id PSCAD reference case digest is stale: $path",
        )
        document["frequencies"] == model.problem.frequencies || error(
            "$case_id PSCAD reference frequency grid differs: $path",
        )
        document["port_order"] == model.port_order || error(
            "$case_id PSCAD reference terminal order differs: $path",
        )
        push!(records,
            reference_record(
                :pscad,
                Symbol(document["formula"]),
                Symbol(document["field"]),
                path,
                document
            ))
    end
    count(record -> record.backend === :pscad, records) > 0 || error(
        "no complete PSCAD formulation reference was found for $case_id",
    )
    return records
end

function comparison_path(case_id, selected)
    joinpath(
        COMPARISON_ROOT,
        "cases",
        string(case_id),
        string(selected.id) * ".jld2"
    )
end

function comparison_failure_path(case_id, selected)
    return joinpath(
        dirname(comparison_path(case_id, selected)),
        string(selected.id) * "_failure.jld2"
    )
end

function clear_comparison_failure(case_id, selected)
    path = comparison_failure_path(case_id, selected)
    isfile(path) && rm(path)
    hash_path = path * ".sha256"
    isfile(hash_path) && rm(hash_path)
    return nothing
end

function comparison_hashes(records)
    sort!([record.sha256 for record in records])
end

function valid_comparison(path, model, selected, records)
    isfile(path) || return false
    try
        document = JLD2.load(path)
        return document["schema_version"] == COMPARISON_SCHEMA_VERSION &&
               document["status"] === :complete &&
               document["case_source_sha256"] == model.source_sha256 &&
               document["engine_sha256"] == ENGINE_SHA256 &&
               document["comparison_source_sha256"] ==
               COMPARISON_SOURCE_SHA256 &&
               document["frequencies"] == model.problem.frequencies &&
               document["port_order"] == model.port_order &&
               document["variant_id"] == selected.id &&
               document["formula"] == selected.identifier &&
               document["reference_sha256"] == comparison_hashes(records)
    catch
        return false
    end
end

function comparison_record(record, candidate)
    return (;
        reference_backend = record.backend,
        reference_formula = record.formula,
        reference_field = record.field,
        reference_path = record.path,
        reference_sha256 = record.sha256,
        errors = LineCableModels.Engine.compare(record.parameters, candidate)
    )
end

function record_comparison(
        case_id,
        model,
        selected,
        candidate,
        comparisons,
        elapsed
)
    path = comparison_path(case_id, selected)
    mkpath(dirname(path))
    temporary = path * ".new"
    JLD2.jldsave(
        temporary;
        schema_version = COMPARISON_SCHEMA_VERSION,
        kind = :cross_backend_formulation_comparison,
        status = :complete,
        case_id = string(case_id),
        case_source_sha256 = model.source_sha256,
        engine_sha256 = ENGINE_SHA256,
        comparison_source_sha256 = COMPARISON_SOURCE_SHA256,
        frequencies = copy(candidate.f),
        port_order = copy(model.port_order),
        basis = :pul,
        domain = :PhaseDomain,
        variant_id = selected.id,
        variant_kind = selected.kind,
        formula = selected.identifier,
        earth_impedance = selected.earth_impedance,
        earth_admittance = selected.earth_admittance,
        Z = copy(candidate.Z.values),
        Y = copy(candidate.Y.values),
        comparisons,
        reference_sha256 = sort!([comparison.reference_sha256 for comparison in comparisons]),
        elapsed_seconds = elapsed,
        recorded_at_utc = string(now(UTC))
    )
    mv(temporary, path; force = true)
    write_hash(path)
    return path
end

function record_comparison_failure(case_id, model, selected, exception, elapsed)
    path = comparison_failure_path(case_id, selected)
    mkpath(dirname(path))
    temporary = path * ".new"
    JLD2.jldsave(
        temporary;
        schema_version = COMPARISON_SCHEMA_VERSION,
        kind = :cross_backend_formulation_failure,
        status = :execution_failure,
        case_id = string(case_id),
        case_source_sha256 = model.source_sha256,
        engine_sha256 = ENGINE_SHA256,
        comparison_source_sha256 = COMPARISON_SOURCE_SHA256,
        frequencies = copy(model.problem.frequencies),
        port_order = copy(model.port_order),
        variant_id = selected.id,
        variant_kind = selected.kind,
        formula = selected.identifier,
        exception_type = string(typeof(exception)),
        message = sprint(showerror, exception),
        elapsed_seconds = elapsed,
        recorded_at_utc = string(now(UTC))
    )
    mv(temporary, path; force = true)
    write_hash(path)
    return path
end

function element_rows(case_id, selected, comparisons, path)
    rows = NamedTuple[]
    for comparison in comparisons
        errors = comparison.errors
        for row in axes(errors.Z.absolute, 1), column in axes(errors.Z.absolute, 2)

            push!(rows,
                (;
                    case = string(case_id),
                    reference_backend = string(comparison.reference_backend),
                    reference_formula = string(comparison.reference_formula),
                    reference_field = string(comparison.reference_field),
                    candidate_kind = string(selected.kind),
                    candidate_formula = string(selected.identifier),
                    row,
                    column,
                    Z_rms_absolute = errors.Z.absolute[row, column],
                    Z_rms_relative = errors.Z.relative[row, column],
                    Y_rms_absolute = errors.Y.absolute[row, column],
                    Y_rms_relative = errors.Y.relative[row, column],
                    artifact = path
                ))
        end
    end
    return rows
end

function table(path, rows, names)
    temporary = path * ".new"
    mkpath(dirname(path))
    open(temporary, "w") do io
        println(io, join(string.(names), '\t'))
        for row in rows
            values = replace.(
                string.((getproperty(row, name) for name in names)),
                '\t' => ' ',
                '\n' => ' '
            )
            println(io, join(values, '\t'))
        end
    end
    mv(temporary, path; force = true)
    write_hash(path)
    return path
end

function comparison_main(args = ARGS)
    coverage = catalogue()
    variants = variant.(coverage)
    available = sort!(collect(keys(case_index())); by = string)
    requested = isempty(args) ? available : Symbol.(args)
    all(in(available), requested) || error(
        "unknown case requested; available cases are $(join(available, ", "))",
    )

    rows = NamedTuple[]
    skipped = NamedTuple[]
    failures = NamedTuple[]
    println(
        "PLAN\tcases=", length(requested),
        "\tanalytical_variants=", length(variants),
        "\treferences=fem+all_applicable_pscad",
        "\tmetric=elementwise_rms_across_frequency",
        "\tmonte_carlo=false"
    )
    for case_id in requested
        model = reference_case(case_id)
        case_references = references(case_id, model)
        println(
            "BEGIN_CASE\t", case_id,
            "\tfrequencies=", length(model.problem.frequencies),
            "\tterminals=", length(model.port_order),
            "\treferences=", length(case_references)
        )
        for selected in variants
            reason = case_skip_reason(model, selected)
            if reason !== nothing
                push!(skipped,
                    (;
                        case = string(case_id),
                        candidate_kind = string(selected.kind),
                        candidate_formula = string(selected.identifier),
                        reason
                    ))
                println("SKIP\t", case_id, "\t", selected.id, "\t", reason)
                continue
            end

            path = comparison_path(case_id, selected)
            comparisons = nothing
            elapsed = 0.0
            if valid_comparison(path, model, selected, case_references)
                document = JLD2.load(path)
                comparisons = document["comparisons"]
                elapsed = document["elapsed_seconds"]
                write_hash(path)
                clear_comparison_failure(case_id, selected)
                println("REUSE\t", case_id, "\t", selected.id)
            else
                isfile(path) && archive_stale_artifact(path)
                started = time()
                try
                    candidate = compute(model.problem, formulation(selected))
                    candidate.f == model.problem.frequencies || error(
                        "$case_id $(selected.id) changed the reference frequency grid",
                    )
                    size(candidate.Z.values) == model.expected_size || error(
                        "$case_id $(selected.id) returned Z size $(size(candidate.Z.values)); expected $(model.expected_size)",
                    )
                    comparisons = comparison_record.(case_references, Ref(candidate))
                    elapsed = time() - started
                    path = record_comparison(
                        case_id,
                        model,
                        selected,
                        candidate,
                        comparisons,
                        elapsed
                    )
                    clear_comparison_failure(case_id, selected)
                    println(
                        "PASS\t", case_id, "\t", selected.id,
                        "\treferences=", length(comparisons),
                        "\telapsed_s=", round(elapsed; digits = 3)
                    )
                catch exception
                    elapsed = time() - started
                    failure_path = record_comparison_failure(
                        case_id,
                        model,
                        selected,
                        exception,
                        elapsed
                    )
                    push!(failures,
                        (;
                            case = string(case_id),
                            candidate_kind = string(selected.kind),
                            candidate_formula = string(selected.identifier),
                            exception_type = string(typeof(exception)),
                            message = sprint(showerror, exception),
                            elapsed_seconds = elapsed,
                            artifact = failure_path
                        ))
                    println(
                        stderr,
                        "FAIL\t", case_id, "\t", selected.id,
                        "\t", sprint(showerror, exception)
                    )
                    continue
                end
            end
            append!(rows, element_rows(case_id, selected, comparisons, path))
            flush(stdout)
            flush(stderr)
        end
    end

    summary_path = table(
        joinpath(COMPARISON_ROOT, "elementwise_rms.tsv"),
        rows,
        (
            :case,
            :reference_backend,
            :reference_formula,
            :reference_field,
            :candidate_kind,
            :candidate_formula,
            :row,
            :column,
            :Z_rms_absolute,
            :Z_rms_relative,
            :Y_rms_absolute,
            :Y_rms_relative,
            :artifact
        )
    )
    skip_path = table(
        joinpath(COMPARISON_ROOT, "skipped.tsv"),
        skipped,
        (:case, :candidate_kind, :candidate_formula, :reason)
    )
    failure_path = table(
        joinpath(COMPARISON_ROOT, "failures.tsv"),
        failures,
        (
            :case,
            :candidate_kind,
            :candidate_formula,
            :exception_type,
            :message,
            :elapsed_seconds,
            :artifact
        )
    )
    jld_path = joinpath(COMPARISON_ROOT, "summary.jld2")
    temporary = jld_path * ".new"
    JLD2.jldsave(
        temporary;
        schema_version = COMPARISON_SCHEMA_VERSION,
        kind = :cross_backend_formulation_summary,
        metric = :elementwise_rms_across_frequency,
        monte_carlo = false,
        engine_sha256 = ENGINE_SHA256,
        comparison_source_sha256 = COMPARISON_SOURCE_SHA256,
        cases = string.(requested),
        registered_earth_impedance = collect(EARTH_IMPEDANCE.formulas()),
        registered_earth_admittance = collect(EARTH_ADMITTANCE.formulas()),
        rows,
        skipped,
        failures,
        recorded_at_utc = string(now(UTC))
    )
    mv(temporary, jld_path; force = true)
    write_hash(jld_path)
    println(
        "COMPLETE\telement_rows=", length(rows),
        "\tskipped=", length(skipped),
        "\tfailures=", length(failures),
        "\tsummary=", summary_path,
        "\tskips=", skip_path,
        "\tfailure_report=", failure_path,
        "\tjld2=", jld_path
    )
    return isempty(failures) ? 0 : 1
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(comparison_main())
