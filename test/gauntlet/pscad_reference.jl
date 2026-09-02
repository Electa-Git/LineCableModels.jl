using Dates
using JLD2
using LineCableModels
using SHA

module PSCADReferenceCaseLoader
using SHA
import TOML
import LineCableModels
using LineCableModels: AbstractGrid, Grid, Gridspace

const GAUNTLET_ROOT = @__DIR__
include(joinpath(GAUNTLET_ROOT, "case_loader.jl"))
include(joinpath(GAUNTLET_ROOT, "reference_grid.jl"))
end

module GauntletSupport
using LineCableModels
const GAUNTLET_ROOT = @__DIR__
const WORK_ROOT = joinpath(GAUNTLET_ROOT, "benchmarks", ".work")
function benchmark_metadata end
function formulation_record end
include(joinpath(GAUNTLET_ROOT, "pscad", "PSCADBenchmarks.jl"))
end

using .PSCADReferenceCaseLoader: case_index, reference_case
using .GauntletSupport: benchmark_metadata, formulation_record
using .GauntletSupport.PSCADBenchmarks

const SCHEMA_VERSION = 1
const ROOT = joinpath(
    pkgdir(LineCableModels),
    ".linecablemodels",
    "gauntlet",
    "references",
    "pscad"
)

const PSCAD_CATALOGUE = (
    (
        id = :DeriSemlyen1981,
        field = :overhead,
        selector = :DeriSemlyen1981,
        stem = "deri"
    ),
    (
        id = :DirectNumericalIntegration,
        field = :overhead,
        selector = PSCADBenchmarks.DirectNumericalIntegration(:overhead),
        stem = "direct_overhead"
    ),
    (
        id = :WedepohlWilcox1973,
        field = :underground,
        selector = :WedepohlWilcox1973,
        stem = "wedepohl"
    ),
    (
        id = :DirectNumericalIntegration,
        field = :underground,
        selector = PSCADBenchmarks.DirectNumericalIntegration(:underground),
        stem = "direct_ground"
    ),
    (
        id = :Saad1996,
        field = :underground,
        selector = :Saad1996,
        stem = "saad"
    ),
    (
        id = :Ametani2009,
        field = :mixed,
        selector = :Ametani2009,
        stem = "ametani"
    ),
    (
        id = :Lucca1994,
        field = :mixed,
        selector = :Lucca1994,
        stem = "lucca"
    )
)

function placement(model)
    heights = getproperty.(model.problem.system.positions, :y)
    all(<(0), heights) && return :underground
    all(>(0), heights) && return :overhead
    return :mixed
end

function variant_id(selected)
    Symbol(lowercase(string(selected.field)), "_", lowercase(string(selected.id)))
end

function formulation(selected)
    Formulation(:pscad; earth_impedance = selected.selector)
end

function artifact_path(case_id, selected)
    joinpath(
        ROOT,
        string(case_id),
        string(variant_id(selected)),
        "reference.jld2"
    )
end

function source_digest()
    roots = (
        joinpath(pkgdir(LineCableModels), "src", "importexport", "pscad"),
        joinpath(@__DIR__, "pscad")
    )
    paths = reduce(vcat,
        [filter(path -> endswith(path, ".jl") || endswith(path, ".ps1"),
             [joinpath(directory, file)
              for (directory, _, files) in walkdir(root)
              for file in files])
         for root in roots])
    append!(paths, (
        @__FILE__,
        joinpath(@__DIR__, "case_loader.jl"),
        joinpath(@__DIR__, "reference_grid.jl")
    ))
    sort!(unique!(paths))
    contents = UInt8[]
    for path in paths
        append!(contents, read(path))
        push!(contents, 0x00)
    end
    return bytes2hex(sha256(contents))
end

const SOURCE_SHA256 = source_digest()

function write_hash(path)
    open(path * ".sha256", "w") do io
        println(io, bytes2hex(sha256(read(path))), "  ", basename(path))
    end
    return path
end

function archive(path)
    digest = bytes2hex(sha256(read(path)))
    stem, extension = splitext(path)
    destination = "$(stem)_stale_$(first(digest, 12))$(extension)"
    if isfile(destination)
        rm(path)
    else
        mv(path, destination)
    end
    rm(path * ".sha256"; force = true)
    write_hash(destination)
    return destination
end

function valid_existing(path, model, selected)
    isfile(path) || return false
    try
        document = JLD2.load(path)
        return document["schema_version"] == SCHEMA_VERSION &&
               document["status"] === :complete &&
               document["case_source_sha256"] == model.source_sha256 &&
               document["source_sha256"] == SOURCE_SHA256 &&
               document["frequencies"] == model.problem.frequencies &&
               document["port_order"] == model.port_order &&
               document["formula"] == selected.id &&
               document["field"] == selected.field
    catch
        return false
    end
end

function diagnostic(path)
    isfile(path) ? read(path, String) : ""
end

function record_result(case_id, model, selected, parameters, execution, elapsed)
    path = artifact_path(case_id, selected)
    mkpath(dirname(path))
    temporary = path * ".new"
    formulation_value = formulation(selected)
    JLD2.jldsave(
        temporary;
        schema_version = SCHEMA_VERSION,
        kind = :pscad_line_parameters_reference,
        status = :complete,
        case_id = string(case_id),
        case_source_sha256 = model.source_sha256,
        source_sha256 = SOURCE_SHA256,
        frequencies = copy(parameters.f),
        port_order = copy(model.port_order),
        basis = :pul,
        domain = :PhaseDomain,
        formula = selected.id,
        field = selected.field,
        formulation = formulation_record(formulation_value),
        Z = copy(parameters.Z.values),
        Y = copy(parameters.Y.values),
        pscad_version = execution.pscad_version,
        pscad_elapsed_seconds = execution.elapsed_seconds,
        elapsed_seconds = elapsed,
        diagnostics = (
            console = diagnostic(execution.console_path),
            stdout = diagnostic(execution.stdout_path),
            stderr = diagnostic(execution.stderr_path)
        ),
        recorded_at_utc = string(now(UTC))
    )
    mv(temporary, path; force = true)
    write_hash(path)
    failure = joinpath(dirname(path), "failure.jld2")
    rm(failure; force = true)
    rm(failure * ".sha256"; force = true)
    return path
end

function record_failure(case_id, model, selected, exception, elapsed)
    path = joinpath(dirname(artifact_path(case_id, selected)), "failure.jld2")
    mkpath(dirname(path))
    temporary = path * ".new"
    JLD2.jldsave(
        temporary;
        schema_version = SCHEMA_VERSION,
        kind = :pscad_line_parameters_failure,
        status = :execution_failure,
        case_id = string(case_id),
        case_source_sha256 = model.source_sha256,
        source_sha256 = SOURCE_SHA256,
        frequencies = copy(model.problem.frequencies),
        port_order = copy(model.port_order),
        formula = selected.id,
        field = selected.field,
        exception_type = string(typeof(exception)),
        message = sprint(showerror, exception),
        elapsed_seconds = elapsed,
        recorded_at_utc = string(now(UTC))
    )
    mv(temporary, path; force = true)
    write_hash(path)
    return path
end

function write_tsv(path, rows, names)
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

function run_case(case_id, config)
    model = reference_case(case_id)
    case_placement = placement(model)
    rows = NamedTuple[]
    failures = NamedTuple[]
    for selected in PSCAD_CATALOGUE
        if selected.field !== case_placement
            push!(rows,
                (;
                    case = string(case_id),
                    placement = string(case_placement),
                    field = string(selected.field),
                    formula = string(selected.id),
                    status = "not_applicable",
                    elapsed_seconds = NaN,
                    artifact = ""
                ))
            println(
                "SKIP\t", case_id, "\t", variant_id(selected),
                "\tcase placement is ", case_placement
            )
            continue
        end

        path = artifact_path(case_id, selected)
        if valid_existing(path, model, selected)
            write_hash(path)
            failure = joinpath(dirname(path), "failure.jld2")
            rm(failure; force = true)
            rm(failure * ".sha256"; force = true)
            document = JLD2.load(path)
            push!(rows,
                (;
                    case = string(case_id),
                    placement = string(case_placement),
                    field = string(selected.field),
                    formula = string(selected.id),
                    status = "reused",
                    elapsed_seconds = document["elapsed_seconds"],
                    artifact = path
                ))
            println("REUSE\t", case_id, "\t", variant_id(selected), "\t", path)
            continue
        elseif isfile(path)
            println("STALE_ARCHIVED\t", case_id, "\t", archive(path))
        end

        selected_formulation = formulation(selected)
        started = time()
        try
            parameters = compute(
                model.problem,
                selected_formulation;
                options = (
                    remote = config,
                    output_stem = selected.stem,
                    output_basis = :pul,
                    verbosity = (default = 0, PSCAD = 2)
                )
            )
            elapsed = time() - started
            size(parameters.Z.values) == model.expected_size || error(
                "$case_id PSCAD Z size $(size(parameters.Z.values)) does not match $(model.expected_size)",
            )
            parameters.f == model.problem.frequencies || error(
                "$case_id PSCAD changed the selected reference grid",
            )
            execution = benchmark_metadata(model.problem, selected_formulation)
            path = record_result(
                case_id,
                model,
                selected,
                parameters,
                execution,
                elapsed
            )
            push!(rows,
                (;
                    case = string(case_id),
                    placement = string(case_placement),
                    field = string(selected.field),
                    formula = string(selected.id),
                    status = "complete",
                    elapsed_seconds = elapsed,
                    artifact = path
                ))
            println(
                "PASS\t", case_id, "\t", variant_id(selected),
                "\telapsed_s=", round(elapsed; digits = 3),
                "\t", path
            )
        catch exception
            elapsed = time() - started
            path = record_failure(case_id, model, selected, exception, elapsed)
            push!(failures,
                (;
                    case = string(case_id),
                    field = string(selected.field),
                    formula = string(selected.id),
                    exception_type = string(typeof(exception)),
                    message = sprint(showerror, exception),
                    elapsed_seconds = elapsed,
                    artifact = path
                ))
            println(
                stderr,
                "FAIL\t", case_id, "\t", variant_id(selected),
                "\t", sprint(showerror, exception),
                "\t", path
            )
        end
        flush(stdout)
        flush(stderr)
    end
    return rows, failures
end

function main(args = ARGS)
    available = sort!(collect(keys(case_index())); by = string)
    requested = isempty(args) ? available : Symbol.(args)
    all(in(available), requested) || error(
        "unknown case requested; available cases are $(join(available, ", "))",
    )

    config = PSCADBenchmarks._load_config()
    rows = NamedTuple[]
    failures = NamedTuple[]
    println(
        "PLAN\tcases=", length(requested),
        "\tpscad_formulations=", length(PSCAD_CATALOGUE),
        "\tpolicy=all_applicable_pscad_earth_fields"
    )
    for case_id in requested
        case_rows, case_failures = run_case(case_id, config)
        append!(rows, case_rows)
        append!(failures, case_failures)
    end

    summary_path = write_tsv(
        joinpath(ROOT, "summary.tsv"),
        rows,
        (:case, :placement, :field, :formula, :status, :elapsed_seconds, :artifact)
    )
    failure_path = write_tsv(
        joinpath(ROOT, "failures.tsv"),
        failures,
        (:case, :field, :formula, :exception_type, :message, :elapsed_seconds, :artifact)
    )
    println(
        "COMPLETE\treferences=", count(
            row -> row.status in ("complete", "reused"), rows
        ),
        "\tskipped=", count(row -> row.status == "not_applicable", rows),
        "\tfailures=", length(failures),
        "\tsummary=", summary_path,
        "\tfailure_report=", failure_path
    )
    return isempty(failures) ? 0 : 1
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main())
