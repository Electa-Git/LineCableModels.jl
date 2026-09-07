using JLD2
using LineCableModels
using SHA
import TOML

const GAUNTLET_ROOT = @__DIR__
include(joinpath(GAUNTLET_ROOT, "case_loader.jl"))
include(joinpath(GAUNTLET_ROOT, "provenance.jl"))

if !isempty(ARGS) && first(ARGS) in ("run", "resume", "status")
    include(joinpath(GAUNTLET_ROOT, "runner.jl"))
end

function usage(io::IO = stdout)
    print(io, """
Usage: lcm gauntlet case <command> [options]
       lcm gauntlet run --directory DIR [--cases ID,ID] [--backends coaxial,fem,pscad]
                       [--formulas catalogue|default] [--dielectric default|Ametani2004]
                       [--select SLOT=ID,ID ...] [--combine product|zip]
                       [--propagation deterministic,linear_error,monte_carlo]
                       [--uncertainty-percent P --uncertainty-tags geometry,cable_layer]
                       [--trials N --seed INTEGER]
       lcm gauntlet resume --directory DIR
       lcm gauntlet status --directory DIR

run uses the case frequency grid, normalised to at least 0.1 Hz, and never runs FEM Monte Carlo.
--select maps any Formulation keyword slot to a Grid; use --formulas default with explicit axes.
UQ needs explicit uncertainty settings; Monte Carlo also needs an explicit seed.
Linear propagation requires coaxial; Monte Carlo on FEM is disabled.
resume preserves completed calculations; changed inputs require a new campaign directory.
status checks the live execution lock and distinguishes running from interrupted jobs.
Commands:
  import       Import a trusted Julia file that returns LineParametersProblem
  list         List indexed cases
  show ID      Display one materialised case
  validate     Materialise and validate one case or the complete catalogue
  catalogue    Write or check the detached documentation manifest

Import options:
  --id ID --source FILE [--project DIR] [--description TEXT] [--dry-run] [--force]

Catalogue options:
  --output FILE [--check]
""")
end

function option(args, name; default = nothing)
    index = findfirst(==(name), args)
    index === nothing && return default
    index == length(args) && throw(ArgumentError("$name requires a value"))
    return args[index + 1]
end


function required_option(args, name)
    value = option(args, name)
    value === nothing && throw(ArgumentError("$name is required"))
    return value
end

flag(args, name) = name in args

function checked_id(raw::AbstractString)
    occursin(_CASE_IDENTIFIER, raw) || throw(ArgumentError(
        "case identifiers must match $(_CASE_IDENTIFIER); got $(repr(raw))",
    ))
    return Symbol(raw)
end

function port_order(problem)
    assignments = problem.system.connection_order
    terminals = problem.system.terminal_order
    length(assignments) == length(terminals) || throw(DimensionMismatch(
        "system terminal and connection orders differ",
    ))
    any(iszero, assignments) && throw(ArgumentError(
        "Gauntlet cases require every declared terminal to have a positive phase",
    ))
    sort(assignments) == collect(1:length(assignments)) || throw(ArgumentError(
        "Gauntlet phase assignments must be unique and contiguous from one",
    ))
    ports = Vector{String}(undef, length(assignments))
    for (terminal, phase) in zip(terminals, assignments)
        ports[phase] = "cable:$(terminal.cable):$(terminal.terminal)"
    end
    return ports
end

function validate_problem(problem, id::Symbol, ports)
    problem.system.system_id == string(id) || throw(ArgumentError(
        "imported system ID differs from requested case ID",
    ))
    definition = case_definition(
        _ -> problem,
        id,
        (;),
        ports
    )
    _validate_loaded_problem(definition, problem)
    T = eltype(problem)
    for design in problem.system.designs
        LineCableModels.Engine.flatten(LineCableModels.LineCableModelsCoaxial(), design, T)
    end
    return problem
end

function write_index(index_path, entries)
    temporary = index_path * ".new"
    open(temporary, "w") do io
        println(io, "[cases]")
        for (id, path) in sort!(collect(entries); by = first)
            println(io, string(id), " = ", repr(path))
        end
    end
    mv(temporary, index_path; force = true)
    return index_path
end

function wrapper_text(id, description, ports)
    return """case_definition(
    :$(id),
    (;),
    $(repr(ports));
    description = $(repr(description)),
    assets = (\"data/$(id).json\",)
) do _
    LineCableModels.import_data(
        :json,
        LineCableModels.Engine.LineParametersProblem;
        file_name = joinpath(@__DIR__, \"data\", \"$(id).json\")
    )
end
"""
end

function import_case(args)
    raw_id = required_option(args, "--id")
    source = required_option(args, "--source")
    id = checked_id(raw_id)
    source = abspath(source)
    isfile(source) || throw(ArgumentError("case source does not exist: $source"))
    description = option(
        args,
        "--description";
        default = string(id)
    )
    source_project = abspath(option(args, "--project"; default = GAUNTLET_ROOT))
    isdir(source_project) || throw(ArgumentError(
        "source project directory does not exist: $source_project",
    ))
    isfile(joinpath(source_project, "Project.toml")) || throw(ArgumentError(
        "source project has no Project.toml: $source_project",
    ))
    root = abspath(option(args, "--case-root"; default = CASE_ROOT))
    index_path = joinpath(root, "index.toml")
    wrapper = joinpath(root, "$(id).jl")
    data = joinpath(root, "data", "$(id).json")
    exists = isfile(wrapper) || isfile(data) ||
             (isfile(index_path) && haskey(_case_index_document(index_path), string(id)))
    exists && !flag(args, "--force") && throw(ArgumentError(
        "case :$id already exists; pass --force to replace it",
    ))

    mktempdir() do staging
        staged_data = joinpath(staging, "$(id).json")
        worker = joinpath(GAUNTLET_ROOT, "import_worker.jl")
        command = `$(Base.julia_cmd()) --project=$source_project --startup-file=no $worker $source $staged_data $(string(id))`
        run(command)
        problem = LineCableModels.import_data(
            :json,
            LineCableModels.Engine.LineParametersProblem;
            file_name = staged_data
        )
        ports = port_order(problem)
        validate_problem(problem, id, ports)
        wrapper_contents = wrapper_text(id, description, ports)
        println("case=:$(id)")
        println("ports=", join(ports, ','))
        println("input_sha256=", numerical_input_sha256(problem))
        flag(args, "--dry-run") && return

        mkpath(joinpath(root, "data"))
        staged_wrapper = joinpath(staging, "$(id).jl")
        write(staged_wrapper, wrapper_contents)
        entries = isfile(index_path) ? Dict(
            Symbol(key) => String(value)
            for (key, value) in _case_index_document(index_path)
        ) : Dict{Symbol, String}()
        entries[id] = "$(id).jl"
        prior = Dict(
            path => (isfile(path) ? read(path) : nothing)
            for path in (data, wrapper, index_path)
        )
        try
            cp(staged_data, data; force = true)
            cp(staged_wrapper, wrapper; force = true)
            write_index(index_path, entries)
        catch
            for (path, contents) in prior
                if contents === nothing
                    rm(path; force = true)
                else
                    write(path, contents)
                end
            end
            rethrow()
        end
        println("written=", wrapper)
    end
    return nothing
end

function summary_rows(model)
    problem = model.nominal_problem
    system = problem.system
    earth = problem.earth_props
    frequencies = problem.frequencies
    return NamedTuple[
        (property = "Cables", value = string(length(system.designs))),
        (property = "Terminals", value = string(length(system.terminal_order))),
        (property = "Positive phases", value = string(LineCableModels.nphases(system))),
        (property = "Line length [m]", value = string(system.line_length)),
        (property = "Temperature [°C]", value = string(problem.temperature)),
        (property = "Frequencies", value = string(length(frequencies))),
        (property = "Frequency range [Hz]", value = "$(first(frequencies)) – $(last(frequencies))"),
        (property = "Earth layers including air", value = string(length(earth.layers))),
        (property = "Vertical earth layers", value = string(earth.vertical_layers))
    ]
end

function design_rows(model)
    rows = NamedTuple[]
    designs = model.nominal_problem.system.designs
    retained = Int[]
    for (index, design) in pairs(designs)
        any(other -> designs[other] == design, retained) && continue
        push!(retained, index)
        push!(rows, (
            cable = index,
            cable_id = design.cable_id,
            instances = count(==(design), designs),
            terminals = length(design.terminal_order),
            regions = length(design.geometry.regions),
            outer_diameter_m = 2 * LineCableModels.outer_radius(design)
        ))
    end
    return rows
end

function catalogue_records()
    records = NamedTuple[]
    for id in sort!(collect(keys(case_index())); by = string)
        model = load_case(id)
        push!(records, (
            id = string(id),
            description = model.definition.description,
            problem = LineCableModels.ImportExport.serialize_value(model.nominal_problem),
            port_order = model.port_order,
            parameters = parameter_manifest(model),
            summary = summary_rows(model),
            designs = design_rows(model),
            source_sha256 = model.source_sha256,
            input_sha256 = numerical_input_sha256(model.nominal_problem)
        ))
    end
    return records
end

function catalogue_digest(records)
    semantic_sha256([
        (id = record.id, source = record.source_sha256, input = record.input_sha256)
        for record in records
    ])
end

function catalogue(args)
    output = abspath(required_option(args, "--output"))
    records = catalogue_records()
    digest = catalogue_digest(records)
    if flag(args, "--check")
        isfile(output) || throw(ArgumentError("catalogue manifest is missing: $output"))
        document = JLD2.load(output)
        document["schema_version"] == 1 || throw(ArgumentError(
            "unsupported case-catalogue schema in $output",
        ))
        document["catalogue_sha256"] == digest || throw(ArgumentError(
            "case catalogue is stale; regenerate it with `lcm gauntlet case catalogue`",
        ))
        println("catalogue is current: $output")
        return
    end
    mkpath(dirname(output))
    temporary = output * ".new"
    provenance = repository_provenance()
    JLD2.jldsave(
        temporary;
        schema_version = 1,
        kind = :gauntlet_case_catalogue,
        catalogue_sha256 = digest,
        repository_commit = provenance.commit,
        repository_dirty = provenance.dirty,
        cases = records
    )
    mv(temporary, output; force = true)
    println("catalogue=", output)
end

function main(args = ARGS)
    isempty(args) && return usage()
    args[1] in ("--help", "-h", "help") && return usage()
    if args[1] in ("run", "resume", "status")
        directory = required_option(args, "--directory")
        if args[1] == "status"
            println("job\tstate\tcompleted\trequested\tskipped\tmessage")
            for row in GauntletSupport.campaign_status(directory)
                println(join((row.id, row.state, row.completed, row.requested,
                    row.skipped, replace(row.message, '\n'=>' ', '\t'=>' ')), '\t'))
            end
            return
        end
        if args[1] == "resume"
            GauntletSupport.resume_campaign(directory) || error("campaign contains failed calculations; inspect status")
            return
        end
        selected = option(args, "--cases")
        ids = selected === nothing ? sort!(collect(keys(case_index())); by=string) :
            checked_id.(split(selected, ','))
        backends = Symbol.(split(option(args, "--backends"; default="coaxial,fem,pscad"), ','))
        formulas = option(args, "--formulas"; default="catalogue")
        formulas in ("catalogue", "default") || throw(ArgumentError("--formulas must be catalogue or default"))
        dielectric = Symbol(option(args, "--dielectric"; default="default"))
        choices = GauntletSupport.parse_selections(args)
        combine = Symbol(option(args, "--combine"; default="product"))
        propagation = Symbol.(split(option(args, "--propagation"; default="deterministic"), ','))
        "--trials" in args && :monte_carlo ∉ propagation && throw(ArgumentError(
            "--trials only applies to Monte Carlo propagation"))
        percent = option(args, "--uncertainty-percent")
        tags = option(args, "--uncertainty-tags")
        (percent === nothing) == (tags === nothing) || throw(ArgumentError(
            "supply --uncertainty-percent and --uncertainty-tags together"))
        uncertainty = percent === nothing ? nothing : GauntletSupport.RelativeStandardUncertainty(
            parse(Float64, percent); tags=Symbol.(split(tags, ',')))
        trials = parse(Int, option(args, "--trials"; default=string(GauntletSupport.UQ_MONTE_CARLO_TRIALS)))
        raw_seed = option(args, "--seed")
        seed = raw_seed === nothing ? nothing : parse(UInt64, raw_seed)
        GauntletSupport.run_campaign(directory, ids;
            backends, catalogue=formulas == "catalogue", dielectric, choices, combine,
            propagation, uncertainty, trials, seed) ||
            error("campaign contains failed calculations; inspect status")
        return
    end
    args[1] == "case" || throw(ArgumentError(
        "unknown gauntlet command $(repr(args[1])); expected case, run, resume or status",
    ))
    length(args) >= 2 || return usage()
    command = args[2]
    rest = args[3:end]
    command == "import" && return import_case(rest)
    if command == "list"
        for id in sort!(collect(keys(case_index())); by = string)
            model = load_case(id)
            println(id, '\t', model.definition.description)
        end
        return
    end
    if command == "show"
        length(rest) == 1 || throw(ArgumentError("case show requires one ID"))
        model = load_case(checked_id(only(rest)))
        println(model.definition.description)
        show(stdout, MIME("text/plain"), model.nominal_problem)
        println()
        return
    end
    if command == "validate"
        ids = isempty(rest) ? sort!(collect(keys(case_index())); by = string) :
              [checked_id(only(rest))]
        for id in ids
            model = load_case(id)
            validate_problem(model.nominal_problem, id, model.port_order)
            println("valid\t", id)
        end
        return
    end
    command == "catalogue" && return catalogue(rest)
    throw(ArgumentError("unknown gauntlet case command $(repr(command))"))
end

try
    main()
catch exception
    showerror(stderr, exception)
    println(stderr)
    exit(2)
end
