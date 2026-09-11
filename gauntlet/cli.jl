using LineCableModels, JLD2, SHA
import TOML
include(joinpath(@__DIR__, "Gauntlet.jl"))
using .Gauntlet
const GAUNTLET_ROOT=@__DIR__

function usage(io::IO = stdout)
    print(
        io, """
Usage: lcm gauntlet case import|list|show|validate|catalogue [options]
       lcm gauntlet run --definition FILE.jl --directory DIR [--configuration FILE.jl]
                       [--frequencies FILE.toml | --grid log|linear --count N --bounds LOW,HIGH]
                       [--on-error continue|fail] [--resume] [--recover-solvers]
                       [--progress auto|plain|off]
       lcm gauntlet resume --directory DIR [--recover-solvers] [--progress auto|plain|off]
       lcm gauntlet status --directory DIR [--watch]
       lcm gauntlet compare --definition FILE.toml --output DIR
       lcm gauntlet lock --directory DIR --output DIR [--benchmark ID[,ID]]
                        [--expected SNAPSHOT] [--note TEXT] [--illustrations FILE.toml]
       lcm gauntlet package --definition RELEASE.toml --output DIR
       lcm gauntlet bind --package DIR --url URL [--artifacts FILE.toml] [--current]

A definition file returns a BenchmarkDefinition, a vector of definitions, or a
constructor accepting declared keyword arguments. Configuration returns a named
tuple of those arguments (for example reference_options=(remote=station,)).
Frequency overrides are passed to the definition constructor before materialization.
Frequency files contain `frequencies = [ ... ]` in Hz; grids require spacing, count
and both bounds. Omitted overrides preserve the benchmark's declared frequencies.
Case import: --id ID --source FILE [--project DIR] [--description TEXT] [--dry-run] [--force]
Case catalogue: --output FILE [--check]
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
    occursin(Gauntlet._CASE_IDENTIFIER, raw) || throw(ArgumentError(
        "case identifiers must match $(Gauntlet._CASE_IDENTIFIER); got $(repr(raw))",
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
    Gauntlet._validate_loaded_problem(definition, problem)
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
    return """function build_case(::Val{:$(id)}, _)
    LineCableModels.import_data(:json,LineCableModels.Engine.LineParametersProblem;
        file_name=joinpath(@__DIR__,"data","$(id).json"))
end

case_definition(Base.Fix1(build_case,Val(:$(id))),:$(id),(;),$(repr(ports));
    description=$(repr(description)),assets=("data/$(id).json",))
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
    root = abspath(option(args, "--case-root"; default = Gauntlet.CASE_ROOT))
    index_path = joinpath(root, "index.toml")
    wrapper = joinpath(root, "$(id).jl")
    data = joinpath(root, "data", "$(id).json")
    exists = isfile(wrapper) || isfile(data) ||
             (isfile(index_path) &&
              haskey(Gauntlet._case_index_document(index_path), string(id)))
    exists && !flag(args, "--force") &&
        throw(ArgumentError(
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
        entries = isfile(index_path) ?
                  Dict(
            Symbol(key) => String(value)
        for (key, value) in Gauntlet._case_index_document(index_path)
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
        (property = "Frequency range [Hz]",
            value = "$(first(frequencies)) – $(last(frequencies))"),
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
        push!(rows,
            (
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
        push!(records,
            (
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
    Gauntlet.semantic_sha256([(id = record.id, source = record.source_sha256,
                                  input = record.input_sha256)
                              for record in records])
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
    revision = repository_revision()
    JLD2.jldsave(
        temporary;
        schema_version = 1,
        kind = :gauntlet_case_catalogue,
        catalogue_sha256 = digest,
        repository_commit = revision.commit,
        repository_dirty = revision.dirty,
        cases = records
    )
    mv(temporary, output; force = true)
    println("catalogue=", output)
end

function main(args = ARGS)
    isempty(args) && return usage()
    args[1] in ("--help", "-h", "help") && return usage()
    if args[1] in ("lock","package","bind")
        allowed=args[1] == "lock" ? ("--directory","--output","--benchmark","--expected","--note","--illustrations") :
            args[1] == "package" ? ("--definition","--output") : ("--package","--url","--artifacts","--current")
        seen=Set{String}()
        index=2
        while index <= length(args)
            name=args[index]
            name in allowed && !(name in seen) || throw(ArgumentError("unknown or repeated option: $name"))
            push!(seen,name)
            index+=name == "--current" ? 1 : 2
            index <= length(args)+1 || throw(ArgumentError("$name requires a value"))
        end
    end
    if args[1] == "compare"
        foreach(println,
            compare_saved(required_option(args, "--definition");
                directory = required_option(args, "--output")))
        return
    elseif args[1] == "package"
        println(package_collection(required_option(args,"--definition");output=required_option(args,"--output")))
        return
    elseif args[1] == "bind"
        println(bind_published_artifact(required_option(args,"--package"),required_option(args,"--url");
            artifacts_toml=option(args,"--artifacts";default=Gauntlet.ARTIFACTS_TOML),current=flag(args,"--current")))
        return
    elseif args[1] == "lock"
        ids=option(args,"--benchmark")
        figure_file=option(args,"--illustrations")
        illustrations=figure_file === nothing ? () : [(path=normpath(joinpath(dirname(abspath(figure_file)),entry["path"])),
            benchmark=entry["benchmark"],caption=entry["caption"],selection=entry["selection"]) for entry in TOML.parsefile(figure_file)["illustrations"]]
        println(lock_campaign(required_option(args,"--directory"),required_option(args,"--output");
            benchmarks=ids === nothing ? nothing : split(ids,','),expected=option(args,"--expected"),
            note=option(args,"--note";default=""),illustrations))
        return
    elseif args[1] == "status"
        flag(args,"--watch") && return watch_campaign(required_option(args,"--directory"))
        for row in campaign_status(required_option(args, "--directory"))
            println(row.id, '\t', row.state, '\t', row.identity, '\t', row.previous ? "previous result retained" : "", '\t', row.message)
        end
        return
    elseif args[1] == "resume"
        outcomes=resume_campaign(required_option(args, "--directory");
            recover_solvers=flag(args,"--recover-solvers"),
            progress=Symbol(option(args,"--progress";default="auto")))
        all(row -> row.state === :complete, outcomes) ||
            error("campaign contains failed calculations; inspect status")
        return
    elseif args[1] == "run"
        source=abspath(required_option(args, "--definition"))
        configuration=option(args, "--configuration")
        settings=configuration === nothing ? (;) :
                 Base.include(Main, abspath(configuration))
        settings isa NamedTuple ||
            throw(ArgumentError("configuration must return explicit constructor keyword arguments"))
        frequency_file=option(args, "--frequencies")
        grid=option(args, "--grid")
        frequency_file === nothing || grid === nothing ||
            throw(ArgumentError("choose a frequency file or a full grid"))
        (grid === nothing && ("--count" in args || "--bounds" in args)) &&
            throw(ArgumentError("grid requires --grid, --count and --bounds"))
        if frequency_file !== nothing
            settings=merge(settings, (frequencies = TOML.parsefile(frequency_file)["frequencies"],))
        elseif grid !== nothing
            grid in ("log", "linear") ||
                throw(ArgumentError("--grid must be log or linear"))
            count=parse(Int, required_option(args, "--count"))
            bounds=parse.(Float64, split(required_option(args, "--bounds"), ','))
            count > 0 && length(bounds)==2 && 0 < bounds[1] < bounds[2] ||
                throw(ArgumentError("positive increasing bounds and a positive count are required"))
            frequencies=grid == "log" ?
                        10.0 .^ range(log10(bounds[1]), log10(bounds[2]); length = count) :
                        collect(range(bounds...; length = count))
            settings=merge(settings, (; frequencies))
        end
        value=Base.include(Gauntlet, source)
        if value isa Union{BenchmarkDefinition, AbstractVector{<:BenchmarkDefinition}}
            isempty(settings) ||
                throw(ArgumentError("overrides require a declaration constructor"))
        else
            value=Base.invokelatest(value; settings...)
        end
        definitions=value isa BenchmarkDefinition ? [value] : value
        definitions isa AbstractVector{<:BenchmarkDefinition} ||
            throw(ArgumentError("definition must produce benchmark definitions"))
        outcomes=Base.invokelatest(
            run_campaign, required_option(args, "--directory"), definitions;
            on_error = Symbol(option(args, "--on-error"; default = "continue")),
            resume=flag(args,"--resume"),recover_solvers=flag(args,"--recover-solvers"),
            progress=Symbol(option(args,"--progress";default="auto")),
            execution_sources = configuration === nothing ?
                                [(path = source, module_name = :Gauntlet)] :
                                [(path = abspath(configuration), module_name = :Main),
                (path = source, module_name = :Gauntlet)])
        all(row -> row.state === :complete, outcomes) ||
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
