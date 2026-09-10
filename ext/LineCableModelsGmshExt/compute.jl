struct FEMTeeLogger{A <: AbstractLogger, B <: AbstractLogger} <: AbstractLogger
    first::A
    second::B
end

function Logging.min_enabled_level(logger::FEMTeeLogger)
    min(
        Logging.min_enabled_level(logger.first),
        Logging.min_enabled_level(logger.second)
    )
end
function Logging.catch_exceptions(logger::FEMTeeLogger)
    Logging.catch_exceptions(logger.first) || Logging.catch_exceptions(logger.second)
end
function Logging.shouldlog(logger::FEMTeeLogger, level, module_, group, id)
    return Logging.shouldlog(logger.first, level, module_, group, id) ||
           Logging.shouldlog(logger.second, level, module_, group, id)
end
function Logging.handle_message(
        logger::FEMTeeLogger,
        level,
        message,
        module_,
        group,
        id,
        file,
        line;
        kwargs...
)
    if Logging.shouldlog(logger.first, level, module_, group, id)
        Logging.handle_message(
            logger.first, level, message, module_, group, id, file, line; kwargs...
        )
    end
    if Logging.shouldlog(logger.second, level, module_, group, id)
        Logging.handle_message(
            logger.second, level, message, module_, group, id, file, line; kwargs...
        )
    end
    return nothing
end

struct FEMGmshSession
    owned::Bool
    previous_model::String
    initial_models::Set{String}
    initial_views::Set{Int}
    terminal_option::Float64
    verbosity_option::Float64
    onelab::FEMOnelabSnapshot
end

const FEM_SESSION_LOCK = ReentrantLock()

# Keep this file in place: unlinking a locked inode permits another coordinator
# to lock a replacement file. OS ownership is released when the stream closes
# or its Julia process exits, including after a hard crash.
function _claim_run(run::FEMRun)
    io = open(joinpath(run.path, "coordinator.lock"), "a+")
    seekstart(io)
    status = if Sys.iswindows()
        ccall(:_locking, Cint, (Cint, Cint, Clong), fd(io), 2, 1)
    else
        ccall(:flock, Cint, (Cint, Cint), fd(io), 6)
    end
    if status != 0
        close(io)
        _fem_error(:execution, "FEM run", :ownership,
            "another coordinator owns this run, or the filesystem cannot lock it";
            run_directory=run.path)
    end
    return io
end

function _release_run(io::IOStream)
    if Sys.iswindows()
        seekstart(io)
        ccall(:_locking, Cint, (Cint, Cint, Clong), fd(io), 0, 1)
    end
    close(io)
    return nothing
end

function _start_gmsh(verbosity::Int)
    owned = !Bool(gmsh.is_initialized())
    # Backend-owned meshes must not depend on an unrecorded user gmshrc.
    owned && gmsh.initialize(String[], false, false)
    previous_model = try
        gmsh.model.get_current()
    catch
        ""
    end
    initial_models = Set(String.(gmsh.model.list()))
    initial_views = Set(Int.(gmsh.view.get_tags()))
    terminal_option = gmsh.option.get_number("General.Terminal")
    verbosity_option = gmsh.option.get_number("General.Verbosity")
    onelab = _snapshot_onelab()
    gmsh.option.set_number("General.Terminal", verbosity > 0 ? 1 : 0)
    gmsh.option.set_number("General.Verbosity", verbosity)
    return FEMGmshSession(
        owned,
        previous_model,
        initial_models,
        initial_views,
        terminal_option,
        verbosity_option,
        onelab
    )
end

function _finish_gmsh(session::FEMGmshSession)
    if session.owned
        Gmsh.finalize()
        return nothing
    end
    for tag in setdiff(Set(Int.(gmsh.view.get_tags())), session.initial_views)
        try
            gmsh.view.remove(tag)
        catch
        end
    end
    for name in setdiff(Set(String.(gmsh.model.list())), session.initial_models)
        try
            gmsh.model.set_current(name)
            gmsh.model.remove()
        catch
        end
    end
    isempty(session.previous_model) || try
        gmsh.model.set_current(session.previous_model)
    catch
    end
    gmsh.option.set_number("General.Terminal", session.terminal_option)
    gmsh.option.set_number("General.Verbosity", session.verbosity_option)
    _restore_onelab(session.onelab)
    return nothing
end

function _runtime_root()
    return joinpath(pkgdir(LineCableModels), ".linecablemodels", "fem")
end

function _create_run(runtime_root::String)
    runs = joinpath(runtime_root, "runs")
    mkpath(runs)
    path = mktempdir(runs; prefix = "run-", cleanup = false)
    for directory in ("input", "mesh", "raw", "maps", "logs")
        mkpath(joinpath(path, directory))
    end
    run = FEMRun(path, created, "run created", :none, "")
    _transition!(run, created, "run created")
    return run
end

function _resume_value_matches(existing::AbstractDict, expected::AbstractDict)
    length(existing) == length(expected) || return false
    return all(keys(expected)) do key
        haskey(existing, key) &&
            _resume_value_matches(existing[key], expected[key])
    end
end

function _resume_value_matches(existing::AbstractVector, expected::AbstractVector)
    length(existing) == length(expected) || return false
    return all(splat(_resume_value_matches), zip(existing, expected))
end

_resume_value_matches(existing, expected) = isequal(existing, expected)

function _resume_inputs_match(path::String, model::FEMResolvedModel, inputs::NamedTuple)
    snapshot = joinpath(path, "input", "problem.json")
    computation = joinpath(path, "input", "computation.json")
    state = joinpath(path, "run.json")
    all(isfile, (snapshot, computation, state)) || return false
    existing, recorded, run_state = try
        (JSON3.read(read(snapshot, String)), JSON3.read(read(computation, String)),
            JSON3.read(read(state, String)))
    catch
        return false
    end
    # External Gmsh sessions can carry arbitrary caller-owned meshing settings.
    inputs.owned_gmsh || return false
    if String(run_state.state) == string(completed)
        inputs.getdp_identity === nothing && return false
        inputs.execution.ui && return false
        inputs.execution.mesh_policy === :remesh && return false
        isfile(joinpath(path, "raw", "checksums.json")) || return false
    end
    expected = ImportExport.serialize_value(model.problem)
    comparable = Dict(String(key)=>value for (key,value) in pairs(recorded))
    requested = Dict(String(key)=>value for (key,value) in pairs(JSON3.read(JSON3.write(inputs))))
    # A solver-input schema change is an intentional restart boundary.
    get(comparable, "schema_version", 0) == inputs.schema_version == 6 || return false
    get(comparable, "solver_protocol", 0) == inputs.solver_protocol == 3 || return false
    # Scheduling and executable location do not change the numerical problem.
    for record in (comparable, requested)
        execution = Dict(String(key)=>value for (key,value) in pairs(record["execution"]))
        pop!(execution, "frequency_workers", nothing)
        pop!(execution, "getdp_executable", nothing)
        record["execution"] = execution
        pop!(record, "getdp_provenance", nothing)
    end
    return _resume_value_matches(existing, expected) &&
        _resume_value_matches(comparable, requested)
end

function _resume_run(
        runtime_root::String,
        requested::Union{Nothing, Symbol, String},
        model::FEMResolvedModel,
        inputs::NamedTuple
)
    requested === nothing && return _create_run(runtime_root)
    runs = joinpath(runtime_root, "runs")
    mkpath(runs)
    runs_path = realpath(runs)
    candidate = if requested === :latest
        directories = filter(isdir, readdir(runs_path; join = true))
        sort!(directories; by = path -> stat(path).mtime, rev = true)
        index = findfirst(
            path -> _resume_inputs_match(path, model, inputs), directories
        )
        index === nothing ? nothing : directories[index]
    elseif requested isa String
        path = abspath(requested)
        isdir(path) || throw(ArgumentError(
            "resume_run_directory does not exist: $path",
        ))
        _resume_inputs_match(path, model, inputs) || throw(ArgumentError(
            "resume_run_directory needs matching problem, constitutive inputs, " *
            "solver and adapter identities, and backend-owned meshing settings: $path; " *
            "start a new run when these inputs changed or the run predates input snapshots",
        ))
        path
    else
        throw(ArgumentError(
            "resume_run_directory must be nothing, :latest, or a path string",
        ))
    end
    candidate === nothing && return _create_run(runtime_root)
    path = realpath(candidate)
    dirname(path) == runs_path || throw(ArgumentError(
        "resume_run_directory must be an immediate child of $runs_path",
    ))
    for directory in ("input", "mesh", "raw", "maps", "logs")
        isdir(joinpath(path, directory)) || throw(ArgumentError(
            "resume_run_directory is missing $directory/: $path",
        ))
    end
    document = try
        JSON3.read(read(joinpath(path, "run.json"), String))
    catch
        nothing
    end
    mesh_source = document === nothing ? :none :
                  Symbol(String(document.mesh_source))
    mesh_fingerprint = document === nothing ? "" :
                       String(document.mesh_fingerprint)
    run = FEMRun(
        path,
        created,
        "resuming interrupted run",
        mesh_source,
        mesh_fingerprint,
        document === nothing ? 0 : Int(document.getdp_invocations)
    )
    if document !== nothing
        run.completed_columns = Int(get(document, :completed_columns, 0))
        run.completed_frequencies = Int(get(document, :completed_frequencies, 0))
        if String(document.state) == string(completed)
            run.state = completed
            run.message = "reading completed run"
            return run
        end
    end
    @info "Resuming compatible FEM run" run_directory=path
    return run
end

function _transition!(run::FEMRun, state::FEMRunState, message::AbstractString)
    run.state = state
    run.message = String(message)
    _write_json_atomic(joinpath(run.path, "run.json"),
        (
            schema = "LineCableModels.FEMRun",
            version = 2,
            state = string(state),
            message = run.message,
            mesh_source = String(run.mesh_source),
            mesh_fingerprint = run.mesh_fingerprint,
            getdp_invocations = run.getdp_invocations,
            completed_columns = run.completed_columns,
            completed_frequencies = run.completed_frequencies,
            updated_unix_seconds = time()
        ))
    return state
end

function _prepare_run_inputs!(run::FEMRun, model::FEMResolvedModel)
    asset_directory = joinpath(run.path, "input", "getdp")
    mkpath(asset_directory)
    for (name, path) in pairs(_getdp_assets(asset_directory))
        contents = getproperty(FEM_GETDP_SOURCES, name)
        if isfile(path)
            read(path, String) == contents || _fem_error(
                :getdp, "GetDP", :assets,
                "stored solver sources differ from this implementation; start a new run";
                run_directory = run.path)
        else
            write(path, contents)
        end
    end
    problem_path = joinpath(run.path, "input", "problem.json")
    model_data_path = joinpath(run.path, "input", "model_data.pro")
    _write_problem_snapshot(problem_path, model.problem)
    _write_model_data(model_data_path, model)
    mkpath(joinpath(run.path, "raw", "jobs"))
    open(joinpath(run.path, "raw", "Z.tsv"), "w") do io
        println(io, join(FEM_RAW_HEADER, '\t'))
    end
    open(joinpath(run.path, "raw", "P.tsv"), "w") do io
        println(io, join(FEM_RAW_HEADER, '\t'))
    end
    open(joinpath(run.path, "raw", "scan_complete.tsv"), "w") do io
        println(io, join(FEM_COMPLETE_HEADER, '\t'))
    end
    return model_data_path
end

function _fem_computation_options(options::NamedTuple)
    allowed = (
        :verbosity,
        :output_basis,
        :trace,
        :on_result,
        :log_file,
        :resume_run_directory
    )
    unknown = filter(key -> key ∉ allowed, keys(options))
    isempty(unknown) || throw(ArgumentError(
        "unknown LineCableModelsFEM computation options: $(sort!(collect(unknown)))",
    ))
    log_file = get(options, :log_file, nothing)
    log_file isa Union{Nothing, AbstractString} || throw(ArgumentError(
        "log_file must be a path string or nothing",
    ))
    log_file === "" && throw(ArgumentError("log_file cannot be empty"))
    resume = get(options, :resume_run_directory, nothing)
    resume isa Union{Nothing, Symbol, AbstractString} || throw(ArgumentError(
        "resume_run_directory must be nothing, :latest, or a path string",
    ))
    resume isa Symbol && resume !== :latest &&
        throw(ArgumentError(
            "the only symbolic resume_run_directory is :latest",
        ))
    resume === "" && throw(ArgumentError(
        "resume_run_directory cannot be empty",
    ))
    standard_keys = filter(
        key -> key ∉ (:log_file, :resume_run_directory), keys(options)
    )
    standard_values = map(key -> getproperty(options, key), standard_keys)
    standard = NamedTuple{Tuple(standard_keys)}(Tuple(standard_values))
    execution = computation_options(
        LineCableModels.LineCableModelsCoaxial, standard
    )
    return (;
        execution...,
        log_file = log_file === nothing ? nothing : String(log_file),
        resume_run_directory = resume isa AbstractString ? String(resume) : resume
    )
end

function _headless_solve!(
        run::FEMRun,
        model::FEMResolvedModel,
        formulation::LineCableModelsFEM,
        execution::NamedTuple,
        runtime_root::String,
        inputs::NamedTuple
)
    geometry = _build_geometry!(
        model, "LineCableModelsFEM-$(basename(run.path))"
    )
    _transition!(run, geometry_ready, "geometry ready")
    @info "FEM geometry ready" run_directory=run.path
    mesh_paths = _select_meshes!(
        run, model, geometry, formulation, runtime_root
    )
    _transition!(run, mesh_ready, "mesh ready")
    @info "FEM mesh ready" source=run.mesh_source fingerprint=run.mesh_fingerprint
    model_data_path = _prepare_run_inputs!(run, model)
    _publish_transport!(
        run, model, model_data_path, last(mesh_paths), formulation
    )
    _transition!(run, running, "GetDP frequency batches running")
    @info "Starting isolated GetDP frequency batches" workers=formulation.execution.frequency_workers
    _run_getdp!(run, model, formulation, mesh_paths)
    scan = _parse_scan(run, model, formulation)
    _write_scan_checksums(run, scan)
    gmsh.onelab.set_number(_onelab_name("completion_status"), [1.0])
    _transition!(run, completed, "results validated")
    parameters = _line_parameters(run, model, formulation, execution, scan, inputs)
    @info "FEM scan completed successfully"
    return parameters
end

function _ui_solve!(
        run::FEMRun,
        model::FEMResolvedModel,
        formulation::LineCableModelsFEM,
        execution::NamedTuple,
        runtime_root::String,
        inputs::NamedTuple
)
    geometry = _build_geometry!(
        model, "LineCableModelsFEM-$(basename(run.path))"
    )
    _transition!(run, geometry_ready, "geometry ready")
    state = :geometry_ready
    _publish_ui!(run, model, state)
    Bool(gmsh.fltk.is_available()) || gmsh.fltk.initialize()
    gmsh.fltk.wait(0.05)
    yield()
    mesh_paths = nothing
    model_data_path = nothing
    parameters = nothing
    while true
        available = Bool(gmsh.fltk.is_available())
        transition = _ui_transition(state, _take_ui_action(), available)
        if transition === :closed_before_mesh
            _transition!(run, not_executed, "UI closed before mesh generation")
            _fem_error(
                :not_executed,
                model.problem.system.system_id,
                :ui,
                "Gmsh UI closed before mesh generation; no FEM solve was executed";
                run_directory = run.path
            )
        elseif transition === :closed_before_solve
            _transition!(run, not_executed, "UI closed after mesh generation before solve")
            _fem_error(
                :not_executed,
                model.problem.system.system_id,
                :ui,
                "Gmsh UI closed after mesh generation but before Run model; " *
                "no FEM solve was executed";
                run_directory = run.path
            )
        elseif !available
            parameters === nothing || return parameters
            _transition!(run, not_executed, "UI closed without a solve")
            _fem_error(
                :not_executed,
                model.problem.system.system_id,
                :ui,
                "Gmsh UI closed before a FEM solve was executed";
                run_directory = run.path
            )
        elseif transition === :mesh_required
            _set_ui_status(:mesh_required)
        elseif transition === :mesh_requested
            mesh_paths = _select_meshes!(
                run, model, geometry, formulation, runtime_root
            )
            model_data_path = _prepare_run_inputs!(run, model)
            _publish_transport!(
                run, model, model_data_path, last(mesh_paths), formulation
            )
            _transition!(run, mesh_ready, "mesh ready")
            state = _set_ui_status(:mesh_ready)
        elseif transition === :solve_requested
            mesh_paths === nothing && begin
                _set_ui_status(:mesh_required)
                gmsh.fltk.wait(0.05)
                continue
            end
            state = _set_ui_status(:running)
            _transition!(run, running, "GetDP frequency batches running")
            progress = (-1, -1)
            _run_getdp!(run, model, formulation, mesh_paths; pump=() -> begin
                # fltk.wait can initialize a GUI again after it was closed.
                Bool(gmsh.fltk.is_available()) || return false
                current = (run.completed_frequencies, run.completed_columns)
                if current != progress
                    gmsh.onelab.set_number(_onelab_name("ui/completed_frequencies"), [current[1]])
                    gmsh.onelab.set_number(_onelab_name("ui/completed_columns"), [current[2]])
                    progress = current
                end
                gmsh.fltk.wait(0.01)
                Bool(gmsh.fltk.is_available())
            end)
            scan = _parse_scan(run, model, formulation)
            _write_scan_checksums(run, scan)
            formulation.execution.plot_field_maps && _merge_maps!(scan.map_paths)
            gmsh.onelab.set_number(_onelab_name("completion_status"), [1.0])
            gmsh.onelab.set_number(_onelab_name("ui/completed_frequencies"), [run.completed_frequencies])
            gmsh.onelab.set_number(_onelab_name("ui/completed_columns"), [run.completed_columns])
            _transition!(run, completed, "results ready")
            parameters = _line_parameters(run, model, formulation, execution, scan, inputs)
            state = _set_ui_status(:results_ready)
        end
        gmsh.fltk.wait(0.05)
        yield()
    end
end

function _attach_run_directory(
        exception::LineCableModelsFEMError,
        run::FEMRun
)
    exception.run_directory !== nothing && return exception
    return LineCableModelsFEMError(
        exception.category,
        exception.object_id,
        exception.field,
        exception.message;
        run_directory = run.path
    )
end

function _compute_fem(
        problem::LineParametersProblem{Float64},
        formulation::LineCableModelsFEM,
        execution::NamedTuple,
        model::FEMResolvedModel = _resolved_fem_model(problem, formulation)
)
    runtime_root = _runtime_root()
    inputs = _fem_input_record(model, formulation)
    run = _resume_run(
        runtime_root, execution.resume_run_directory, model, inputs
    )
    if run.state === completed
        # Read-only reuse: no Gmsh session, scratch reset, log append, state
        # transition, or successful-run cleanup may touch historical evidence.
        scan = _parse_scan(run, model, formulation)
        _check_scan_checksums(run, scan)
        @info "FEM reuses completed resolved inputs" run_directory=run.path
        return _line_parameters(run, model, formulation, execution, scan, inputs)
    end
    ownership = _claim_run(run)
    parameters = try
        # Another coordinator can finish between candidate selection and our
        # lock acquisition. Refresh counters/state while holding ownership,
        # and preserve a now-completed run as a read-only result.
        if isfile(joinpath(run.path, "input", "computation.json")) &&
                isfile(joinpath(run.path, "input", "problem.json"))
            run = _resume_run(runtime_root, run.path, model, inputs)
            if run.state === completed
                scan = _parse_scan(run, model, formulation)
                _check_scan_checksums(run, scan)
                return _line_parameters(run, model, formulation, execution, scan, inputs)
            end
        end
        _assert_no_live_attempts(run)
        _compute_owned_fem(problem, formulation, execution, model, run, runtime_root, inputs)
    finally
        _release_run(ownership)
    end
    if !formulation.execution.keep_run_directory
        expected_parent = realpath(joinpath(runtime_root, "runs"))
        realpath(dirname(run.path)) == expected_parent || error(
            "refusing to remove FEM run outside the runtime root")
        rm(run.path; recursive=true, force=true)
    end
    return parameters
end

function _compute_owned_fem(problem, formulation, execution, model, run, runtime_root, inputs)
    _write_json_atomic(joinpath(run.path, "input", "computation.json"), inputs)
    session = nothing
    try
        session = _start_gmsh(formulation.execution.gmsh_verbosity)
        parameters = formulation.execution.ui ?
                     _ui_solve!(run, model, formulation, execution, runtime_root, inputs) :
                     _headless_solve!(run, model, formulation, execution, runtime_root, inputs)
        return parameters
    catch exception
        if run.state ∉ (not_executed, cancelled)
            _transition!(run, failed, sprint(showerror, exception))
        end
        if exception isa LineCableModelsFEMError
            throw(_attach_run_directory(exception, run))
        end
        _fem_error(
            :execution,
            problem.system.system_id,
            :backend,
            sprint(showerror, exception);
            run_directory = run.path
        )
    finally
        session === nothing || _finish_gmsh(session)
    end
end

const FEM_ADAPTER_SOURCES = let files = ("model.jl", "geometry.jl", "mesh.jl",
        "onelab.jl", "getdp.jl", "workers.jl", "results.jl", "compute.jl")
    digests = map(files) do file
        path = joinpath(@__DIR__, file)
        Base.include_dependency(path)
        bytes2hex(sha256(read(path)))
    end
    NamedTuple{Symbol.(files)}(digests)
end

function _fem_input_record(model::FEMResolvedModel, formulation::LineCableModelsFEM)
    selection = _getdp_selection(formulation)
    getdp_identity = _getdp_identity(selection.path)
    getdp_provenance = selection
    mesh_path = formulation.execution.mesh_path
    return (
        schema_version = 6,
        solver_protocol = 3,
        mesh_fingerprint = _mesh_fingerprint(model, gmsh.GMSH_API_VERSION),
        materials = [(kind=material.kind, tag=material.physical_tag,
            mu_r=material.mu_r, sigma=real.(material.admittivity),
            omega_epsilon=imag.(material.admittivity)) for material in model.material_plans],
        earth_materials = [(rho=state.rho, eps_r=state.eps_r, mu_r=state.mu_r)
            for state in model.earth_materials],
        air = (eps_r=model.problem.earth_props.layers[1].eps_r,
            mu_r=model.problem.earth_props.layers[1].mu_r),
        mesh_plans = model.mesh_plans,
        region_mesh_sizes = getproperty.(model.region_plans, :mesh_size),
        cable_outer_mesh_sizes = model.cable_outer_mesh_sizes,
        mesh_growth_factor = model.mesh_growth_factor,
        options = formulation.options,
        execution = formulation.execution,
        supplied_mesh = mesh_path === nothing || !isfile(mesh_path) ? nothing :
            bytes2hex(open(sha256, mesh_path)),
        owned_gmsh = !Bool(gmsh.is_initialized()),
        getdp_identity,
        getdp_provenance,
        gmsh_version = gmsh.GMSH_API_VERSION,
        gmsh_library = String(gmsh.lib),
        julia_version = string(VERSION),
        adapter_sources = FEM_ADAPTER_SOURCES,
        solver_sources = map(source -> bytes2hex(sha256(source)), FEM_GETDP_SOURCES)
    )
end

function compute(
        problem::LineParametersProblem,
        formulation::Union{LineCableModelsFEM, AbstractVector{<:LineCableModelsFEM}};
        options::NamedTuple = (;)
)
    return lock(FEM_SESSION_LOCK) do
        _compute_request(problem, formulation; options)
    end
end

function _compute_request(problem, formulation; options)
    problem = _preflight_fem_problem(problem)
    execution = _fem_computation_options(options)
    # Scalar and collection calls share completion notification and reuse rules.
    formulations = formulation isa LineCableModelsFEM ? [formulation] : formulation
    console = ConsoleLogger(stderr, Logging.Debug)
    logger = Engine.ConsoleVerbosityLogger(console, execution.verbosity)
    values = if execution.log_file === nothing
        with_logger(logger) do
            _compute_fem(problem, formulations, execution)
        end
    else
        mkpath(dirname(abspath(execution.log_file)))
        open(execution.log_file, "a") do io
            file_logger = SimpleLogger(io, Logging.Debug)
            with_logger(FEMTeeLogger(logger, file_logger)) do
                _compute_fem(problem, formulations, execution)
            end
        end
    end
    return formulation isa LineCableModelsFEM ? first(values) : values
end

function _compute_fem(
        problem::LineParametersProblem{Float64},
        formulations::AbstractVector{<:LineCableModelsFEM},
        execution::NamedTuple
)
    isempty(formulations) && throw(ArgumentError(
        "FEM formulation collections cannot be empty"))
    # Resolve and validate all requests before opening Gmsh or starting GetDP.
    models = [_resolved_fem_model(problem, formulation) for formulation in formulations]
    # The problem and loaded solver source are common to this batch. Only actual
    # material/mesh inputs and execution settings distinguish its calculations.
    keys = [JSON3.write(_fem_input_record(model, formulation))
        for (model, formulation) in zip(models, formulations)]
    first_result = _compute_fem(problem, first(formulations), execution, first(models))
    values = Vector{typeof(first_result)}(undef, length(formulations))
    values[1] = first_result
    execution.on_result === nothing || execution.on_result(problem, 1, first_result)
    completed = Dict(first(keys) => 1)
    for index in 2:length(formulations)
        formulation = formulations[index]
        previous = get(completed, keys[index], nothing)
        value = if previous === nothing || formulation.execution.ui ||
                   formulation.execution.mesh_policy === :remesh
            _compute_fem(problem, formulation, execution, models[index])
        else
            source = values[previous]
            @info "FEM reuses identical resolved inputs" formulation=index source_formulation=previous
            # Results remain independently mutable and each request keeps its own
            # selection record. The shared run record identifies the actual solve.
            metadata = merge(deepcopy(source.details),
                (formulations=formulation_record(formulation),))
            LineParameters(PhaseDomain,
                SeriesImpedance(copy(source.Z.values); basis=Engine.basis(source)),
                ShuntAdmittance(copy(source.Y.values); basis=Engine.basis(source)),
                copy(source.f), metadata)
        end
        typeof(value) === eltype(values) || throw(ArgumentError(
            "FEM formulations produced inconsistent result types"))
        values[index] = value
        completed[keys[index]] = index
        execution.on_result === nothing || execution.on_result(problem, index, value)
    end
    return values
end
