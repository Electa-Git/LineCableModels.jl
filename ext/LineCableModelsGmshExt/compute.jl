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

function _start_gmsh(verbosity::Int)
    owned = !Bool(gmsh.is_initialized())
    owned && Gmsh.initialize(String[]; finalize_atexit = false)
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

function _transition!(run::FEMRun, state::FEMRunState, message::AbstractString)
    run.state = state
    run.message = String(message)
    _write_json_atomic(joinpath(run.path, "run.json"),
        (
            schema = "LineCableModels.FEMRun",
            version = 1,
            state = string(state),
            message = run.message,
            mesh_source = String(run.mesh_source),
            mesh_fingerprint = run.mesh_fingerprint,
            getdp_invocations = run.getdp_invocations,
            updated_unix_seconds = time()
        ))
    return state
end

function _prepare_run_inputs!(run::FEMRun, model::FEMResolvedModel)
    problem_path = joinpath(run.path, "input", "problem.json")
    model_data_path = joinpath(run.path, "input", "model_data.pro")
    _write_problem_snapshot(problem_path, model.problem)
    _write_model_data(model_data_path, model)
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
    allowed = (:verbosity, :output_basis, :trace, :log_file)
    unknown = filter(key -> key ∉ allowed, keys(options))
    isempty(unknown) || throw(ArgumentError(
        "unknown LineCableModelsFEM computation options: $(sort!(collect(unknown)))",
    ))
    log_file = get(options, :log_file, nothing)
    log_file isa Union{Nothing, AbstractString} || throw(ArgumentError(
        "log_file must be a path string or nothing",
    ))
    log_file === "" && throw(ArgumentError("log_file cannot be empty"))
    standard_keys = filter(!=(:log_file), keys(options))
    standard_values = map(key -> getproperty(options, key), standard_keys)
    standard = NamedTuple{Tuple(standard_keys)}(Tuple(standard_values))
    execution = computation_options(
        Val(LineCableModels.LineCableModelsCoaxial), standard
    )
    return (; execution..., log_file = log_file === nothing ? nothing : String(log_file))
end

function _headless_solve!(
        run::FEMRun,
        model::FEMResolvedModel,
        formulation::LineCableModelsFEM,
        execution::NamedTuple,
        runtime_root::String
)
    geometry = _build_geometry!(
        model, "LineCableModelsFEM-$(basename(run.path))"
    )
    _transition!(run, geometry_ready, "geometry ready")
    @info "FEM geometry ready" run_directory=run.path
    mesh_path = _select_mesh!(run, model, geometry, formulation, runtime_root)
    _transition!(run, mesh_ready, "mesh ready")
    @info "FEM mesh ready" source=run.mesh_source fingerprint=run.mesh_fingerprint
    model_data_path = _prepare_run_inputs!(run, model)
    _publish_transport!(run, model, model_data_path, mesh_path, formulation)
    _transition!(run, running, "GetDP scan running")
    @info "Starting one-pass GetDP frequency/source scan"
    _run_getdp!(run, model, formulation, mesh_path)
    scan = _parse_scan(run, model, formulation)
    parameters = _line_parameters(run, model, formulation, execution, scan)
    gmsh.onelab.set_number(_onelab_name("completion_status"), [1.0])
    _transition!(run, completed, "results validated")
    @info "FEM scan completed successfully"
    return parameters
end

function _ui_solve!(
        run::FEMRun,
        model::FEMResolvedModel,
        formulation::LineCableModelsFEM,
        execution::NamedTuple,
        runtime_root::String
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
    mesh_path = nothing
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
            mesh_path = _select_mesh!(run, model, geometry, formulation, runtime_root)
            model_data_path = _prepare_run_inputs!(run, model)
            _publish_transport!(
                run, model, model_data_path, mesh_path, formulation
            )
            _transition!(run, mesh_ready, "mesh ready")
            state = _set_ui_status(:mesh_ready)
        elseif transition === :solve_requested
            mesh_path === nothing && begin
                _set_ui_status(:mesh_required)
                gmsh.fltk.wait(0.05)
                continue
            end
            state = _set_ui_status(:running)
            _transition!(run, running, "GetDP scan running")
            _run_getdp!(run, model, formulation, mesh_path)
            scan = _parse_scan(run, model, formulation)
            parameters = _line_parameters(run, model, formulation, execution, scan)
            formulation.execution.plot_field_maps && _merge_maps!(scan.map_paths)
            gmsh.onelab.set_number(_onelab_name("completion_status"), [1.0])
            _transition!(run, completed, "results ready")
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
        problem::LineParametersProblem,
        formulation::LineCableModelsFEM,
        execution::NamedTuple
)
    model = _resolved_fem_model(problem, formulation)
    runtime_root = _runtime_root()
    run = _create_run(runtime_root)
    session = nothing
    succeeded = false
    try
        session = _start_gmsh(formulation.execution.gmsh_verbosity)
        parameters = formulation.execution.ui ?
                     _ui_solve!(run, model, formulation, execution, runtime_root) :
                     _headless_solve!(run, model, formulation, execution, runtime_root)
        succeeded = true
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
        if succeeded && !formulation.execution.keep_run_directory
            expected_parent = realpath(joinpath(runtime_root, "runs"))
            run_parent = realpath(dirname(run.path))
            run_parent == expected_parent || error(
                "refusing to remove FEM run outside the runtime root"
            )
            rm(run.path; recursive = true, force = true)
        end
    end
end

function compute(
        problem::LineParametersProblem,
        formulation::LineCableModelsFEM;
        options::NamedTuple = (;)
)
    execution = _fem_computation_options(options)
    console = ConsoleLogger(stderr, Logging.Debug)
    logger = Engine.ConsoleVerbosityLogger(console, execution.verbosity)
    if execution.log_file === nothing
        return with_logger(logger) do
            _compute_fem(problem, formulation, execution)
        end
    end
    mkpath(dirname(abspath(execution.log_file)))
    return open(execution.log_file, "a") do io
        file_logger = SimpleLogger(io, Logging.Debug)
        with_logger(FEMTeeLogger(logger, file_logger)) do
            _compute_fem(problem, formulation, execution)
        end
    end
end
