function _performance_identity()
    return (
        julia_version = string(VERSION),
        kernel = string(Sys.KERNEL),
        architecture = string(Sys.ARCH),
        cpu = Sys.CPU_NAME,
        cpu_threads = Sys.CPU_THREADS,
        threads = Threads.nthreads(),
        optimization_level = Int(Base.JLOptions().opt_level),
        check_bounds = Int(Base.JLOptions().check_bounds),
        blas_threads = BLAS.get_num_threads(),
        blas = sprint(show, BLAS.get_config())
    )
end

function _compute_calculation(calculation; options=calculation.options)
    problem=calculation.problem
    if problem isa LineCableModels.ParametricProblem
        isempty(options) || (problem=LineCableModels.ParametricProblem(
            problem.space,merge(problem.options,options)))
        return compute(problem,calculation.formulation)
    end
    return isempty(options) ? compute(problem,calculation.formulation) :
        compute(problem,calculation.formulation;options)
end

function _source_timings(result)
    result isa ParametricResult && return (points=map(_source_timings,collect(result)),)
    retained=details(result)
    if haskey(retained,:execution) && haskey(retained.execution,:source_elapsed_seconds)
        source=retained.execution
        return (backend=:pscad,scope=source.source_elapsed_scope,
            seconds=source.source_elapsed_seconds,reused=get(source,:reused,false),
            adapter_wall_seconds=get(source,:wall_seconds,missing))
    elseif haskey(retained,:fem) && haskey(retained.fem,:timing)
        return retained.fem.timing
    end
    return (;)
end

function _result_reused(result; partial=true)
    result isa ParametricResult && return any(value->_result_reused(value;partial),result)
    retained=details(result)
    haskey(retained,:execution) && get(retained.execution,:reused,false) && return true
    return haskey(retained,:fem) && haskey(retained.fem,:timing) &&
        (get(retained.fem.timing,:reused,false) ||
            partial && get(retained.fem.timing,:recovered_columns,0)>0)
end

function _external_formulation(formulation)
    formulation isa Union{Engine.LineCableModelsFEM,PSCAD.PSCADFormulation} && return true
    formulation isa Union{Gridspace,AbstractVector} && return any(_external_formulation,formulation)
    formulation isa Union{LineCableModels.MonteCarlo,LineCableModels.LinearError,
        LineCableModels.Combinatorial} && return _external_formulation(formulation.inner)
    return false
end

function _performance_calculation(calculation)
    without_callback(options)=Base.structdiff(options,
        (;on_result=get(options,:on_result,nothing)))
    options=without_callback(calculation.options)
    problem=calculation.problem
    if problem isa LineCableModels.ParametricProblem
        problem=LineCableModels.ParametricProblem(problem.space,
            merge(without_callback(problem.options),options))
        options=(;)
    end
    return BenchmarkCalculation(calculation.id,problem,calculation.formulation;options)
end

function benchmark_local(
        case;
        options::NamedTuple = (;),
        samples::Int = 10,
        seconds::Real = 10
)
    calculation=_performance_calculation(BenchmarkCalculation(:local,
        case.problem,case.formulation;options))
    external=_external_formulation(calculation.formulation)
    reused=Ref(false)
    result=Ref{Any}()
    trial=_performance_span(;sample="local",samples,role=:candidate) do
        external || _compute_calculation(calculation)
        benchmark = BenchmarkTools.@benchmarkable $result[]=_compute_calculation($calculation) teardown=($reused[] |= _result_reused($result[])) samples=samples seconds=seconds evals=1
        Base.run(benchmark;warmup=false)
    end
    times = Float64.(trial.times) .* 1.0e-9
    return (
        minimum_seconds = minimum(times),
        median_seconds = median(times),
        bytes = trial.memory,
        allocations = trial.allocs,
        samples = length(times),
        scope = :compute_call_wall,
        reused = reused[],
        policy = (progress=false,diagnostics=:quiet,callbacks=false,
            warmup=external ? :native_not_repeated : :owned_call,
            allocation_scope=:julia, settings=_selection_value(calculation.options),
            workload=_numerical_record(calculation_record(calculation))),
        environment = _performance_identity()
    )
end

function performance_comparison(
        accepted,
        current,
        tolerance;
        instrumented::Bool = gauntlet_instrumented()
)
    same_environment = accepted.environment == current.environment && !instrumented
    comparable = same_environment &&
        get(accepted,:scope,nothing) === get(current,:scope,nothing) &&
        get(accepted,:policy,nothing) == get(current,:policy,nothing) &&
        !get(accepted,:reused,false) && !get(current,:reused,false)
    ratios = (
        median_time = current.median_seconds / accepted.median_seconds,
        bytes = accepted.bytes == 0 ? (current.bytes == 0 ? 1.0 : Inf) :
                current.bytes / accepted.bytes,
        allocations = accepted.allocations == 0 ?
                      (current.allocations == 0 ? 1.0 : Inf) :
                      current.allocations / accepted.allocations
    )
    passes = comparable ?
             ratios.median_time <= tolerance.median_time_ratio &&
             ratios.bytes <= tolerance.bytes_ratio &&
             ratios.allocations <= tolerance.allocations_ratio : nothing
    return (; comparable, ratios, passes)
end
