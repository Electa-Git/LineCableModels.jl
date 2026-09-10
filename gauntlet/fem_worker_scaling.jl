# Run from the repository root with --project=gauntlet.
# Uses eight frozen detailed-cable meshes; writes only fresh validation runs.
using LineCableModels, Gmsh, JSON3, Printf, Dates, LinearAlgebra
const E = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
const source=joinpath(pwd(), ".linecablemodels/fem/runs/run-orkY1H")
const original=LineCableModels.ImportExport.deserialize_value(JSON3.read(read(
    joinpath(source, "input/problem.json"), String)))
const recorded=JSON3.read(read(joinpath(source, "input/computation.json"), String))
const indices=[1, 15, 29, 43, 58, 72, 86, 101]
const problem=LineParametersProblem(
    original.system; frequencies = original.frequencies[indices],
    temperature = original.temperature, earth_props = original.earth_props)
const options=(; (Symbol(k)=>v for (k, v) in pairs(recorded.options))...)
const root=joinpath(pwd(), ".linecablemodels/fem/scaling-"*Dates.format(now(), "yyyymmdd-HHMMSS"))
const reports=Any[]
reference=nothing
for workers in (1, 2, 4, 8)
    form=Formulation(:LineCableModelsFEM; options,
        fem_options = (
            getdp_verbosity = 4, gmsh_verbosity = 0, frequency_workers = workers,
            solver_threads = 1, keep_run_directory = true))
    model=E._resolved_fem_model(problem, form)
    for (current, old) in zip(model.material_plans, recorded.materials)
        @assert real.(current.admittivity)==Float64.(old.sigma[indices])
        @assert imag.(current.admittivity)==Float64.(old.omega_epsilon[indices])
    end
    for (plan, index) in zip(model.mesh_plans, indices)
        @assert plan.domain_radius==recorded.mesh_plans[index].domain_radius
        @assert plan.shell_outer_radius==recorded.mesh_plans[index].shell_outer_radius
    end
    run=E._create_run(root)
    E._write_json_atomic(joinpath(run.path, "input/computation.json"), E._fem_input_record(model, form))
    E._prepare_run_inputs!(run, model)
    meshes=[joinpath(source, "mesh", i==101 ? "model.msh" :
                                     @sprintf("frequency_%04d.msh", i)) for i in indices]
    println("SCALING_START workers=", workers, " run=", run.path);
    flush(stdout)
    peak_worker=Ref(0.0);
    peak_aggregate=Ref(0.0);
    sampled=Ref(0.0)
    function sample_memory()
        Sys.islinux() && time()-sampled[]>0.1 || return true
        sampled[]=time();
        aggregate=0.0
        attempts=joinpath(run.path, "attempts")
        isdir(attempts) || return true
        for directory in readdir(attempts; join = true)
            rss=try
                manifest=JSON3.read(read(joinpath(directory, "attempt.json"), String))
                manifest.state=="running" || continue
                text=read("/proc/$(manifest.pid)/status", String)
                parse(Float64, match(r"VmRSS:\s+(\d+)\s+kB", text)[1])/1024
            catch
                continue
            end
            peak_worker[]=max(peak_worker[], rss);
            aggregate+=rss
        end
        peak_aggregate[]=max(peak_aggregate[], aggregate)
        return true
    end
    started=time()
    E._run_getdp!(run, model, form, meshes; pump = sample_memory)
    elapsed=time()-started
    scan=E._parse_scan(run, model, form)
    global reference=reference===nothing ? scan : reference
    @assert scan.Z==reference.Z && scan.P==reference.P
    E._write_scan_checksums(run, scan);
    E._transition!(run, E.completed, "worker scaling validation completed")
    peak=0.0
    events=Dict(key=>0
    for key in ("MatLUFactorSym", "MatLUFactorNum", "MatSolve", "KSPSolve"))
    for attempt in readdir(joinpath(run.path, "attempts"); join = true)
        for found in
            eachmatch(r"Mem\s*=\s*([0-9.]+)\s*Mb", read(joinpath(attempt, "getdp.log"), String))
            peak=max(peak, parse(Float64, found[1]))
        end
        for line in readlines(joinpath(attempt, "petsc.log"))
            fields=split(line)
            isempty(fields) && continue
            haskey(events, fields[1]) && (events[fields[1]]+=parse(Int, fields[2]))
        end
    end
    @assert events["MatLUFactorSym"]==events["MatLUFactorNum"]==8
    @assert events["MatSolve"]==72
    report=(; workers, solver_threads = 1, indices, elapsed_seconds = elapsed,
        exact_serial_agreement = true, peak_reported_worker_mb = peak,
        sampled_peak_worker_rss_mib = peak_worker[], sampled_peak_aggregate_rss_mib = peak_aggregate[], events, run = run.path)
    push!(reports, report);
    E._write_json_atomic(joinpath(root, "report.json"), reports)
    println(JSON3.write(report));
    flush(stdout)
end
