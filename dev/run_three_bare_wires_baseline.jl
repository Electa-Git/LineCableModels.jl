# Load in the IDE or run with julia --startup-file=no dev/run_three_bare_wires_baseline.jl.
# Uses the active environment, without activating/installing anything. Every run
# creates a new attempt; retained results are never replaced. No plotting needed.
using LineCableModels
using Gmsh
using Printf, TOML, SHA, Dates, Logging, Pkg

isdefined(@__MODULE__, :CurrentScenarios) ||
    include(joinpath(@__DIR__, "..", "test", "support", "scenarios.jl"))

function baseline_toml(path, record)
    open(io -> TOML.print(io, record; sorted = true), path, "w")
end

function baseline_checksums(directory)
    return Dict(relpath(joinpath(root, file), directory) => open(bytes2hex ∘ sha256, joinpath(root, file))
    for (root, _, files) in walkdir(directory) for file in files)
end

function baseline_write_result(directory, result, expected_frequencies)
    f = LineCableModels.frequencies(result)
    labels = String.(details(result).data.coordinates)
    @assert f == expected_frequencies
    @assert labels == ["cable:$i:core" for i in 1:3]
    @assert length(unique(labels)) == 3
    for (name, selector, unit) in (("Z", Z, "ohm_per_m"), ("Y", Y, "S_per_m"))
        values = observe(result, selector)
        @assert size(values) == (3, 3, 9)
        @assert all(isfinite, values)
        path = joinpath(directory, "$name.csv")
        open(path, "w") do io
            println(io, "frequency_hz,response_terminal,basis_terminal,real_$unit,imaginary_$unit")
            for k in eachindex(f), i in 1:3, j in 1:3
                value = values[i, j, k]
                @printf(io, "%.17g,%s,%s,%.17g,%.17g\n",
                    f[k], labels[i], labels[j], real(value), imag(value))
            end
        end
        # Check every stored entry and its coordinates, not only the dimensions.
        rows = readlines(path)[2:end]
        @assert length(rows) == 81
        for (row, (k, i, j)) in
            zip(rows, ((k, i, j) for k in eachindex(f) for i in 1:3 for j in 1:3))
            fs, li, lj, re, im = split(row, ',')
            @assert parse(Float64, fs) == f[k] && li == labels[i] && lj == labels[j]
            @assert complex(parse(Float64, re), parse(Float64, im)) == values[i, j, k]
        end
    end
    return labels
end

function compare_three_bare_wires_baseline(attempt)
    summary = Dict{String, Any}("status"=>"differences only; no acceptance verdict",
        "near_zero_rule"=>"abs(component) <= eps(Float64)*maximum(abs, FEM channel); relative error is NaN there")
    for placement in keys(CurrentScenarios.three_bare_wires_layouts), quantity in ("Z", "Y")

        name = "three_bare_wires_$placement"
        paths = [joinpath(attempt, name, backend, "$quantity.csv")
                 for backend in ("unified", "fem")]
        if !all(isfile, paths)
            summary["$name/$quantity"] = Dict("status"=>"unavailable; inspect scan status in manifest")
            continue
        end
        rows = map(path -> split.(readlines(path)[2:end], ','), paths)
        @assert length(rows[1]) == length(rows[2]) == 81
        values = map(rs -> [complex(parse(Float64, r[4]), parse(Float64, r[5])) for r in rs], rows)
        candidate, reference = values
        near_zero = eps(Float64)*maximum(abs, reference)
        relative(delta, ref) = abs(ref) <= near_zero ? NaN : abs(delta)/abs(ref)
        errors = abs.(candidate-reference)
        worst = argmax(errors)
        relative_errors = [relative(c-r, r) for (c, r) in zip(candidate, reference)]
        resolved_relative = filter(isfinite, relative_errors)
        worst_relative = isempty(resolved_relative) ? nothing :
                         argmax([isfinite(x) ? x : -Inf for x in relative_errors])
        summary["$name/$quantity"] = Dict("max_absolute_error"=>maximum(errors),
            "max_relative_error"=>isempty(resolved_relative) ? NaN :
                                  maximum(resolved_relative),
            "worst_absolute_entry"=>join(rows[1][worst][1:3], ","),
            "relative_at_worst_absolute_entry"=>relative_errors[worst],
            "reference_magnitude_at_worst_absolute_entry"=>abs(reference[worst]),
            "worst_relative_entry"=>worst_relative === nothing ? "unresolved" :
                                    join(rows[1][worst_relative][1:3], ","),
            "absolute_at_worst_relative_entry"=>worst_relative === nothing ? NaN :
                                                errors[worst_relative],
            "reference_magnitude_at_worst_relative_entry"=>worst_relative === nothing ?
                                                           NaN :
                                                           abs(reference[worst_relative]),
            "near_zero_threshold"=>near_zero)
        # Near-zero components retain their absolute difference; no arbitrary
        # denominator floor turns them into apparently small relative errors.
        open(joinpath(attempt, name, "$(quantity)-differences.csv"), "w") do io
            println(io,
                "frequency_hz,response_terminal,basis_terminal,absolute_error,relative_error,real_absolute_error,real_relative_error,imaginary_absolute_error,imaginary_relative_error,real_reference_near_zero,imaginary_reference_near_zero")
            for (i, (c, r)) in enumerate(zip(candidate, reference))
                @assert rows[1][i][1:3] == rows[2][i][1:3]
                d = c-r
                @printf(io, "%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%s,%s\n",
                    join(rows[1][i][1:3], ","), abs(d), relative(d, r),
                    abs(real(d)), relative(real(d), real(r)),
                    abs(imag(d)), relative(imag(d), imag(r)),
                    string(abs(real(r))<=near_zero), string(abs(imag(r))<=near_zero))
            end
        end
    end
    baseline_toml(joinpath(attempt, "comparison.toml"), summary)
    return summary
end

function run_three_bare_wires_baseline(; analytical = true, fem = true)
    repository = normpath(joinpath(@__DIR__, ".."))
    root = joinpath(repository, "test", "fixtures", "reference", "three_bare_wires")
    mkpath(root)
    attempt = mktempdir(root; prefix = "capture-$(Dates.format(now(UTC), "yyyymmddTHHMMSS"))-", cleanup = false)
    println("Baseline directory: ", attempt)
    flush(stdout)
    write(joinpath(attempt, "source.patch"), read(`git -C $repository diff HEAD --binary`, String))
    write(joinpath(attempt, "source-status.txt"), read(`git -C $repository status --short`, String))
    # Before: capture required an unrelated historical local plan. Now retain
    # the executable fixture and runner; source.patch records the actual code.
    for source in ("test/support/scenarios.jl", "dev/run_three_bare_wires_baseline.jl")
        target = joinpath(attempt, "sources", source)
        mkpath(dirname(target))
        cp(joinpath(repository, source), target)
    end
    packages = Dict(string(uuid) => Dict("name"=>dep.name,
                        "version"=>string(dep.version), "source"=>something(dep.source, ""))
    for (uuid, dep) in Pkg.dependencies())
    manifest = Dict{String, Any}("status"=>"capturing", "julia"=>string(VERSION),
        "project"=>something(Base.active_project(), ""),
        "revision"=>strip(read(`git -C $repository rev-parse HEAD`, String)),
        "packages"=>packages, "runs"=>Dict{String, Any}(),
        "matching_criteria"=>"Pending human review; no numerical acceptance tolerance imposed.")
    baseline_toml(joinpath(attempt, "manifest.toml"), manifest)
    reductions = (
        reduce_bundle = false, kron_reduction = false, ideal_transposition = false)
    selections = Pair{String, Any}[]
    analytical && push!(selections,
        "unified" => (
            Formulation(earth_impedance = formula(:unified),
                earth_admittance = formula(:unified);
                options = reductions), (trace = true, output_basis = :pul)))
    fem && push!(selections,
        "fem" => (
            Formulation(:LineCableModelsFEM;
                options = merge((physics = :quasi_tem,), reductions)),
            (ui = false, mesh_policy = :remesh, resume_run_directory = nothing,
                keep_run_directory = true, trace = true, output_basis = :pul,
                verbosity = (default = 1,), gmsh_verbosity = 2, getdp_verbosity = 4,
                frequency_workers = 2, solver_threads = 4, plot_field_maps = false)))
    # Each public call is a complete fresh scan, not a permutation or a resumed solve.
    for (backend, (formulation, options)) in selections,
        (placement, heights) in pairs(CurrentScenarios.three_bare_wires_layouts)

        name = "three_bare_wires_$placement"
        directory = joinpath(attempt, name, backend)
        mkpath(directory)
        problem = CurrentScenarios.three_bare_wires_problem(; heights, name)
        export_data(:json, problem; file_name = joinpath(directory, "problem.json"))
        write(joinpath(directory, "formulation.txt"), repr(NamedTuple(formulation)))
        write(joinpath(directory, "computation-options.txt"), repr(options))
        record = Dict{String, Any}("status"=>"running", "directory"=>relpath(directory, attempt))
        manifest["runs"]["$name/$backend"] = record
        baseline_toml(joinpath(attempt, "manifest.toml"), manifest)
        println("Starting $name / $backend / 9 frequencies")
        flush(stdout)
        started = time()
        open(joinpath(directory, "compute.log"), "w") do log
            redirect_stdio(stdout = log, stderr = log) do
                with_logger(ConsoleLogger(log)) do
                    try
                        result = @time compute(problem, formulation; options)
                        record["seconds"] = time()-started
                        record["terminal_order"] = baseline_write_result(directory, result, problem.frequencies)
                        write(joinpath(directory, "details.txt"), repr(details(result).data))
                        if backend == "fem"
                            run_directory = details(result).data.fem.run.run_directory
                            record["fem_run_directory"] = run_directory
                            baseline_toml(joinpath(directory, "fem-artifacts.toml"),
                                baseline_checksums(run_directory))
                        end
                        record["status"] = "captured"
                    catch exception
                        record["seconds"] = time()-started
                        record["status"] = "failed"
                        record["error"] = sprint(showerror, exception)
                        showerror(log, exception, catch_backtrace())
                        println(log)
                        # Retain the backend's failed work as well as successful solves.
                        if hasproperty(exception, :run_directory) &&
                           exception.run_directory !== nothing
                            record["fem_run_directory"] = exception.run_directory
                            baseline_toml(joinpath(directory, "fem-artifacts.toml"),
                                baseline_checksums(exception.run_directory))
                        end
                    end
                end
            end
        end
        record["sha256"] = baseline_checksums(directory)
        baseline_toml(joinpath(attempt, "manifest.toml"), manifest)
        println("$(record["status"]): $name / $backend ($(round(record["seconds"]; digits=2)) s)")
        flush(stdout)
    end
    manifest["status"] = all(run["status"] == "captured"
    for run in values(manifest["runs"])) ?
                         "captured; pending comparison and human review" :
                         "incomplete; failures retained"
    baseline_toml(joinpath(attempt, "manifest.toml"), manifest)
    compare_three_bare_wires_baseline(attempt)
    return attempt
end

if abspath(PROGRAM_FILE) == (@__FILE__) || isinteractive()
    baseline_directory = run_three_bare_wires_baseline()
end
