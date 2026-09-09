# Development validation against the four retained 18 kV campaigns.
# Run from the repository root with --project=gauntlet and root on LOAD_PATH.
# Existing campaign inputs and meshes are read only; every execution gets a fresh root.
using LineCableModels, Gmsh, JSON3, LinearAlgebra, Printf, Dates, SHA
const E = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
const ROOT = joinpath(pwd(), ".linecablemodels", "fem", "validation-" *
                                                        Dates.format(now(), "yyyymmdd-HHMMSS"))
mkpath(ROOT)
println("VALIDATION_ROOT ", ROOT);
flush(stdout)
const REPORTS = Any[]

relative(a, b) = norm(a-b) / max(norm(b), eps(Float64))
function read_frequency(run, f, n)
    matrices = [zeros(ComplexF64, n, n) for _ in 1:2]
    for basis in 1:n
        paths = E._column_paths(run.path, f, basis, false)
        for (matrix, path) in zip(matrices, (paths.Z, paths.P))
            for row in filter(!isempty, strip.(readlines(path)))
                fields=split(row, '\t')
                @assert parse(Int, fields[1])==f && parse(Int, fields[4])==basis
                matrix[parse(Int, fields[3]), basis]=complex(
                    parse(Float64, fields[5]), parse(Float64, fields[6]))
            end
        end
    end
    return matrices
end

function validate_campaign(name, method)
    source=joinpath(pwd(), ".linecablemodels/fem/runs", name)
    problem=LineCableModels.ImportExport.deserialize_value(JSON3.read(read(
        joinpath(source, "input/problem.json"), String)))
    recorded=JSON3.read(read(joinpath(source, "input/computation.json"), String))
    options=(; (Symbol(k)=>v for (k, v) in pairs(recorded.options))...)
    form=Formulation(:LineCableModelsFEM; insulation_admittance = method,
        semicon_admittance = method, options,
        fem_options = (getdp_verbosity = 4, gmsh_verbosity = 0, frequency_workers = 4,
            solver_threads = 1, keep_run_directory = true))
    model=E._resolved_fem_model(problem, form)
    run=E._create_run(ROOT)
    inputs=E._fem_input_record(model, form)
    E._write_json_atomic(joinpath(run.path, "input/computation.json"), inputs)
    E._prepare_run_inputs!(run, model)
    @assert read(joinpath(run.path, "input/model_data.pro"))==read(joinpath(source, "input/model_data.pro"))
    @assert E._resume_value_matches(recorded.mesh_plans, JSON3.read(JSON3.write(model.mesh_plans)))
    nf=length(problem.frequencies);
    n=length(model.terminal_ids)
    meshes=[joinpath(source, "mesh", f==nf ? "model.msh" :
                                     @sprintf("frequency_%04d.msh", f)) for f in 1:nf]
    oldrun=E.FEMRun(source, E.completed, "retained reference", :retained, "")
    reference=[E._parse_raw_matrix(Float64, joinpath(source, "raw", q*".tsv"),
                   problem.frequencies, n, oldrun) for q in ("Z", "P")]
    checked=falses(nf)
    errors=Dict(key=>0.0
    for key in ("Z_relative", "P_relative", "Y_primitive_relative",
        "Z_reduced_relative", "P_reduced_relative", "Y_reduced_relative",
        "Z_absolute", "P_absolute", "Y_reduced_absolute"))
    residual=Ref(0.0);
    condition=Ref(0.0);
    entry_count=Ref(0)
    started=time();
    last_print=Ref(0.0)
    function check_progress()
        for f in 1:nf
            checked[f] && continue
            all(b->isfile(E._column_paths(run.path, f, b, false).checkpoint), 1:n) ||
                continue
            z, p=read_frequency(run, f, n);
            zr=reference[1][:, :, f];
            pr=reference[2][:, :, f]
            for (key, a, b) in (("Z", z, zr), ("P", p, pr))
                errors[key * "_relative"]=max(errors[key * "_relative"], relative(a, b))
                errors[key * "_absolute"]=max(errors[key * "_absolute"], maximum(abs.(a-b)))
            end
            errors["Y_primitive_relative"]=max(errors["Y_primitive_relative"], relative(inv(p), inv(pr)))
            reduced=LineCableModels.Engine.reduce_primitive_matrices(
                reshape(z, n, n, 1), reshape(p, n, n, 1),
                problem.system.connection_order, form.options)
            ref_reduced=LineCableModels.Engine.reduce_primitive_matrices(
                reshape(zr, n, n, 1), reshape(pr, n, n, 1),
                problem.system.connection_order, form.options)
            inversion=LineCableModels.Engine.potential_to_admittance(reduced.P; diagnostics = true)
            ref_inversion=LineCableModels.Engine.potential_to_admittance(ref_reduced.P; diagnostics = true)
            for (key, a, b) in (("Z_reduced", reduced.Z, ref_reduced.Z),
                ("P_reduced", reduced.P, ref_reduced.P),
                ("Y_reduced", inversion.Y, ref_inversion.Y))
                errors[key * "_relative"]=max(errors[key * "_relative"], relative(a, b))
            end
            errors["Y_reduced_absolute"]=max(
                errors["Y_reduced_absolute"], maximum(abs.(inversion.Y-ref_inversion.Y)))
            residual[]=max(residual[], maximum(inversion.residuals))
            condition[]=max(condition[], maximum(inversion.condition_numbers))
            checked[f]=true;
            entry_count[]+=2n*n+length(inversion.Y)
        end
        if time()-last_print[]>=15 || all(checked)
            status=(source = name, run = run.path, checked_frequencies = count(checked),
                total_frequencies = nf, checked_entries = entry_count[], errors,
                max_inversion_residual = residual[], max_reduced_P_condition = condition[], elapsed_seconds = time()-started,
                launches = run.getdp_invocations, complete = all(checked),
                passed = all(checked)&&maximum(v
                for (k, v) in errors if endswith(k, "relative"))<1e-8)
            E._write_json_atomic(joinpath(ROOT, name*"-progress.json"), status)
            println(JSON3.write(status));
            flush(stdout);
            last_print[]=time()
        end
        return true
    end
    E._transition!(run, E.running, "full retained-mesh numerical validation")
    E._run_getdp!(run, model, form, meshes; pump = check_progress)
    check_progress()
    scan=E._parse_scan(run, model, form)
    parameters=E._line_parameters(
        run, model, form, E._fem_computation_options((trace = true,)), scan, inputs)
    E._write_scan_checksums(run, scan);
    E._transition!(run, E.completed, "full numerical validation completed")
    @assert all(checked)
    report=JSON3.read(read(joinpath(ROOT, name*"-progress.json"), String), Dict{
        String, Any})
    report["elapsed_seconds"]=time()-started
    report["mesh_sha256"]=[bytes2hex(open(sha256, path)) for path in meshes]
    push!(REPORTS, report)
    E._write_json_atomic(joinpath(ROOT, "report.json"), REPORTS)
    @assert report["passed"] "Numerical mismatch against retained campaign"
    println("CAMPAIGN_PASS ", name, " ", JSON3.write(errors));
    flush(stdout)
end

for (name, method) in (("run-YiqIRh", :default), ("run-aZVWtq", :Ametani2004),
    ("run-Pbj2uz", :Ametani2004), ("run-orkY1H", :default))
    validate_campaign(name, method)
end
println("ALL_FOUR_CAMPAIGNS_PASS ", ROOT);
flush(stdout)
