# Run from the REPL: include("dev/run_quasi_full.jl")
# Edit the inputs below, then include again. Every run gets a fresh directory.
import Pkg
Base.active_project() == normpath(joinpath(@__DIR__, "..", "gauntlet", "Project.toml")) ||
    Pkg.activate(joinpath(@__DIR__, "..", "gauntlet"))
using LineCableModels, Gmsh, LinearAlgebra, Printf
include(joinpath(@__DIR__, "quasi_full_paths.jl"))
FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)

frequencies = 10.0 .^ (-1:6)                    # Hz
radius = 0.0425                                # m; no insulation
positions = [(0.0, -1.0), (1.0, -1.0)]          # m
metal = Material(kind=:conductor, rho=1e-12)     # geometry/material input; PEC CLI below excludes its interior
wire = build(CableDesign, "bare", terminal(:core, core(metal; r=radius)))
system = build(LineCableSystem, [wire, wire], positions;
    connections=[Dict(:core=>1), Dict(:core=>2)], line_length=1.0)
problem = LineParametersProblem(system; frequencies,
    earth_props=homogeneous(rho=0.1, eps_r=1.0, mu_r=1.0))
formulation = Formulation(:LineCableModelsFEM;
    options=(physics=:quasi_fw, reduce_bundle=false, kron_reduction=false, ideal_transposition=false))
execution = computation_options(LineCableModelsFEM, (mesh_policy=:remesh, gmsh_verbosity=2, getdp_verbosity=3,
        plot_field_maps=false, solver_threads=1, keep_run_directory=true))

# Reuse the package's material resolver and mesher. No production FEM solve is
# called: the only equations executed below are the new quasi-full.pro.
model = FEM._resolved_fem_model(FEM._preflight_fem_problem(problem), formulation)
runtime_root = joinpath(pkgdir(LineCableModels), ".linecablemodels", "quasi-full")
manual_run = FEM._create_run(runtime_root)
run_directory = manual_run.path
path_files = String[]
lock(FEM.FEM_SESSION_LOCK) do
    session = FEM._start_gmsh(execution.gmsh_verbosity)
    try
        geometry = FEM._build_geometry!(model, "quasi-full-$(basename(run_directory))")
        global mesh_paths = FEM._select_meshes!(manual_run, model, geometry,
            execution, runtime_root)
        FEM._prepare_run_inputs!(manual_run, model)
        for (mesh, plan) in zip(mesh_paths, model.mesh_plans)
            path = joinpath(run_directory, "input", @sprintf("paths-f%04d.pro", plan.frequency_index))
            write_quasi_full_paths(path, mesh, plan, model, positions, radius)
            push!(path_files, path)
        end
    finally
        FEM._finish_gmsh(session)
    end
end

# model.pro selects the coupled file through its ONELAB Physics constant.
pro_file = joinpath(run_directory, "input", "getdp", "model.pro")
model_data = joinpath(run_directory, "input", "model_data.pro")
basis_file = joinpath(run_directory, "input", "bases.pro")
write(basis_file, "RequestedBases() = {1,2};\n")
getdp = FEM._getdp_selection(execution).path
Zqf = zeros(ComplexF64, 2, 2, length(frequencies))  # ohm/m
Mqf = similar(Zqf)                               # ohm m; raw inverse admittance
Peqf = similar(Zqf)                              # m/F
Yqf = similar(Zqf)                               # S/m

for (index, plan) in enumerate(model.mesh_plans)
    job = joinpath(run_directory, @sprintf("f%04d", index))
    mkpath(job)
    mesh, paths = mesh_paths[index], path_files[index]
    prefix = joinpath(job, "solver")
    maps = Int(execution.plot_field_maps)
    command = `$getdp $pro_file -solve LineCableModelsFEMScan
        -msh $mesh -name $prefix -v $(execution.getdp_verbosity)
        -setstring ModelDataPath $model_data -setstring RunDirectory $job
        -setstring BasisListPath $basis_file -setstring PathDataPath $paths
        -setnumber FrequencyIndex $index -setnumber FrequencyHz $(plan.frequency)
        -setnumber Val_Rint $(plan.domain_radius) -setnumber Val_Rext $(plan.shell_outer_radius)
        -setnumber PlotFieldMaps $maps -setnumber ReuseFactorization 1
        -setnumber Physics 1 -setnumber PerfectConductors 1`
    println("\n", command)
    # Native progress is visible in the REPL; raw columns/maps stay in job.
    run(addenv(Cmd(command; dir=job), "OMP_NUM_THREADS"=>"1", "OPENBLAS_NUM_THREADS"=>"1"))
    for (quantity, matrix) in (("Z", Zqf), ("P", Mqf)), basis in 1:2
        path = joinpath(job, "raw", "jobs", @sprintf("getdp-f%04d-b%04d-%s.tsv", index, basis, quantity))
        FEM._valid_job_raw(path, 2, index, plan.frequency, basis) || error("Invalid GetDP column: $path")
        for line in eachline(path)
            row = split(line)
            matrix[parse(Int, row[3]), basis, index] = complex(parse(Float64, row[5]), parse(Float64, row[6]))
        end
    end
    Yqf[:,:,index] = Mqf[:,:,index] \ Matrix{ComplexF64}(I, 2, 2)
    Peqf[:,:,index] = (2pi*im*plan.frequency) .* Mqf[:,:,index]
    println("f = ", plan.frequency, " Hz")
    for (label, matrix) in (("Z [ohm/m]", Zqf), ("Pe [m/F]", Peqf), ("Y [S/m]", Yqf))
        println(label)
        show(stdout, MIME"text/plain"(), matrix[:,:,index]); println()
    end
end
println("\nRun directory: ", run_directory)
# REPL results: system, problem, formulation, Zqf, Mqf, Peqf, Yqf, run_directory.
# PerfectConductors=1 gives earth/exterior matrices only. Set it to 0 to
# exercise the existing finite-metal a/u block with the chosen material.
nothing
