# Manual fixed-system diagnostic; never imported by feature code.
# Exactly two retained meshes, two solver configurations, two sequential sources.
using SHA, TOML, Dates

const REPO = normpath(joinpath(@__DIR__, "../../../.."))
const SAVED = joinpath(REPO, ".linecablemodels/fem/pml-corner-mesh/air-f0.1-gamma0.0")
const OUT = joinpath(REPO, ".linecablemodels/fem/pml-corner-algebraic")
const GETDP = "/home/amartins/.julia/artifacts/51e049a32eeb2ffebcc8a5945bb00f70470a3d29/getdp-3.5.0-Linux64/bin/getdp"
const MESHES = ("baseline", "metric-1.35")
const CONFIGURATIONS = ("original", "refine3")
# Diagnostics are identical in both configurations. ICNTL(11)=2 does not correct x.
const DIAGNOSTICS = ["-ksp_view", "-mat_mumps_icntl_4", "3", "-mat_mumps_icntl_11", "2"]
# Negative ICNTL(10) requests exactly three native corrections, without a stopping
# tolerance. CNTL(2) is left unchanged, recorded by MUMPS, and inactive in this mode.
const REFINEMENT = ["-mat_mumps_icntl_10", "-3"]
digest(p) = open(p) do io; bytes2hex(sha256(io)); end
say(xs...) = (println(Dates.now(), " ", xs...); flush(stdout))
record(p, x) = open(io -> TOML.print(io, x; sorted=true), p, "w")

function original_files()
    paths = String[]
    for mesh in MESHES
        for (dir, _, files) in walkdir(joinpath(SAVED, mesh))
            append!(paths, joinpath.(dir, files))
        end
    end
    for (dir, _, files) in walkdir(joinpath(REPO, "ext/LineCableModelsGmshExt"))
        append!(paths, joinpath.(dir, files))
    end
    push!(paths, joinpath(REPO, "test/manual/calculations/run_two_bare_wires_fem.jl"))
    return Dict(p => digest(p) for p in paths)
end

function prepare(mesh, config)
    source = joinpath(SAVED, mesh)
    dest = joinpath(OUT, mesh, config)
    ispath(dest) && error("Diagnostic directory already exists: $dest; no implicit rerun")
    mkpath(dest)
    for f in ("study.pro", "study_data.pro", "study.pre", "study.msh", "formulations")
        cp(joinpath(source, f), joinpath(dest, f))
    end
    # Save exact binary RHS and solution vectors in .res as b1,x1,b2,x2. The
    # temporary solution swap is restored before Solve/SolveAgain. Native Print
    # additionally writes each sparse matrix as a binary PETSc file.
    before = """
  CopySolution[Sys_FEM, "diagnostic_saved_x"];
  SetRightHandSideAsSolution[Sys_FEM];
  SaveSolution[Sys_FEM];
  CopySolution["diagnostic_saved_x", Sys_FEM];
  Test[$(raw"$FEMBasisTerminal") == 1]{
    Print[Sys_FEM, File "before1"];
  }{
    Print[Sys_FEM, File "before2"];
  }
"""
    after = """
SaveSolution[Sys_FEM];
Test[$(raw"$FEMBasisTerminal") == 1]{
  Print[Sys_FEM, File "after1"];
}{
  Print[Sys_FEM, File "after2"];
}
"""
    path = joinpath(dest, "formulations/quasi-full.pro")
    text = read(path, String)
    text = replace(text,
        "  Generate[Sys_FEM];\n" => "  Generate[Sys_FEM];\n" * before,
        "  GenerateRHSGroup[Sys_FEM, Terminals];\n" => "  GenerateRHSGroup[Sys_FEM, Terminals];\n" * before,
        "GetResidual[Sys_FEM, \$FEMResidualNorm];" => after * "GetResidual[Sys_FEM, \$FEMResidualNorm];")
    write(path, text)
    # -cal consumes the retained .pre; it does not reconstruct the gauge tree.
    cmd = Cmd(Cmd([GETDP, joinpath(dest, "study.pro"), "-msh", joinpath(dest, "study.msh"),
        "-cal", "-bin", "-setnumber", "FrequencyIndex", "1", "-setnumber", "PlotFieldMaps", "0",
        "-v", "4", "-nt", "1", "-ksp_diagonal_scale", "-ksp_diagonal_scale_fix",
        DIAGNOSTICS..., (config == "refine3" ? REFINEMENT : String[])...]); dir=dest)
    write(joinpath(dest, "command.txt"), repr(cmd) * "\n")
    record(joinpath(dest, "inputs.toml"), Dict(
        "mesh" => digest(joinpath(dest, "study.msh")),
        "preprocessing" => digest(joinpath(dest, "study.pre")),
        "data" => digest(joinpath(dest, "study_data.pro")),
        "saved_formulation" => digest(joinpath(source, "formulations/quasi-full.pro")),
        "instrumented_formulation" => digest(path),
        "configuration" => config, "vector_record_order" => ["b1", "x1", "b2", "x2"]))
    return dest, cmd
end

function main()
    mkpath(OUT)
    originals = original_files()
    record(joinpath(OUT, "preserved-inputs.toml"), originals)
    write(joinpath(OUT, "prescription.txt"), """
Two retained meshes: baseline and metric-1.35; air, 0.1 Hz, Gamma=0.
Two configurations: original defaults and ICNTL(10)=-3 (three fixed corrections).
Both source columns retain GenerateRHSGroup/SolveAgain factorization reuse.
Identical diagnostic settings: ICNTL(4)=3, ICNTL(11)=2, KSP view.
No other numerical settings, mesh changes, fixtures, retries or promotion.
Sparse residuals use exact saved double coefficients and binary solution vectors.
Higher precision accumulation does not improve assembly or solution precision.
Completion permits an inconclusive outcome; no residual-based EM certification.
""")
    jobs = [(mesh, config, prepare(mesh, config)...) for mesh in MESHES for config in CONFIGURATIONS]
    for (mesh, config, dest, cmd) in jobs
        say("BEGIN ", mesh, "/", config)
        t = time()
        # Unbuffer native C and Fortran output in the single parent live log.
        timed = Cmd(Cmd(vcat(["/usr/bin/time", "-f", "%e,%M", "-o",
            joinpath(dest, "resources.csv"), "stdbuf", "-oL", "-eL"], cmd.exec)); dir=dest)
        run(addenv(timed, "GFORTRAN_UNBUFFERED_ALL" => "y"))
        record(joinpath(dest, "execution.toml"), Dict("seconds_including_dumps" => time()-t,
            "completed_at" => string(Dates.now()), "preprocessing_unchanged" =>
            digest(joinpath(dest, "study.pre")) == digest(joinpath(SAVED, mesh, "study.pre"))))
        say("DONE ", mesh, "/", config, " seconds including dumps=", time()-t)
    end
    changed = [p for (p, hash) in originals if digest(p) != hash]
    record(joinpath(OUT, "preservation.toml"), Dict("changed_files" => changed, "files_checked" => length(originals)))
    isempty(changed) || error("Original files changed: $(join(changed, ", "))")
    say("COMPLETE four fixed native executions; eight source columns. Assess saved systems next.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
