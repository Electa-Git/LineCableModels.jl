# Diagnose linear-solver sensitivity on the EXACT frozen reference mesh.
# No meshing, field equations, extraction paths or reference files are changed.
isdefined(@__MODULE__,:ROOT) || include("qualify.jl")
problem, form = fixture(:mixed,1e6,.99)
dir = joinpath(ROOT,"mixed-f1.0e6-gamma0.99")
baseline = joinpath(dir,"baseline")
controls = (
    ("solver-no-reuse",["-setnumber","ReuseFactorization","0"]),
    ("solver-amf",["-mat_mumps_icntl_7","2"]),
    ("solver-gmres",["-ksp_type","gmres","-ksp_rtol","1e-12","-ksp_atol","1e-14",
        "-ksp_max_it","20","-ksp_norm_type","unpreconditioned","-ksp_monitor_true_residual",
        "-ksp_converged_reason"]),
    ("solver-equilibrated",["-mat_mumps_icntl_6","7","-mat_mumps_icntl_8","7",
        "-mat_mumps_icntl_10","5","-mat_mumps_icntl_11","2","-ksp_view"]),
    ("solver-refined",["-mat_mumps_icntl_10","10","-mat_mumps_cntl_2","1e-16",
        "-mat_mumps_icntl_11","1","-ksp_view"]))
for (label,args) in controls
    isempty(ARGS) || label in ARGS || continue
    target = joinpath(dir,label)
    if !isfile(joinpath(target,"sources.toml"))
        for file in [readlines(joinpath(baseline,".onelab-export-files"));
                ".onelab-export-files";"sources.toml";"study.msh";"mesh.toml"]
            mkpath(dirname(joinpath(target,file)))
            cp(joinpath(baseline,file),joinpath(target,file))
        end
    end
    digest(joinpath(target,"study.msh"))==digest(joinpath(baseline,"study.msh")) || error("Control mesh changed")
    say("FINITE GAMMA FIXED-MESH CONTROL ",label)
    solve!(target;extra_args=args)
    compare!(dir,problem,form;candidate=label,tolerance=.02,record_prefix="diagnostic-")
    for m in eachmatch(r"FEM algebraic residual:[^\n]+",read(joinpath(target,"solver.log"),String))
        say(label," ",m.match)
    end
end
say("COMPLETE finite-Gamma fixed-mesh solver controls")
