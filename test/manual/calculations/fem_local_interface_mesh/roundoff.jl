# Identical mesh, PDE, assembly and integration; only LU ordering changes.
include("qualify.jl")
problem, form = fixture(:air,.1,0.)
dir = joinpath(ROOT,"air-f0.1-gamma0.0")
baseline = joinpath(dir,"baseline")
candidate = mkpath(joinpath(dir,"ordering-amf"))
if !isfile(joinpath(candidate,"sources.toml"))
    for file in [readlines(joinpath(baseline,".onelab-export-files"));
            ".onelab-export-files";"sources.toml";"study.msh";"mesh.toml"]
        mkpath(dirname(joinpath(candidate,file)))
        cp(joinpath(baseline,file),joinpath(candidate,file);force=true)
    end
end
digest(joinpath(candidate,"study.msh"))==digest(joinpath(baseline,"study.msh")) || error("Mesh differs")
say("ROUND-OFF CONTROL: exact original mesh; MUMPS ordering 2 (AMF)")
solve!(candidate;ordering=2)
compare!(dir,problem,form;candidate="ordering-amf")
say("COMPLETE fixed-mesh ordering control")
