include("qualify.jl")
problem, form = fixture(:mixed,1e6,.99)
dir = joinpath(ROOT,"mixed-f1.0e6-gamma0.99")
baseline = joinpath(dir,"baseline")
candidate = mkpath(joinpath(dir,"localized-wave-active"))
if !isfile(joinpath(candidate,"sources.toml"))
    for file in [readlines(joinpath(baseline,".onelab-export-files"));".onelab-export-files";"sources.toml"]
        mkpath(dirname(joinpath(candidate,file)))
        cp(joinpath(baseline,file),joinpath(candidate,file);force=true)
    end
end
a = TOML.parsefile(joinpath(baseline,"mesh.toml"))
b = mesh!(candidate,problem,true;media=:wave)
for key in ("pml_coordinates","path_coordinates","contour_coordinates")
    a[key]==b[key] || error("Preservation failure: $key")
end
say("CONSTANT-FIELD CONTROL: baseline nodes=",a["nodes"]," candidate nodes=",b["nodes"],
    " identical mesh bytes=",a["sha256"]==b["sha256"])
solve!(candidate)
compare!(dir,problem,form;candidate="localized-wave-active")
say("COMPLETE constant-field preservation control")
