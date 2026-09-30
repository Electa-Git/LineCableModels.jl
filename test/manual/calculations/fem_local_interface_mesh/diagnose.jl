include("qualify.jl")
problem, form = fixture(:air,.1,0.)
dir = joinpath(ROOT,"air-f0.1-gamma0.0")
baseline = joinpath(dir,"baseline")
a = TOML.parsefile(joinpath(baseline,"mesh.toml"))
media = isempty(ARGS) ? (:air,:soil) : (Symbol(only(ARGS)),)
for medium in media
    label = "localized-$medium"
    candidate = mkpath(joinpath(dir,label))
    if !isfile(joinpath(candidate,"sources.toml"))
        for file in [readlines(joinpath(baseline,".onelab-export-files"));".onelab-export-files";"sources.toml"]
            mkpath(dirname(joinpath(candidate,file)))
            cp(joinpath(baseline,file),joinpath(candidate,file);force=true)
        end
    end
    b = mesh!(candidate,problem,true;media=medium)
    for key in ("pml_coordinates","path_coordinates","contour_coordinates")
        a[key] == b[key] || error("Preservation failure: $key")
    end
    b["nodes"] < a["nodes"] || error("Candidate does not reduce mesh size; no solve")
    solve!(candidate)
    compare!(dir,problem,form;candidate=label)
end
say("COMPLETE diagnostic ",media,"; baseline and rejected candidates retained")
