# Qualification only. Frozen references are read, never remeshed or solved.
# Run serially; the wider prescription is tested on a failing case first.
include("qualify.jl")
const BOTH_CASES = [(:mixed,.1,0.); [(layout,f,0.) for layout in (:air,:soil,:mixed)
    for f in (.1,1e3,1e6) if !(layout==:mixed && f==.1)];
    [(layout,f,.99) for layout in (:mixed,:air,:soil) for f in (.1,1e6)]]

function both_case!(layout, frequency, fraction, media)
    dir = joinpath(ROOT,"$layout-f$(frequency)-gamma$(fraction)")
    problem, form = fixture(layout,frequency,fraction)
    baseline = joinpath(dir,"baseline")
    a = TOML.parsefile(joinpath(baseline,"mesh.toml"))
    digest(joinpath(baseline,"study.msh")) == a["sha256"] || error("Frozen baseline mesh changed")
    for (file,hash) in TOML.parsefile(joinpath(baseline,"sources.toml"))
        digest(joinpath(baseline,file)) == hash || error("Frozen baseline source changed: $file")
    end
    isfile(joinpath(baseline,"results/f0001-quasi-fw-b0000/completed.txt")) || error("Missing frozen result")
    candidate = media==:both ? "localized" : "localized-decay"
    target = joinpath(dir,candidate)
    if !isfile(joinpath(target,"sources.toml"))
        for file in [readlines(joinpath(baseline,".onelab-export-files"));".onelab-export-files";"sources.toml"]
            mkpath(dirname(joinpath(target,file)))
            cp(joinpath(baseline,file),joinpath(target,file))
        end
    end
    b = mesh!(target,problem,true;media)
    for key in ("pml_coordinates","path_coordinates","contour_coordinates")
        a[key] == b[key] || error("Preservation failure $(basename(dir)): $key")
    end
    solve!(target)
    passed = compare!(dir,problem,form;candidate,tolerance=.02,record_prefix="two-percent-")
    say("BOTH MEDIA ",basename(dir)," candidate=",candidate," passed=",passed,
        " nodes ",a["nodes"]," -> ",b["nodes"]," physical triangles ",
        a["air_triangles"]+a["soil_triangles"]," -> ",b["air_triangles"]+b["soil_triangles"])
    return passed
end

function qualify_both()
say("BEGIN both-media qualification; 2% per R/X/G/B entry; no new signs; frozen reference FEM")
failed = nothing
for case in BOTH_CASES
    if !both_case!(case...,:both)
        failed = case
        break
    end
end
if failed === nothing
    record(joinpath(ROOT,"both-media-selection.toml"),Dict("candidate"=>"localized","passed"=>true,"cases"=>length(BOTH_CASES),"relative_tolerance"=>.02))
    say("COMPLETE compact both-media qualification; all ",length(BOTH_CASES)," prescribed cases passed")
else
    say("COMPACT FAILED ",failed,"; wider footprint checked here before remaining cases")
    if both_case!(failed...,:decay)
        passed = true
        for case in BOTH_CASES
            case==failed && continue
            if !both_case!(case...,:decay)
                passed = false
                break
            end
        end
        record(joinpath(ROOT,"both-media-selection.toml"),Dict("candidate"=>"localized-decay","passed"=>passed,"relative_tolerance"=>.02))
        say(passed ? "COMPLETE wider both-media qualification" : "BLOCKING wider candidate failed; no production change")
    else
        record(joinpath(ROOT,"both-media-selection.toml"),Dict("candidate"=>"none","passed"=>false,"failed_case"=>string(failed),"relative_tolerance"=>.02))
        say("BLOCKING both candidates failed at ",failed,"; no production change")
    end
end
end

if abspath(PROGRAM_FILE) == @__FILE__
    qualify_both()
end
