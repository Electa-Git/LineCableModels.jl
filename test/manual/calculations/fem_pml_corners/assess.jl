# Read completed pilot results. This script never invokes Gmsh or GetDP.
using TOML, Printf, Statistics
const ROOT = normpath(joinpath(@__DIR__,"../../../../.linecablemodels/fem/pml-corner-mesh"))
cases = ("air-f0.1-gamma0.0","mixed-f0.1-gamma0.99","mixed-f1.0e6-gamma0.99")
candidates = ("metric-1","metric-1.35")
selection = Dict{String,Any}()
open(joinpath(ROOT,"pilot-assessment.csv"),"w") do io
    println(io,"case,candidate,accuracy_pass,maximum_component_change,new_signs,baseline_corner_triangles,candidate_corner_triangles,baseline_dofs,candidate_dofs,baseline_native_seconds,candidate_native_seconds,candidate_mesh_seconds,candidate_mesh_compile_seconds,baseline_peak_rss_kib,candidate_peak_rss_kib")
    for candidate in candidates
        complete = true; passed = true; ratios = Float64[]
        for case in cases
            dir = joinpath(ROOT,case)
            comparison = joinpath(dir,"comparison-$candidate.toml")
            if !isfile(comparison)
                complete=false
                println("PENDING $case/$candidate")
                continue
            end
            c = TOML.parsefile(comparison)
            a,b = [TOML.parsefile(joinpath(dir,name,"solve.toml")) for name in ("baseline",candidate)]
            ma = TOML.parsefile(joinpath(dir,"baseline/corner-audit.toml"))
            mb = TOML.parsefile(joinpath(dir,candidate,"mesh.toml"))
            passed &= c["passed"]
            push!(ratios,b["seconds"]/a["seconds"])
            println(io,join((case,candidate,c["passed"],c["maximum_component_relative_change"],
                c["new_component_sign_changes"],ma["corner_triangles"],mb["corner_triangles"],
                a["dofs"],b["dofs"],a["seconds"],b["seconds"],mb["mesh_seconds"],mb["compile_seconds"],
                a["peak_rss_kib"],b["peak_rss_kib"]),','))
            @printf("%-27s %-12s pass=%-5s maxchange=%.6g%% signs=%d corners=%d -> %d native=%.2f -> %.2f s mesh=%.3f s\n",
                case,candidate,c["passed"],100c["maximum_component_relative_change"],c["new_component_sign_changes"],
                ma["corner_triangles"],mb["corner_triangles"],a["seconds"],b["seconds"],mb["mesh_seconds"])
        end
        selection[candidate] = Dict("complete"=>complete,"accuracy_pass"=>complete && passed,
            "single_native_run_median_ratio"=>isempty(ratios) ? NaN : median(ratios),
            "timing_qualified"=>false)
    end
end
all_complete = all(v["complete"] for v in values(selection))
survivors = [k for (k,v) in selection if v["accuracy_pass"]]
selection["status"] = !all_complete ? "pending" : isempty(survivors) ? "stop-no-accurate-candidate" : "warmed-timing-required"
selection["production_changes"] = false
selection["relative_tolerance"] = .02
selection["required_warmed_elapsed_reduction"] = .20
open(io->TOML.print(io,selection),joinpath(ROOT,"pilot-selection.toml"),"w")
println("PILOT DECISION: ",selection["status"])
