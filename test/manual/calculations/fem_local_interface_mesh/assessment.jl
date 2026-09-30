# Read completed qualification records only; no solver or plotting dependencies.
using TOML, Printf
root = normpath(joinpath(@__DIR__,"../../../../.linecablemodels/fem/local-interface-mesh"))
open(joinpath(root,"assessment.csv"),"w") do io
    println(io,"case,passed,maximum_component_change,sign_changes,baseline_nodes,candidate_nodes,baseline_seconds,candidate_seconds,baseline_peak_rss_kib,candidate_peak_rss_kib")
    for dir in sort(readdir(root;join=true))
        comparison = joinpath(dir,"comparison-localized-wave.toml")
        isfile(comparison) || continue
        result = TOML.parsefile(comparison)
        a,b = [TOML.parsefile(joinpath(dir,variant,"mesh.toml")) for variant in ("baseline","localized-wave")]
        sa,sb = [TOML.parsefile(joinpath(dir,variant,"solve.toml")) for variant in ("baseline","localized-wave")]
        println(io,join((basename(dir),result["passed"],result["maximum_component_relative_change"],
            result["new_component_sign_changes"],a["nodes"],b["nodes"],sa["seconds"],sb["seconds"],
            sa["peak_rss_kib"],sb["peak_rss_kib"]),','))
        @printf("%-27s pass=%-5s change=%7.4f%% nodes=%7d -> %7d seconds=%6.2f -> %6.2f\n",
            basename(dir),result["passed"],100result["maximum_component_relative_change"],
            a["nodes"],b["nodes"],sa["seconds"],sb["seconds"])
    end
end
