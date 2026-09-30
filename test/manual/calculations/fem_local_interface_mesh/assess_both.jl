# Read completed data only. Keep analytical accuracy separate from FEM preservation.
using TOML, Printf
root = normpath(joinpath(@__DIR__,"../../../../.linecablemodels/fem/local-interface-mesh"))
open(joinpath(root,"both-media-assessment.csv"),"w") do io
    println(io,"case,candidate,quantity,max_change_from_fem,max_baseline_analytical_error,max_candidate_analytical_error,sign_changes,baseline_nodes,candidate_nodes,baseline_seconds,candidate_seconds")
    for dir in sort(readdir(root;join=true)), candidate in ("localized","localized-decay")
        suffix = candidate=="localized" ? "" : "-$candidate"
        marker = joinpath(dir,"two-percent-comparison$suffix.toml")
        isfile(marker) || continue
        result = TOML.parsefile(marker)
        a,b = [TOML.parsefile(joinpath(dir,v,"mesh.toml")) for v in ("baseline",candidate)]
        sa,sb = [TOML.parsefile(joinpath(dir,v,"solve.toml")) for v in ("baseline",candidate)]
        rows = split.(readlines(joinpath(dir,"components$suffix.csv"))[2:end],',')
        @printf("%-27s %-16s passed=%-5s maximum=%.5f%% signs=%d nodes=%d->%d native=%.2f->%.2f s\n",
            basename(dir),candidate,result["passed"],100result["maximum_component_relative_change"],
            result["new_component_sign_changes"],a["nodes"],b["nodes"],sa["seconds"],sb["seconds"])
        for q in ("R","X","G","B")
            selected = filter(r->r[1]==q,rows)
            values = [parse.(Float64,r[4:6]) for r in selected]
            relative(x,y) = iszero(y) ? (iszero(x) ? 0. : Inf) : abs((x-y)/y)
            delta = maximum(relative(b,a) for (a,b,r) in values)
            ea = maximum(relative(a,r) for (a,b,r) in values)
            eb = maximum(relative(b,r) for (a,b,r) in values)
            flips = count(sign(a)!=sign(b) for (a,b,r) in values)
            println(io,join((basename(dir),candidate,q,delta,ea,eb,flips,a["nodes"],b["nodes"],sa["seconds"],sb["seconds"]),','))
        end
    end
end
