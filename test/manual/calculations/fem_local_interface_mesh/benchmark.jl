# One warmed repeat of the representative native pair, on the frozen meshes.
# Julia compilation is outside /usr/bin/time's native process measurements.
isdefined(@__MODULE__,:ROOT) || include("qualify.jl")
candidate = isempty(ARGS) ? "localized-wave" : only(ARGS)
folder = candidate=="localized-wave" ? "warmed-repeat" : "warmed-repeat-$candidate"
for variant in ("baseline",candidate)
    source = joinpath(ROOT,"mixed-f0.1-gamma0.0",variant)
    dir = mkpath(joinpath(ROOT,folder,variant))
    if !isfile(joinpath(dir,"sources.toml"))
        for file in [readlines(joinpath(source,".onelab-export-files"));
                ".onelab-export-files";"sources.toml";"study.msh";"mesh.toml"]
            mkpath(dirname(joinpath(dir,file)))
            cp(joinpath(source,file),joinpath(dir,file);force=true)
        end
    end
    digest(joinpath(source,"study.msh"))==digest(joinpath(dir,"study.msh")) || error("Benchmark mesh changed")
    data = solve!(dir)
    timing = last(collect(eachmatch(r"Wall = ([0-9.]+)s, CPU = ([0-9.]+)s",read(joinpath(dir,"solver.log"),String))))
    data["getdp_cpu_seconds"] = parse(Float64,timing[2])
    record(joinpath(dir,"solve.toml"),data)
end
say("COMPLETE warmed native repeat; compare with original pair, not Julia compilation")
