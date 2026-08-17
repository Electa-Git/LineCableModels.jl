using Gauntlet
using TOML

length(ARGS) in (1, 2) || error(
    "usage: calibrate.jl <normalized-corpus> [evidence.toml]",
)
records = Gauntlet.calibrate_fit(Corpus(ARGS[1]))
payload = Dict(
    "schema_version" => 1,
    "note" => "evidence only; ordinary validation never rewrites tolerances.toml",
    "records" => [Dict(String(key) => value for (key, value) in pairs(record))
     for record in records]
)
if length(ARGS) == 2
    open(ARGS[2], "w") do io
        TOML.print(io, payload; sorted = true)
    end
else
    TOML.print(stdout, payload; sorted = true)
end
